#!/bin/bash

bamfile=$1
outputdir=$2
scripts=$3
blastdb=$4
blacklist=$5
taxid2names=$6
header=$7

# output directory
mkdir -p ${outputdir}
if ! cd ${outputdir}; then exit; fi
mkdir -p logs

# get unmapped reads
samtools view -b -f 4 ../${bamfile} > unassembled.bam

# transform to fasta and split
samtools view unassembled.bam | cut -f1,10 > tmp
split -d -l 1000 -a 3 tmp
num=0
for i in x*; do ((num++)); cat $i | tr ":" "_" | awk '{print ">"$1; print $2;}' > ${num}.fasta; done
echo "num: "${num}
# filecount=$( ls -1 *fasta | wc -l )

# blast
# message looks like this: Your job-array 3147839.1-2:1 ("wait") has been submitted
jid=$( qsub -V -N bunassembled -e logs -o logs -cwd -l mem=10G,time=6:: -t 1-${num} ${scripts}/scripts/blast.sh ${blastdb} | cut -f3 -d' ' | cut -f1 -d'.' )
echo "job id: "${jid}
# pause
qsub -V -b y -o logs -e logs -cwd -N wait -hold_jid ${jid} -sync y echo wait_here

# get top hit, concatenate
for i in *.result; do ${scripts}/scripts/filterblast.py ${i} ${i}.tophit ${scripts}; done
cat *.result > all.concat.txt
cat *.tophit > top.concat.txt
head -1000 logs/* > all.logs.txt
# clean up
echo 'clean up'
rm x*
rm tmp
rm *.fasta
rm *.result
rm *.tophit
rm -r logs
# generate report
${scripts}/scripts/makereport.py --header ${header} --scripts ${scripts} -i top.concat.txt --outputdir . --id 1 --verbose 1 --blacklist ${blacklist} --taxid2names ${taxid2names} --hpc 1
