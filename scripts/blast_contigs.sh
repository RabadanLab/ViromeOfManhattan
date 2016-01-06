#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=4G,time=8::

# This script blasts the assembled contigs, 
# provided they're over the length threshold

contigthreshold=${1}	# contig threshold
blastdb=${2}		# blast db
id=${3}			# identifier
d=${4}			# directory where the parent script resides

echo "------------------------------------------------------------------"
echo BLAST START [[ `date` ]]

mkdir -p blast logs_blast

# get number of contigs
NUM_CONTIGS=$(( `wc -l assembly/contigs_trinity.fasta | cut -f1 -d" "` / 2 ))

# counter for contigs above threshold
j=1

for i in $(seq 1 $NUM_CONTIGS); do 

    line1=$(( 2 * i - 1 ));
    line2=$(( 2 * i ));

    NUM_CHARS=`sed -n "$line2,$line2"p assembly/contigs_trinity.fasta | wc -c`;

    if [ $NUM_CHARS -gt ${contigthreshold} ]; then
        echo SAVE contig $i of $NUM_CONTIGS as ${j}, length is $NUM_CHARS;
        sed -n "$line1,$line2"p assembly/contigs_trinity.fasta > blast/contig_${j}.fasta;
        # increment counter
        j=$((${j}+1))
    else 
        echo contig $i of $NUM_CONTIGS too short, length is $NUM_CHARS;
    fi;
done

jid=$( qsub -N bc_${id} -t 1-${j} ${d}/scripts/blast.sh ${blastdb} | cut -f3 -d' ' | cut -f1 -d'.' )
# message should be like: 'Your job-array 8388982.1-256:1 ("Pan_bc_5") has been submitted'
# hold the script up here, until all the blast jobs finish
qsub -cwd -b y -N wait_${id} -e logs_blast -o logs_blast -hold_jid ${jid} -sync y echo hold

echo BLAST END [[ `date` ]]
echo "------------------------------------------------------------------"
