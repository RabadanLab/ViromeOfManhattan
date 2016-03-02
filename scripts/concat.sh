#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=1G,time=1::

# This script concatenates blast files, logs, etc so as 
# not to leave many files messily scattered on the filesystem

outputdir=${1}	# output directory
logsdir=${2}	# logs directory
noclean=${3}	# no clean boolean
seqtype=${4}	# sequence type (id)

# ids that didn't blast
noblast=""

for f in ${outputdir}/*.result; do
	base=`echo $f | cut -d"." -f1`
	len=`cat $f | wc -l`
	if [ $len -eq 0 ]; then
		noblast=${noblast}","$( basename $base )
		cat ${base}.fasta | sed s/X//g
	fi
done > ${outputdir}/${seqtype}s_no_blastn.fa

if [ ! -z $noblast ]; then
	echo no blastn hits for $( echo $noblast | sed 's/,//' )
fi

echo "concatenate blast results"
for i in ${outputdir}/*.result; do head -1 $i; done > ${outputdir}/top.concat.txt
# cat ${outputdir}/*.result | gzip > ${outputdir}/concat.txt.gz
# cat ${outputdir}/*.fasta | gzip > ${outputdir}/all_${seqtype}s.fa.gz
cat ${outputdir}/*.result > ${outputdir}/concat.txt
cat ${outputdir}/*.fasta > ${outputdir}/${seqtype}s_above_threshold.fa

# concatenate blast logs and remove folder
echo "concatenate blast logs"
head ${logsdir}/* > log.blast
rm -rf ${logsdir} 

if [ ${noclean} -eq 0 ]; then
	echo clean up
	rm ${outputdir}/*.result
	rm ${outputdir}/*.fasta
fi
