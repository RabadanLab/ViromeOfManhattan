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

# ids that didn't blast
noblast=""

for f in ${outputdir}/*.result; do
	base=`echo $f | cut -d"." -f1`
	len=`cat $f | wc -l`
	if [ $len -eq 0 ]; then
		noblast=${noblast}","$( basename $base )
		cat ${base}.fasta | sed s/X//g
	fi
done > ${outputdir}/no_blastn.fa

if [ ! -z $noblast ]; then
	echo no blastn hits for $( echo $noblast | sed 's/,//' )
fi

echo "concatenate blast results"
# just get top hit
for i in ${outputdir}/*.result; do head -1 $i; done > ${outputdir}/top.concat.txt
# concat all
cat ${outputdir}/*.result > ${outputdir}/concat.txt
cat ${outputdir}/*.fasta > ${outputdir}/above_threshold.fa

# concat blast logs and remove folder
echo "concatenate blast logs"
head -100 ${logsdir}/* > ${outputdir}/log.blast
rm -rf ${logsdir} 

if [ ${noclean} -eq 0 ]; then
	echo clean up
	rm ${outputdir}/*.result
	rm ${outputdir}/*.fasta
fi
