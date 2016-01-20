#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=1G,time=1::

# This script concatenates blast files, logs, etc so as 
# not to leave many files messily scattered on the filesystem

noclean=${1}	# no clean boolean

# contig ids that didn't blast
noblast=""

for f in blast/*.result; do
	base=`echo $f | cut -d"." -f1`
	len=`cat $f | wc -l`
	if [ $len -eq 0 ]; then
		noblast=${noblast}","$( basename $base )
		cat ${base}.fasta | sed s/X//g
	fi
done > blast/contigs_no_blastn.fa

if [ ! -z $noblast ]; then
	echo no blastn hits for $( echo $noblast | sed 's/,//' )
fi

echo "concatenate blast results"
for i in blast/*.result; do head -1 $i; done > blast/top.concat.txt
cat blast/*.result | gzip > blast/concat.txt.gz
cat blast/*.fasta | gzip > blast/all_contigs.fa.gz

# concatenate blast logs and remove folder
echo "concatenate blast logs"
head logs_blast/* > log.blast
rm -rf logs_blast

if [ ${noclean} -eq 0 ]; then
	echo clean up
	rm blast/*.result
	rm blast/*.fasta
fi
