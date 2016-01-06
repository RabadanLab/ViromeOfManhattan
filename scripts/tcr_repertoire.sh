#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=24G,time=4::

# This script finds T Cell Receptor Motifs

# directory where this script resides
# d=`dirname $( readlink -m $0 )`

mate1=${1}	# mate 1 
mate2=${2}	# mate 2
d=${3}		# directory where the parent script resides

echo "------------------------------------------------------------------"
echo TCR CLONALITY START [[ `date` ]]

mkdir -p tcr

echo zcat [ `date` ] 
zcat ${mate1} ${mate2} > tcr/both_mates.fastq

echo alpha chain [ `date` ] 
java -Xmx16g -jar ${d}/resources/mitcr.jar -pset flex -level 3 -gene TRA tcr/both_mates.fastq tcr/alpha_chain.txt

echo beta chain [ `date` ] 
java -Xmx16g -jar ${d}/resources/mitcr.jar -pset flex -level 3 -gene TRB tcr/both_mates.fastq tcr/beta_chain.txt

echo clean up
rm tcr/both_mates.fastq

echo TCR CLONALITY END [[ `date` ]]
echo "------------------------------------------------------------------"
