#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=24G,time=4::

echo "------------------------------------------------------------------"
echo TCR CLONALITY START [[ `date` ]]

mkdir tcr
echo zcat [ `date` ] 
zcat $1 $2 > tcr/both_mates.fastq
echo alpha chain [ `date` ] 
java -Xmx16g -jar /ifs/home/c2b2/rr_lab/siz2102/apps/Pandora/mitcr.jar -pset flex -level 3 -gene TRA tcr/both_mates.fastq tcr/alpha_chain.txt
echo beta chain [ `date` ] 
java -Xmx16g -jar /ifs/home/c2b2/rr_lab/siz2102/apps/Pandora/mitcr.jar -pset flex -level 3 -gene TRB tcr/both_mates.fastq tcr/beta_chain.txt
rm tcr/both_mates.fastq

echo TCR CLONALITY END [[ `date` ]]
echo "------------------------------------------------------------------"
