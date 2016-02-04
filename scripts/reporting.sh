#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=1G,time=1::

# This script generates the report

d=${1}		# directory where the parent script resides
blacklist=${2}	# text file of blacklist taxids

echo "------------------------------------------------------------------"
echo REPORTING START [[ `date` ]]

# if blast result doesn't exist
if [ ! -e blast/top.concat.txt ]; then
	echo ERROR blast file not found
	exit
fi

mkdir -p report

echo filtering blast results
${d}/scripts/makereport.py blast/header blast/top.concat.txt $blacklist | grep -v PREDICTED | sort -k 4,4n -k 5,5nr > report/blast.topfilter.txt

echo REPORTING END [[ `date` ]]
echo "------------------------------------------------------------------"
