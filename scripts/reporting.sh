#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=1G,time=1::

# This script generates the report

d=${1}		# directory where the parent script resides
id=${2}		# sample identifier
blacklist=${3}	# text file of blacklist taxids

# exit if previous step produced zero output
if [ ! -s blast/top.concat.txt ]; then exit; fi

echo "------------------------------------------------------------------"
echo REPORTING START [[ `date` ]]

mkdir -p report

echo filtering blast results
# filter PREDICTED; sort by taxids then query sequence length (careful: this line can scramble the header)
${d}/scripts/makereport.py blast/header blast/top.concat.txt ${id} ${blacklist} | grep -v PREDICTED | sort -k5,5n -k6,6nr > report/blast.topfilter.txt

echo REPORTING END [[ `date` ]]
echo "------------------------------------------------------------------"
