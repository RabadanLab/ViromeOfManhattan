#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=2G,time=2::

# This script finds ORFs in the contigs which didn't blast

d=${1}			# directory where the parent script resides
orfthreshold=${2}	# ORF threshold

# exit if previous step produced zero output
if [ ! -s blast/contigs_no_blastn.fa ]; then exit; fi

echo "------------------------------------------------------------------"
echo DISCOVERY START [[ `date` ]]

mkdir -p discovery
${d}/scripts/orf.py -i blast/contigs_no_blastn.fa -t ${orfthreshold} > discovery/orf.fa

echo DISCOVERY END [[ `date` ]]
echo "------------------------------------------------------------------"
