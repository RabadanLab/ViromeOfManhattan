#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=8G,time=2::

# This script finds ORFs in the contigs which didn't blast

# exit if previous step produced zero output
if [ ! -s blast/contigs_no_blastn.fa ]; then exit; fi

echo "------------------------------------------------------------------"
echo DISCOVERY START [[ `date` ]]

mkdir -p discovery
echo prodigal commencing [ `date` ]
prodigal -p meta -f gbk -i blast/contigs_no_blastn.fa -o discovery/prodigal_coords.gbk -a discovery/prodigal_proteins.fasta -s discovery/prodigal_scores.txt
echo prodigal finished [ `date` ]

echo DISCOVERY END [[ `date` ]]
echo "------------------------------------------------------------------"
