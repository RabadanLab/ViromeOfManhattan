#!/bin/sh
#$ -V
#$ -cwd
#$ -o logs_blast/
#$ -e logs_blast/
#$ -l mem=4G,time=4::

# A simple wrapper for blast array job 

blastdb=${1}	# blast db

echo "------------------------------------------------------------------"
echo BLAST ${SGE_TASK_ID} START [[ `date` ]]

blastn -query blast/contig_${SGE_TASK_ID}.fasta -db ${blastdb} > blast/contig_${SGE_TASK_ID}.result;

echo BLAST ${SGE_TASK_ID} END [[ `date` ]]
echo "------------------------------------------------------------------"
