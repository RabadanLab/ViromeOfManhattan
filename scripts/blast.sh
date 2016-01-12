#!/bin/sh
#$ -V
#$ -cwd
#$ -o logs_blast/
#$ -e logs_blast/
#$ -l mem=4G,time=4::

# A simple wrapper for blast array job 

# blast db
blastdb=${1}
# blast format string
fmt=${2}
# blast options
# opts="-task megablast -evalue 0.01"

echo "------------------------------------------------------------------"
echo BLAST ${SGE_TASK_ID} START [[ `date` ]]

blastn -outfmt "6 ${fmt}" -query blast/contig_${SGE_TASK_ID}.fasta -db ${blastdb} > blast/contig_${SGE_TASK_ID}.result;

echo BLAST ${SGE_TASK_ID} END [[ `date` ]]
echo "------------------------------------------------------------------"
