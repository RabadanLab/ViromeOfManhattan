#!/bin/sh
#$ -V
#$ -cwd
#$ -l mem=4G,time=4::

# A simple wrapper for blast array job 

wblast=${1}	# which blast (blastn or blastp)
blastdb=${2}	# blast db
fmt=${3}	# blast format string
seqtype=${4}	# seqtype id
# opts="-task megablast -evalue 0.01"	# blast options

# if a fifth, argument, set SGE_TASK_ID by hand (for the case where qsub is turned off)
if [ $# -eq 5 ]; then
	SGE_TASK_ID=${5}
fi

echo "------------------------------------------------------------------"
echo BLAST ${SGE_TASK_ID} START [[ `date` ]]

${wblast} -outfmt "6 ${fmt}" -query blast/${seqtype}_${SGE_TASK_ID}.fasta -db ${blastdb} > blast/${seqtype}_${SGE_TASK_ID}.result;

echo BLAST ${SGE_TASK_ID} END [[ `date` ]]
echo "------------------------------------------------------------------"
