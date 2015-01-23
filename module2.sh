#!/bin/sh
#$ -V
#$ -cwd
#$ -o out.sge
#$ -e err.sge
#$ -l mem=24G,time=96::
#$ -pe smp 2

echo "------------------------------------------------------------------"
echo "------------------------------------------------------------------"
echo ASSEMBLY START [[ `date` ]]

mkdir assembly

Trinity --seqType fq --JM 20G --output assembly --left host_separation/bwt2_unmapped_1.fastq  --right host_separation/bwt2_unmapped_2.fastq

echo ASSEMBLY END [[ `date` ]]
