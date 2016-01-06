#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=16G,time=4::
#$ -pe smp 4

# This script removes host (human) reads
# using a pass of STAR mapping, followed by
# a pass of bowtie2 mapping

mate1=${1}	# mate 1 
mate2=${2}	# mate 2
refstar=${3}	# STAR ref
refbowtie=${4}	# bowtie ref
noclean=${5}	# no clean boolean

echo "------------------------------------------------------------------"
echo HOST_SEPARATION START [[ `date` ]]

mkdir -p host_separation

echo STAR mapping commenced [ `date` ]
STAR --readFilesCommand zcat --runThreadN 4 --genomeDir ${refstar} --readFilesIn ${mate1} ${mate2} --outFileNamePrefix host_separation/ --outSAMtype BAM Unsorted --outSAMunmapped Within
echo STAR mapping finished [ `date`  ]

echo find unmapped reads
samtools flagstat host_separation/Aligned.out.bam > host_separation/mapping_stats.txt
samtools view -b -f 13 host_separation/Aligned.out.bam | samtools sort -n - host_separation/star_unmapped
bam bam2FastQ --in host_separation/star_unmapped.bam --readname --outBase host_separation/star_unmapped

echo Bowtie2 mapping commenced [ `date` ]
bowtie2 -p 4 -x ${refbowtie} -1 host_separation/star_unmapped_1.fastq -2 host_separation/star_unmapped_2.fastq -S host_separation/bwt2.sam
echo Bowtie2 mapping finished [ `date` ]

echo find unmapped reads
samtools view -S -b -f 13 host_separation/bwt2.sam | samtools sort -n - host_separation/bwt2_unmapped
bam bam2FastQ --in host_separation/bwt2_unmapped.bam --readname --outBase host_separation/bwt2_unmapped

echo clean up
if [ ${noclean} -eq 0 ]; then
	rm -r host_separation/_STARtmp
	rm host_separation/Aligned.out.bam
	rm host_separation/Log.*
	rm host_separation/SJ.out.tab
	rm host_separation/star_unmapped.bam
	rm host_separation/star_unmapped_[12].fastq
	rm host_separation/star_unmapped.fastq
	rm host_separation/bwt2.sam
	rm host_separation/bwt2_unmapped.bam
	rm host_separation/bwt2_unmapped.fastq
fi

mv host_separation/bwt2_unmapped_1.fastq host_separation/unmapped_1.fastq
mv host_separation/bwt2_unmapped_2.fastq host_separation/unmapped_2.fastq

echo gzip commenced [ `date` ]
gzip host_separation/unmapped_1.fastq
gzip host_separation/unmapped_2.fastq
echo gzip finished [ `date` ]

echo HOST_SEPARATION END [[ `date` ]]
echo "------------------------------------------------------------------"
