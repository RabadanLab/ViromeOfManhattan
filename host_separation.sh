#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=16G,time=4::
#$ -pe smp 4

echo "------------------------------------------------------------------"
echo HOST_SEPARATION START [[ `date` ]]

mkdir host_separation
echo STAR mapping commenced [ `date` ]
STAR --readFilesCommand zcat --runThreadN 4 --genomeDir /ifs/scratch/c2b2/rr_lab/siz2102/reference/human/star --readFilesIn $1 $2 --outFileNamePrefix host_separation/ --outSAMtype BAM Unsorted --outSAMunmapped Within
echo STAR mapping finished [ `date`  ]
samtools flagstat host_separation/Aligned.out.bam > host_separation/mapping_stats.txt
samtools view -b -f 13 host_separation/Aligned.out.bam | samtools sort -n - host_separation/star_unmapped
bam bam2FastQ --in host_separation/star_unmapped.bam --readname --outBase host_separation/star_unmapped
echo Bowtie2 mapping commenced [ `date` ]
bowtie2 -p 4 -x /ifs/scratch/c2b2/rr_lab/siz2102/reference/human/bowtie2/hg19 -1 host_separation/star_unmapped_1.fastq -2 host_separation/star_unmapped_2.fastq -S host_separation/bwt2.sam
echo Bowtie2 mapping finished [ `date` ]
samtools view -S -b -f 13 host_separation/bwt2.sam | samtools sort -n - host_separation/bwt2_unmapped
bam bam2FastQ --in host_separation/bwt2_unmapped.bam --readname --outBase host_separation/bwt2_unmapped

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
mv host_separation/bwt2_unmapped_1.fastq host_separation/unmapped_1.fastq
mv host_separation/bwt2_unmapped_2.fastq host_separation/unmapped_2.fastq

echo gzip commenced [ `date` ]
gzip host_separation/unmapped_1.fastq
gzip host_separation/unmapped_2.fastq
echo gzip finished [ `date` ]

echo HOST_SEPARATION END [[ `date` ]]
echo "------------------------------------------------------------------"
