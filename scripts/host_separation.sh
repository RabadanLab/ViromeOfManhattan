#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=16G,time=12::
#$ -pe smp 4

# This script removes host (human) reads
# using a pass of STAR mapping, followed by
# a pass of bowtie2 mapping

mate1=${1}	# mate 1 
mate2=${2}	# mate 2
refstar=${3}	# STAR ref
refbowtie=${4}	# bowtie ref
d=${5}		# directory where the parent script resides
gz=${6}		# gzip boolean
noclean=${7}	# no clean boolean

echo "------------------------------------------------------------------"
echo HOST_SEPARATION START [[ `date` ]]

mkdir -p host_separation

# flags for STAR
starflag=""
# if input files are gzipped
if [ ${gz} -eq 1 ]; then
	starflag="--readFilesCommand zcat"
fi

echo STAR mapping commenced [ `date` ]
STAR --runThreadN 4 --genomeDir ${refstar} --readFilesIn ${mate1} ${mate2} --outFileNamePrefix host_separation/ --outSAMtype BAM Unsorted --outSAMunmapped Within ${starflag}
echo STAR mapping finished [ `date`  ]

echo find unmapped reads
samtools flagstat host_separation/Aligned.out.bam > host_separation/mapping_stats.txt
samtools view -b -f 13 host_separation/Aligned.out.bam | samtools sort -n - host_separation/star_unmapped
samtools view host_separation/star_unmapped.bam | ${d}/scripts/sam2fastq.py host_separation/star_unmapped

echo Bowtie2 mapping commenced [ `date` ]
bowtie2 -p 4 -x ${refbowtie} -1 host_separation/star_unmapped_1.fastq -2 host_separation/star_unmapped_2.fastq -S host_separation/bwt2.sam
echo Bowtie2 mapping finished [ `date` ]

echo find unmapped reads
samtools view -S -b -f 13 host_separation/bwt2.sam | samtools sort -n - host_separation/bwt2_unmapped
samtools view host_separation/bwt2_unmapped.bam | ${d}/scripts/sam2fastq.py host_separation/bwt2_unmapped

if [ ${noclean} -eq 0 ]; then
	echo clean up
	rm -rf host_separation/_STARtmp
	rm host_separation/Aligned.out.bam
	rm host_separation/Log.*
	rm host_separation/SJ.out.tab
	rm host_separation/star_unmapped.bam
	rm host_separation/star_unmapped_[12].fastq
	rm host_separation/bwt2.sam
	rm host_separation/bwt2_unmapped.bam
fi

mv host_separation/bwt2_unmapped_1.fastq host_separation/unmapped_1.fastq
mv host_separation/bwt2_unmapped_2.fastq host_separation/unmapped_2.fastq

echo gzip commenced [ `date` ]
gzip host_separation/unmapped_1.fastq
gzip host_separation/unmapped_2.fastq
echo gzip finished [ `date` ]

echo HOST_SEPARATION END [[ `date` ]]
echo "------------------------------------------------------------------"
