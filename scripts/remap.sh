#!/bin/bash
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=4G,time=4::
#$ -pe smp 4

# This script maps reads back onto the contigs

noclean=${1}	# no clean boolean

echo "------------------------------------------------------------------"
echo REMAP START [[ `date` ]]

contigs="blast/contigs_above_threshold.fa"

# if fasta file of contigs doesn't exist
if [ ! -e ${contigs} ]; then
	echo ERROR contigs file not found
	exit
fi

mkdir -p remap/ref

refbowtie="remap/ref/contigref"

echo Bowtie2 ref build commenced [ `date` ]
bowtie2-build ${contigs} ${refbowtie}
echo Bowtie2 ref build finished [ `date` ]

echo Bowtie2 mapping commenced [ `date` ]
bowtie2 -p 4 -x ${refbowtie} -1 host_separation/unmapped_1.fastq.gz -2 host_separation/unmapped_2.fastq.gz -S remap/unmappedreads2contigs.sam
echo Bowtie2 mapping finished [ `date` ]

# convert to bam
samtools view -bS remap/unmappedreads2contigs.sam > remap/unmappedreads2contigs.bam
rm remap/unmappedreads2contigs.sam

if [ ${noclean} -eq 0 ]; then
	echo clean up
	rm -r remap/ref 
fi

echo REMAP END [[ `date` ]]
echo "------------------------------------------------------------------"

# blastdb=${2}	# blast db
# noclean=${3}	# no clean boolean
#
# mate1="host_separation/unmapped_1.fastq.gz"	# mate 1
# mate2="host_separation/unmapped_2.fastq.gz"	# mate 2
#
# if [ ${remap} -eq 1 ]; then
# 	echo creating a pathogen reference
# 	# get seqids of hits
# 	cat report/blast.topfilter.txt | sed '1d' | cut -f2 | sort -u > report/entry.batch.txt
# 	# create a fasta file of these sequences
# 	blastdbcmd -db ${blastdb} -entry_batch report/entry.batch.txt -outfmt '%f' > report/pathogen.ref.fa
# 	echo runnning STAR genome generate
# 	mkdir -p report/star_ref report/star_map
# 	# "[seg fault] is a known problem for very small genomes. At the genome generation step, please try to reduce
# 	# the value of --genomeSAindexNbases to <=8, and then re-run the mapping step. Generally, --genomeSAindexNbases
# 	# needs to be scaled with the genome length, as ~min(14,log2(ReferenceLength)/2 - 1).
# 	# I will need to incorporate this scaling in the future releases." -alexdobin
# 	STAR --runMode genomeGenerate --genomeDir report/star_ref --genomeFastaFiles report/pathogen.ref.fa --genomeSAindexNbases 8
# 	echo STAR mapping commenced [ `date` ]
# 	STAR --readFilesCommand zcat --runThreadN 4 --genomeDir report/star_ref --readFilesIn ${mate1} ${mate2} --outFileNamePrefix report/star_map/ --outSAMtype BAM Unsorted --outSAMunmapped Within
# 	echo STAR mapping finished [ `date`  ]
# fi
#
# if [ ${noclean} -eq 0 ]; then
# 	echo clean up
# 	rm -f report/entry.batch.txt
# 	rm -f report/pathogen.ref.fa
# 	rm -rf report/star_ref
# fi
