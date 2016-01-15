#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=14G,time=2::
#$ -pe smp 4

# This script generates the report

blacklist=${1}	# blacklist Python pickle file
d=${2}		# directory where the parent script resides
remap=${3}	# boolean: create pathogen fasta and map reads back onto this reference
blastdb=${4}	# blast db
noclean=${5}	# no clean boolean
mate1="host_separation/unmapped_1.fastq.gz"	# mate 1 
mate2="host_separation/unmapped_2.fastq.gz"	# mate 2

echo "------------------------------------------------------------------"
echo REPORTING START [[ `date` ]]

mkdir -p report

echo filtering blast results
cat blast/header > report/blast.topfilter.txt
${d}/scripts/makereport.py $blacklist blast/top.concat.txt >> report/blast.topfilter.txt

if [ ${remap} -eq 1 ]; then
	echo creating a pathogen reference
	# get seqids of hits
	cat report/blast.topfilter.txt | sed '1d' | cut -f2 | sort -u > report/entry.batch.txt
	# create a fasta file of these sequences
	blastdbcmd -db ${blastdb} -entry_batch report/entry.batch.txt -outfmt '%f' > report/pathogen.ref.fa
	echo runnning STAR genome generate
	mkdir -p report/star_ref report/star_map
	# "[seg fault] is a known problem for very small genomes. At the genome generation step, please try to reduce 
	# the value of --genomeSAindexNbases to <=8, and then re-run the mapping step. Generally, --genomeSAindexNbases 
	# needs to be scaled with the genome length, as ~min(14,log2(ReferenceLength)/2 - 1). 
	# I will need to incorporate this scaling in the future releases." -alexdobin
	STAR --runMode genomeGenerate --genomeDir report/star_ref --genomeFastaFiles report/pathogen.ref.fa --genomeSAindexNbases 8
	echo STAR mapping commenced [ `date` ]
	STAR --readFilesCommand zcat --runThreadN 4 --genomeDir report/star_ref --readFilesIn ${mate1} ${mate2} --outFileNamePrefix report/star_map/ --outSAMtype BAM Unsorted --outSAMunmapped Within
	echo STAR mapping finished [ `date`  ]
fi

if [ ${noclean} -eq 0 ]; then
	echo clean up
	rm report/entry.batch.txt
	rm report/pathogen.ref.fa
	rm -r report/star_ref
fi

echo REPORTING END [[ `date` ]]
echo "------------------------------------------------------------------"
