#!/bin/bash
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=12G,time=12::
#$ -pe smp 8

# This script performs assembly on the reads leftover after host removal

# defaults
noclean=0		# no clean boolean
output="assembly/contigs_trinity.fasta"

while [[ $# > 0 ]]; do

	flag=${1}

	case $flag in
		-o|--outputdir)	# the output directory
		outputdir="${2}"
		shift ;;

		-l|--logsdir)	# the logs directory
		logsdir="${2}"
		shift ;;

		-d|--scripts)	# the git repository directory
		d="${2}"
		shift ;;

		--noclean)	# noclean bool
		noclean="${2}"
		shift ;;

		-v|--verbose)	# verbose
		verbose=true ;;

		*)
				# unknown option
		;;
	esac
	shift
done

# function to check if file is zero size or doesn't exist
function iszero {
	if [ ! -s $1 ]; then exit; fi
} 

# exit if previous step produced zero output
iszero host_separation/unmapped_1.fastq.gz
iszero host_separation/unmapped_2.fastq.gz

echo "------------------------------------------------------------------"
echo ASSEMBLY START [[ `date` ]]

mkdir -p assembly_trinity assembly

# perform Trinity assembly
Trinity --seqType fq \
 --max_memory 50G \
 --CPU 8 \
 --normalize_reads \
 --output assembly_trinity \
 --left host_separation/unmapped_1.fastq.gz \
 --right host_separation/unmapped_2.fastq.gz \

# exit if no output
iszero assembly_trinity/Trinity.fasta

# rename Trinity contigs, join sequence portion of fasta 
cat assembly_trinity/Trinity.fasta | awk 'BEGIN{f=0; counter=1}{if ($0~/^>/) {if (f) {printf "\n"; counter++}; print ">contig_"counter; f=1} else printf $0}END{printf "\n"}' > ${output}

# get number of contigs
NUM_CONTIGS=$(( $( cat assembly/contigs_trinity.fasta | wc -l )/2 ))

# compute simple distribution
cat ${output} | paste - - | awk '{print length($2)}' | sort -nr | ${d}/scripts/tablecount | awk -v tot=${NUM_CONTIGS} 'BEGIN{x=0}{x+=$2; print $1"\t"$2"\t"x"/"tot"\t"int(100*x/tot)"%"}' > assembly/contigs.distrib.txt

if [ ${noclean} -eq 0 ]; then
	echo clean up
	rm -r assembly_trinity
fi

echo ASSEMBLY END [[ `date` ]]
echo "------------------------------------------------------------------"

echo "------------------------------------------------------------------"
echo REMAP START [[ `date` ]]

mkdir -p assembly/ref_remap
refbowtie="assembly/ref_remap/ref"

echo Bowtie2 ref build commenced [ `date` ]
bowtie2-build ${output} ${refbowtie}
echo Bowtie2 ref build finished [ `date` ]

echo Bowtie2 mapping commenced [ `date` ]
bowtie2 -p 4 -x ${refbowtie} -1 host_separation/unmapped_1.fastq.gz -2 host_separation/unmapped_2.fastq.gz -S assembly/reads2contigs.sam
echo Bowtie2 mapping finished [ `date` ]

# convert to bam
samtools view -bS assembly/reads2contigs.sam > assembly/reads2contigs.bam
rm assembly/reads2contigs.sam 

if [ ${noclean} -eq 0 ]; then
	echo clean up
	rm -r assembly/ref_remap
fi

echo REMAP END [[ `date` ]]
echo "------------------------------------------------------------------"
