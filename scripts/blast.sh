#!/bin/sh
#$ -V
#$ -cwd
#$ -l mem=4G,time=4::

# A simple wrapper for blast array job 

# defaults
outputdir="blast"

while [[ $# > 0 ]]; do

	flag=${1}

	case $flag in
		-o|--outputdir)	# the output directory
		outputdir="${2}"
		shift ;;

		-d|--scripts)	# the git repository directory
		d="${2}"
		shift ;;

		--whichblast)	# which blast (blastn or blastp)
		wblast="${2}"
		shift ;;

		--db)		# blast db
		blastdb="${2}"
		shift ;;

		--fmt)		# blast format string
		fmt="${2}"
		shift ;;

		--sgeid)	# set SGE_TASK_ID by hand (for the case where qsub is turned off)
		SGE_TASK_ID="${2}"
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

# set input
input="${outputdir}/blast_${SGE_TASK_ID}.fasta"

echo "------------------------------------------------------------------"
echo BLAST ${SGE_TASK_ID} START [[ `date` ]]

# tmp hack
if [ "${wblast}" == "blastp" ]; then
	wblast=${wblast}" -task blastp-fast"
fi

${wblast} -outfmt "6 ${fmt}" -query ${input} -db ${blastdb} > ${outputdir}/blast_${SGE_TASK_ID}.result;

echo BLAST ${SGE_TASK_ID} END [[ `date` ]]
echo "------------------------------------------------------------------"
