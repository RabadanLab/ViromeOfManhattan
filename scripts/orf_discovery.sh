#!/bin/bash
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=2G,time=2::

# This script finds ORFs in the contigs which didn't blast

# defaults
input="blast/no_blastn.fa"
doblast=0	# dont blast by default

while [[ $# > 0 ]]; do

	flag=${1}

	case $flag in
		-i|--input)	# the input fasta
		input="${2}"
		shift ;;

		-o|--outputdir)	# the output directory
		outputdir="${2}"
		shift ;;

		-l|--logsdir)	# the logs directory
		logsdir="${2}"
		shift ;;

		-d|--scripts)	# the git repository directory
		d="${2}"
		shift ;;

		--threshold)	# ORF threshold
		orfthreshold="${2}"
		shift ;;

		--db)		# the database prefix
		dbprefix="${2}"
		shift ;;

		--id)		# id
		id="${2}"
		shift ;;

		--blast)	# blast bool
		doblast="${2}"
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
	if [ ! -s $1 ]; then 
		echo "${1} empty. Exiting..."
		exit; 
	fi
} 

# exit if previous step produced zero output
iszero ${input}

echo "------------------------------------------------------------------"
echo DISCOVERY START [[ `date` ]]

mkdir -p discovery
${d}/scripts/orf.py -i ${input} -t ${orfthreshold} > discovery/orf.fa

# blastp to nr, if blast flag AND orf.fa nonempty
if [ ${doblast} -eq 1 -a -s discovery/orf.fa ]; then
	${d}/scripts/blast_wrapper.sh --scripts ${d} --outputdir discovery/blast -i discovery/orf.fa --logsdir logs_blast2 --whichblast blastp --threshold 100 --db ${dbprefix} --id ${id} --noclean ${noclean}
fi

echo DISCOVERY END [[ `date` ]]
echo "------------------------------------------------------------------"
