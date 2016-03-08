#!/bin/bash
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=2G,time=2::

# This script finds ORFs in the contigs which didn't blast

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

		--threshold)	# ORF threshold
		orfthreshold="${2}"
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

# exit if previous step produced zero output
if [ ! -s blast/contigs_no_blastn.fa ]; then exit; fi

echo "------------------------------------------------------------------"
echo DISCOVERY START [[ `date` ]]

mkdir -p discovery
${d}/scripts/orf.py -i blast/contigs_no_blastn.fa -t ${orfthreshold} > discovery/orf.fa

# blastp to nr

echo DISCOVERY END [[ `date` ]]
echo "------------------------------------------------------------------"
