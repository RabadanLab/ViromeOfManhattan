#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=4G,time=8::

# This script blasts the entries of a fasta file,
# provided they're over the length threshold

# defaults
wblast="blastn"		# which blast program to use
input="assembly/contigs_trinity.fasta"
outputdir="blast"	# directories
logsdir="logs_blast"
lenthreshold=0		# length threshold
seqtype="contig"	# kind of seq data (e.g., contig, ORF)
noclean=0		# no clean boolean
noSGE=0			# sge boolean
# blast format string
fmt="qseqid sseqid saccver staxids pident nident length mismatch gapopen gaps qstart qend qlen qframe qcovs sstart send slen sframe sstrand evalue bitscore stitle"

# argparse derived from http://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
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
		-d|--repodir)	# the git repository directory
		d="${2}"
		shift ;;
		--threshold)	# the length threshold
		lenthreshold="${2}"
		shift ;;
		--db)		# the database prefix
		dbprefix="${2}"
		shift ;;
		-b|--blast)	# which blast to use (blastn, blastp, etc)
		wblast="${2}"
		shift ;;
		--seqtype)	# an id string for sequence type
		seqtype="${2}"
		shift ;;
		--nosge)	# no SGE bool
		noSGE="${2}"
		shift ;;
		--noclean)	# noclean bool
		noclean="${2}"
		shift ;;
		--id)		# id
		id="${2}"
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
if [ ! -s ${input} ]; then exit; fi

echo "------------------------------------------------------------------"
echo BLAST START [[ `date` ]]

mkdir -p ${outputdir} ${logsdir}

# return counter for contigs above threshold
# assume fastajoinlines, i.e., only one sequence line per entry
j=$( cat ${input} | paste - - | awk -v cutoff=${lenthreshold} -v dir=${outputdir} -v seq=${seqtype} 'BEGIN{counter=0}{
	if (length($2) >= cutoff) {
		counter++;
		myfile=dir"/"seq"_"counter".fasta";
		print ">"seq"_"counter" (formerly_"substr($1,2)")" > myfile;
		print $2 >> myfile
	}
}END{print counter}' )

echo ${fmt} | sed 's/ /\t/g' > ${outputdir}/header

# if qsub
if [ ${noSGE} -eq 0 ]; then
	message=$( qsub -N bc_${id} -e ${logsdir} -o ${logsdir} -t 1-${j} ${d}/scripts/blast.sh ${wblast} ${dbprefix} "${fmt}" ${seqtype} | grep submitted )
	echo $message
	jid=$( echo $message | cut -f3 -d' ' | cut -f1 -d'.' )
	# message should be like: 'Your job-array 8388982.1-256:1 ("bc_5") has been submitted'
	# hold the script up here, until all the blast jobs finish
	# concat top blast hits; concat log files into one, so as not to clutter the file system
	qsub -N wait_${id} -hold_jid ${jid} -sync y ${d}/scripts/concat.sh ${outputdir} ${logsdir} ${noclean} ${seqtype}
# if no qsub
else

	for i in $( seq ${j} ); do
		${d}/scripts/blast.sh ${wblast} ${dbprefix} "${fmt}" ${seqtype} ${i} > ${logsdir}/bc_${id}.${i}.o 2> ${logsdir}/bc_${id}.${i}.e
	done

	${d}/scripts/concat.sh ${outputdir} ${logsdir} ${noclean} ${seqtype}
fi

echo BLAST END [[ `date` ]]
echo "------------------------------------------------------------------"
