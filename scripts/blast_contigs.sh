#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=4G,time=8::

# This script blasts the assembled contigs, 
# provided they're over the length threshold

contigthreshold=${1}	# contig threshold
blastdb=${2}		# blast db
id=${3}			# identifier
d=${4}			# directory where the parent script resides

echo "------------------------------------------------------------------"
echo BLAST START [[ `date` ]]

mkdir -p blast logs_blast

# return counter for contigs above threshold
# assume fastajoinlines, i.e., only one sequence line per entry
j=$( cat assembly/contigs_trinity.fasta | paste - - | awk -v cutoff=${contigthreshold} 'BEGIN{counter=0}{
	if (length($2) >= cutoff) {
		counter++;
		myfile="blast/contig_"counter".fasta";
		print ">contig_"counter" (formerly "substr($1,2)")" > myfile;
		print $2 >> myfile
	}
}END{print counter}' )

# blast format string
fmt="qseqid sseqid saccver staxids pident nident length mismatch gapopen gaps qstart qend qlen qframe qcovs sstart send slen sframe sstrand evalue bitscore stitle"

echo ${fmt} | sed 's/ /\t/g' > blast/header

jid=$( qsub -N bc_${id} -t 1-${j} ${d}/scripts/blast.sh ${blastdb} "${fmt}" | cut -f3 -d' ' | cut -f1 -d'.' )
# message should be like: 'Your job-array 8388982.1-256:1 ("bc_5") has been submitted'
# hold the script up here, until all the blast jobs finish
# concat top blast hits; concat log files into one, so as not to clutter the file system
qsub -N wait_${id} -hold_jid ${jid} -sync y ${d}/scripts/concat.sh

echo BLAST END [[ `date` ]]
echo "------------------------------------------------------------------"
