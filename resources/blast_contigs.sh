#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=4G,time=160::

# This script blasts the assembled contigs, 
# provided they're over the length threshold

contigthreshold=${1}	# contig threshold
blastdb=${2}		# blast db

echo "------------------------------------------------------------------"
echo BLAST START [[ `date` ]]

mkdir -p blast

NUM_CONTIGS=$(( `wc -l assembly/contigs_trinity.fasta | cut -f1 -d" "` / 2 ))

for i in $(seq 1 $NUM_CONTIGS); do 
    line1=$(( 2 * i - 1 ));
    line2=$(( 2 * i ));
    NUM_CHARS=`sed -n "$line2,$line2"p assembly/contigs_trinity.fasta | wc -c`;
    if [ $NUM_CHARS -gt ${contigthreshold} ]; then
        echo BLAST contig $i of $NUM_CONTIGS, length is $NUM_CHARS;
        sed -n "$line1,$line2"p assembly/contigs_trinity.fasta > blast/contig_$i.fasta;
        blastn -query blast/contig_$i.fasta -db ${blastdb} > blast/contig_$i.result;
    else 
        echo contig $i of $NUM_CONTIGS too short, length is $NUM_CHARS;
        fi;
    done

echo BLAST END [[ `date` ]]
echo "------------------------------------------------------------------"
