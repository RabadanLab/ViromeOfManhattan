#!/bin/sh
#$ -V
#$ -cwd
#$ -e err.sge
#$ -o out.sge
#$ -l mem=4G,time=24::

echo "------------------------------------------------------------------"
echo "------------------------------------------------------------------"
echo BLAST START [[ `date` ]]

mkdir blast

NUM_CONTIGS=$(( `wc -l assembly/chrysalis/bundled_iworm_contigs.fasta | cut -f1 -d" "` / 2 ))

for i in $(seq 1 $NUM_CONTIGS); do 
    line1=$(( 2 * i - 1 ));
    line2=$(( 2 * i ));
    NUM_CHARS=`sed -n "$line1,$line2"p assembly/chrysalis/bundled_iworm_contigs.fasta | wc -c`;
    if [ $NUM_CHARS -gt 1000 ]; then
        echo BLAST contig $i of $NUM_CONTIGS, length is $NUM_CHARS;
        sed -n "$line1,$line2"p assembly/chrysalis/bundled_iworm_contigs.fasta > blast/contig_$i.fasta;
        blastn -query blast/contig_$i.fasta -db nt > blast/contig_$i.result;
    else 
        echo contig $i of $NUM_CONTIGS too short, length is $NUM_CHARS;
        fi;
    done

echo BLAST END [[ `date` ]]
