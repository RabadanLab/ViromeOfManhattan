#!/bin/bash

# simple wrapper for blast

echo [start]
echo [pwd] `pwd`
echo [date] `date`

db=$1

input=${SGE_TASK_ID}.fasta
output=${SGE_TASK_ID}.result

echo "input "${input}
echo "output "${output}
echo "db "${db}

# restrict to top 10 hits
blastn -num_alignments 10 -outfmt "6 qseqid sseqid saccver staxids pident nident length mismatch gapopen gaps qstart qend qlen qframe qcovs sstart send slen sframe sstrand evalue bitscore stitle" -query ${input} -db ${db} > ${output}

echo [finish]
echo [date] `date`
