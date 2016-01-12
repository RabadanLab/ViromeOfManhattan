#!/bin/sh

echo "------------------------------------------------------------------"
echo REPORTING START [[ `date` ]]

mkdir -p report
echo contig_id"	"contig_len"	"accession"	"description > report/top_blast_hits.tsv;

for f in blast/*.result; do 
	ctg_num=`echo $f | cut -d"/" -f2 | cut -d"." -f1 | cut -d"_" -f2`;
	ctg_len=`tail -n +2 blast/contig_$ctg_num.fasta | wc -c`;
	top_hit=`sed -n 21,21p $f`;
	top_acc=`echo $top_hit | cut -d"|" -f2`;
	top_rest=`echo $top_hit | cut -d"|" -f3`;
	echo $1_$ctg_num"	"$ctg_len"	"$top_acc"	"$top_rest >> report/top_blast_hits.tsv;
done

echo REPORTING END [[ `date` ]]
echo "------------------------------------------------------------------"
