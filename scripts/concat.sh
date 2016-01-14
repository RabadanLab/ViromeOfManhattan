#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=1G,time=1::

echo "concatenate top blast hits"
for i in blast/*.result; do head -1 $i; done > blast/top.concat

# concatenate blast logs and remove folder
echo "concatenate blast logs"
head logs_blast/* > log.blast
rm -rf logs_blast
