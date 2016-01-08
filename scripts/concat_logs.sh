#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=1G,time=1::

# concatenate blast logs and remove folder

head logs_blast/* > log.blast
rm -rf logs_blast
