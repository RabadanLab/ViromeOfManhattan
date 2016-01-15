#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=8G,time=2::

echo "------------------------------------------------------------------"
echo DISCOVERY START [[ `date` ]]

mkdir -p discovery

# if file is not zero size
if [ -s blast/contigs_no_blastn.fa ]; then
    echo prodigal commencing [ `date` ]
    prodigal -p meta -f gbk -i blast/contigs_no_blastn.fa -o discovery/prodigal_coords.gbk -a discovery/prodigal_proteins.fasta -s discovery/prodigal_scores.txt
    echo prodigal finished [ `date` ]
else
    echo all contigs have blastn hits
fi

echo DISCOVERY END [[ `date` ]]
echo "------------------------------------------------------------------"

#    echo BLASTX commencing [ `date` ]
#    blastx -query blast/contigs_no_blastn.fa -db nr -out discovery/blastx_result.xml -outfmt 5
#    echo BLASTX finished [ `date` ]

#    echo findorf commencing [ `date` ]
#    findorf join --output discovery/findorf_joined_blastx_dbs.pkl --ref blast/contigs_no_blastn.fa discovery/blastx_result.xml
#    findorf predict -v --gtf discovery/findorf_orfs.gtf --protein discovery/findorf_proteins.fasta --orf discovery/findorf_orfs.fasta --input discovery/findorf_joined_blastx_dbs.pkl
#    echo findorf finished [ `date` ]
