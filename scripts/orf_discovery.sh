#!/bin/sh
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=8G,time=2::

echo "------------------------------------------------------------------"
echo DISCOVERY START [[ `date` ]]

mkdir -p discovery

for f in blast/*.result; do
    base=`echo $f | cut -d"." -f1`
    len=`cat $f | wc -l`
    if [ $len -lt 43 ]; then 
        echo no blastn hits for `echo $base | cut -d"/" -f2`
        cat $base.fasta | sed s/X//g >> discovery/contigs_no_blastn.fasta
    fi
done

if [ -e discovery/contigs_no_blastn.fasta ]; then

    echo prodigal commencing [ `date` ]
    prodigal -p meta -f gbk -i discovery/contigs_no_blastn.fasta -o discovery/prodigal_coords.gbk -a discovery/prodigal_proteins.fasta -s discovery/prodigal_scores.txt
    echo prodigal finished [ `date` ]
else
    echo all contigs have blastn hits
fi

echo DISCOVERY END [[ `date` ]]
echo "------------------------------------------------------------------"

#    echo BLASTX commencing [ `date` ]
#    blastx -query discovery/contigs_no_blastn.fasta -db /ifs/scratch/c2b2/rr_lab/shares/ref/blastdb/nr-2012-04-08/nr -out discovery/blastx_result.xml -outfmt 5
#    echo BLASTX finished [ `date` ]

#    echo findorf commencing [ `date` ]
#    findorf join --output discovery/findorf_joined_blastx_dbs.pkl --ref discovery/contigs_no_blastn.fasta discovery/blastx_result.xml
#    findorf predict -v --gtf discovery/findorf_orfs.gtf --protein discovery/findorf_proteins.fasta --orf discovery/findorf_orfs.fasta --input discovery/findorf_joined_blastx_dbs.pkl
#    echo findorf finished [ `date` ]
