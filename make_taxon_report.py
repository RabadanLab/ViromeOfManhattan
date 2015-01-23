#!/usr/bin/env python

import subprocess
from Bio import SeqIO, Entrez
import pandas as pd

Entrez.email = "Your.Name.Here@example.org"

top_blast = pd.read_csv('report/top_blast_hits.tsv', index_col=0, header=0, sep='\t')
top_blast['taxon'] = top_blast['accession']

for contig in top_blast.index:
    acc_id = top_blast.ix[contig, 'accession']
    if not pd.isnull(acc_id):
        print acc_id
        handle = Entrez.efetch(db="nucleotide", id=acc_id, rettype="gb", retmode="text")
        text = handle.read()
        f_gbk = open('report/tmp.gbk','w')
        f_gbk.write(text)
        f_gbk.close()

        tax_id = subprocess.check_output('grep taxon report/tmp.gbk | cut -d: -f2 | cut -d\\" -f1', shell=True).strip()
        top_blast.ix[contig, 'taxon'] = tax_id

top_blast.to_csv('report/with_taxon_id.tsv', sep='\t')
