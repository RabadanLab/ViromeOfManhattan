#!/usr/bin/env python

import argparse

# The goal of this script is to filter blast results based on taxid
# Don't want human sequences or anything in the taxid blacklist,
# which contains non-pathogens

# defaults
# taxid corresponding to homo sapiens
humantaxid = '9606'
# desired header
desiredfields = ['qseqid', 'sseqid', 'qlen','saccver','staxids','evalue', 'bitscore', 'stitle']

# parse arguments
prog_description = 'Make tsv report'
parser = argparse.ArgumentParser(description=prog_description)
parser.add_argument('-i', '--input', required=True, help='the contig blast file input')
parser.add_argument('-i2', '--input2', help='the ORF blast file input')
parser.add_argument('--header', required=True, help='the blast header file')
parser.add_argument('--sample', required=True, help='sample identifier')
parser.add_argument('--id2reads', required=True, help='the output file of samtools idxstats mapping ids to #reads')
parser.add_argument('--blacklist', help='the file of blacklist taxids')
args = parser.parse_args()

# taxid blacklist
filterlist = []
# header
header = []
# id 2 reads dict (output of samtools idxstats)
idx = {}

# load blacklist if supplied
if args.blacklist:
    with open(args.blacklist, 'r') as f:
        filterlist = f.read().split('\n')[:-1]	# final newline causes empty list elt

# load idx file
with open(args.id2reads, 'r') as f:
    for line in f:
        # map id to #reads
        idx[line.split()[0].strip()] = line.split()[2].strip()

# load header
with open(args.header, 'r') as f:
    header = map(str.strip, f.read().split())

# get indicies,fields in file
myindicies = []
myfields = []
for j,k in enumerate(header):
    if k in desiredfields:
        myindicies.append(j)
        myfields.append(k)

# print header:
print('sampleid\t' + '\t'.join(myfields) + '\tnum_reads')

# get index of taxid, qseqid:
taxidindex = myfields.index('staxids')
qseqidindex = myfields.index('qseqid')

# list of files
myfiles = [args.input]
if args.input2: myfiles.append(args.input2)

for myfile in myfiles:
    with open(myfile, 'r') as f:
        for line in f:
            # get desired fields
            fields = [line.split('\t')[i].strip() for i in myindicies]
            taxid = fields[taxidindex]
            qseqid = fields[qseqidindex]
            # get read counts
            readcounts = '-'
            if qseqid in idx: readcounts = idx[qseqid]
            # bypass human taxids
            if taxid == humantaxid:
                pass
            elif taxid not in filterlist:
                print(args.sample + '\t' + '\t'.join(fields) + '\t' + readcounts)
