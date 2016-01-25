#!/usr/bin/env python

import sys, os

# The goal of this script is to filter blast results based on taxid
# Don't want human sequences or anything in the taxid blacklist,
# which contains non-pathogens

blastheader = sys.argv[1]	# the blast header
blastout = sys.argv[2]		# the blast file
# myblacklist = sys.argv[3]	# the file of blacklist taxids

# taxid blacklist
filterlist = []

# header
header = []
# desired header
desiredfields = ['qseqid', 'sseqid', 'qlen','saccver','staxids','evalue', 'bitscore', 'stitle']

# load blacklist if supplied (hacky - fix later)
if len(sys.argv) > 3 and sys.argv[3] != 'None':
    with open(sys.argv[3], 'r') as f:
        filterlist = f.read().split('\n')[:-1]	# final newline causes empty list elt

# load header
with open(blastheader, 'r') as f:
    header = f.read().split()

# get indicies of desired fields
myindicies = [(j,k) for j,k in enumerate(header) if k in desiredfields]
# looks something like:
# [(0, 'qseqid'), (1, 'sseqid'), (2, 'saccver'), (3, 'staxids'), (12, 'qlen'), (20, 'evalue'), (21, 'bitscore'), (22, 'stitle')]

# print header:
print('\t'.join([i[1] for i in myindicies]))

# get index of taxid:
taxidindex = [i[1] for i in myindicies].index('staxids')

with open(blastout, 'r') as f:
    for line in f:
        # get desired fields
        myfields = [line.split('\t')[i] for i in [j[0] for j in myindicies]] 
        taxid = myfields[taxidindex]
	# bypass human taxids
        if taxid == '9606':
            pass
        elif taxid not in filterlist:
            print('\t'.join(myfields)),
