#!/usr/bin/env python

import sys, pickle, os

# The goal of this script is to filter blast results based on taxid
# Don't want human sequences or anything in the taxid blacklist,
# which contains non-pathogens

mypkl = sys.argv[1]		# arg 1 is the pickle
blastheader = sys.argv[2]	# arg 2 is the blast header
blastout = sys.argv[3]		# arg 3 is the blast file

# taxid blacklist
filterlist = []

# header
header = []
# desired header
desiredfields = ['qseqid', 'sseqid', 'qlen','saccver','staxids','evalue', 'bitscore', 'stitle']

# load pickle
with open(mypkl, 'rb') as handle:
    # taxid to (parent taxid, rank)
    filterlist = pickle.load(handle)

# load header
with open(blastheader, 'rb') as f:
    # taxid to (parent taxid, rank)
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
