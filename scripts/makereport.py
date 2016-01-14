#!/usr/bin/env python

import sys, pickle, os

# The goal of this script is to filter blast results based on taxid
# Don't want human sequences or anything in the taxid blacklist,
# which contains non-pathogens

# arg 1 is the pickle
mypkl = sys.argv[1]
# arg 2 is the blast file
blastout = sys.argv[2]

try:
    os.mkdir('report')
except:
    pass

# taxid blacklist
filterlist = []

# load pickle
with open(mypkl, 'rb') as handle:
    # taxid to (parent taxid, rank)
    filterlist = pickle.load(handle)

# with open('report/blast.topfilter.txt', 'w') as g:
with open(blastout, 'r') as f:
    for line in f:
        taxid = line.split()[3]
	# bypass human taxids
        if taxid == '9606':
            pass
        elif taxid not in filterlist:
            # g.write(line)
            print(line),
