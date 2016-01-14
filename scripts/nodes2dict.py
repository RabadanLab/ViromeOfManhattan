#!/usr/bin/env python

import sys, pickle

# This script associates each taxid with
# its parent taxid of a specific rank (e.g., phylum)
# It can then make a blacklist of uninteresting taxid id 
# (e.g., all taxid whose phylum == chordata or who
# belong to the 'other sequences' lineage)

# arg1 is the NCBI nodes.dmp file
# download it: wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
nodesfile = sys.argv[1]

# taxid to (parent taxid, rank)
id2parentrank = {}

# verbose
verbose = 0

# list of all taxids
alltaxid = []

print('map taxids to parent taxids')

# open the nodes.dmp file
with open(nodesfile, 'r') as f:
    for line in f:
	# eliminate tabs and split on |
	l = line.replace('\t', '').split('|')
        # add taxid to master list
        alltaxid.append(l[0])
	# map taxid to (parent taxid, rank)
	id2parentrank[l[0]] = (l[1], l[2])

print('pickle taxid to parent dict')

# pickle dict for possible later use
with open('taxid2parent.pkl', 'wb') as handle:
    pickle.dump(id2parentrank, handle)

print('generate lists')

# list of taxids in the lineage of other sequences (artificial and synthetic)
# http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=28384
otherseq = []

# list of taxids in the lineage of chordata (Taxonomy ID: 7711, Rank: phylum)
# http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=7711
chordata = []

# make the dict
for i in alltaxid:
    # check if id in dict
    if i in id2parentrank:
	# copy current taxid
	j = i

        if verbose:
	    print('taxid: '),
	    print(i)
	    print('rank: '),
            print(id2parentrank[i])

        # if j not 1, move up the tree
        while j != '1':
            # parent taxid
	    j = id2parentrank[j][0]

	    if verbose: print('parent taxid: ' + j)

            # if member of "other sequences" (Taxonomy ID: 28384), add to list and break
            if j == '28384':
                otherseq.append(i)
		break
	    # if member of chordata
            elif j == '7711':
                chordata.append(i)
		break
	    # if not in dict, break to avoid infinite loop
	    elif not (j in id2parentrank):
		if verbose: print(i + ' not found')
		break

print('pickle')

# chordata and 'other sequences' are not pathogens, so they get blacklisted
with open('blacklist.pkl', 'wb') as handle:
    pickle.dump(chordata + otherseq, handle)
