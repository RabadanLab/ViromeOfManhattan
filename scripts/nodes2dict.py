#!/usr/bin/env python

import sys, pickle

# This script associates each taxid with
# its parent taxid of a specific rank (e.g., phylum)
# It can then make a blacklist of uninteresting taxid id 
# (e.g., all taxid whose phylum == chordata)  

# arg1 is the NCBI nodes.dmp file
# download it: wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
nodesfile = sys.argv[1]

# taxid to (parent taxid, current rank)
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

# Now the goal is to find the level (e.g., phylum)
# of the set of all taxids - i.e., map the taxid to 
# the taxid of its parent node at a certain rank

# set the level
# level = 'order'
# level = 'class'
level = 'phylum'

# taxid to level
id2level = {}

print('map taxids to ' + level)

# make the dict
for i in alltaxid:
    # check if id in dict
    if i in id2parentrank:
	# if so, get the rank
        rank = id2parentrank[i][1]
	# new taxid
	j = i

        if verbose:
	    print('taxid: '),
	    print(i)
	    print('rank: '),
            print(id2parentrank[i])

        # if rank is not final level, move up the tree
        while rank != level and j != '1':
            # parent taxid
	    j = id2parentrank[j][0]
	    if verbose: print('parent taxid: ' + j)
	    # find rank of parent taxid
	    if j in id2parentrank:
                rank = id2parentrank[j][1]
		if verbose: print('parent rank: ' + id2parentrank[j])
	    else:
		if verbose: print(i + ' not found')
		break

        # if rank matches level, update dict
	if (rank == level) and (i not in id2level):
	    id2level[i] = j
	    if verbose:
	        print('id2level: '),
	        print(id2level)

# print('resultant dict')
# print(id2level)

print('pickle taxid to ' + level + ' dict')

# pickle
with open('taxid2' + level + '.pkl', 'wb') as handle:
    pickle.dump(id2level, handle)

# get list of taxids corresponding to Chordata (Taxonomy ID: 7711, Rank: phylum)
chordata = [x for x,y in id2level.iteritems() if y == '7711']

# print(chordata)

with open('chordatalist.pkl', 'wb') as handle:
    pickle.dump(chordata, handle)
