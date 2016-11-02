#!/usr/bin/env python

"""
    Helper functions for the assembly step
    ~~~~~~
"""

import sys
import numpy as np
import pandas as pd
from scipy.stats import entropy

# -------------------------------------

def norm_entropy(count_list):
    count_vector = np.array(count_list)
    prob_vector = count_vector / float(count_vector.sum())
    prob_uniform = np.array([1.0/len(prob_vector)] * len(prob_vector))
    H_norm = entropy(prob_vector) / entropy(prob_uniform)
    return H_norm 

def computedistrib(infile, outfile):
    """compute simple distribution of contigs"""

    # a list of contig lengths
    x = []

    with open(infile, 'r') as g:
        for line in g:
            if line[0] != '>': 
                x.append(len(line.rstrip()))

    from collections import Counter

    # count occurences of each length
    xcount = Counter(x)

    # tot length
    tot = sum(dict(xcount).values())

    # running sum
    runningsum = 0

    with open(outfile, 'w') as f:
        for i in sorted(dict(xcount), reverse=True):
            runningsum += xcount[i]
            f.write(str(i) + '\t')
            f.write(str(xcount[i]) + '\t')
            f.write(str(runningsum) + '/' + str(tot) + '\t')
            f.write(str(int(100*runningsum/tot)) + '%')
            f.write('\n')
        
# -------------------------------------

def formatpileup(infile, idxfile, outfile):
    """format the pileup file for computing entropy"""

    # id 2 length dict (output of samtools idxstats)
    idx = {}

    # load idx file
    with open(idxfile, 'r') as f:
        for line in f:
            # map id to length of contig
            idx[line.split()[0].strip()] = line.split()[1].strip()

    myid = ''    # contig id
    pos = ''    # position

    with open(outfile, 'w') as f:
        with open(infile, 'r') as g:
            for line in g:
                # get number reads 
                numrds = line.split()[3]

                # if beginning of a new contig (id != previous id)
                if line.split()[0] != myid:
                    # if change (and not first contig), check if previous contig was covered until the end
                    # if not covered, pad with zeros
                    if myid: 
                        if int(idx[myid]) > int(pos): 
                            for i in range(int(pos) + 1, int(idx[myid]) + 1): 
                                f.write(myid + '\t' + str(i) + '\t0\n')
 
                    # if contig starts at postion > 1, pad with zeros 
                    if int(line.split()[1]) > 1: 
                        for i in range(1, int(line.split()[1])): 
                            f.write(line.split()[0] + '\t' + str(i) + '\t0\n')

                    # set new id 
                    myid = line.split()[0]

                    # write current line
                    f.write(myid + '\t' + line.split()[1] + '\t' + numrds + '\n')

                # if discontinuity (position - previous position > 1), pad with zeros
                elif (int(line.split()[1]) - int(pos)) > 1:
                    for i in range(int(pos) + 1, int(line.split()[1])):
                        f.write(myid + '\t' + str(i) + '\t0\n')

                    f.write(myid + '\t' + line.split()[1] + '\t' + numrds + '\n')

                # otherwise, simply write line
                else:
                    f.write(myid + '\t' + line.split()[1] + '\t' + numrds + '\n')

                # get position (this will become previous position for next iteration)
                pos = line.split()[1]

        # check if last contig covered until the end
        if int(idx[myid]) > int(pos):
            for i in range(int(pos) + 1, int(idx[myid]) + 1):
                f.write(myid + '\t' + str(i) + '\t0\n')

# -------------------------------------

if __name__ == "__main__":

    # to execute as a stand-alone script, give the name of the function 
    # as the first arg, followed by the args to the function
    if sys.argv[1] in globals():
        try:
            globals()[sys.argv[1]](*sys.argv[2:])
	except:
            print('Error! Are you sure you\'re using the correct arguments?')
    else:
        print('Function not found')
