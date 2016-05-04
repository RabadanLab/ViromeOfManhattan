#!/usr/bin/env python

"""
    Helper functions for the assembly step
    ~~~~~~
"""

# -------------------------------------

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

    pass
