#!/usr/bin/env python

import argparse
import sys

# -------------------------------------

def get_arg():
    """
    Get Arguments
    :rtype: object
    """
    # parse arguments

    prog_description = 'preprocess'
    parser = argparse.ArgumentParser(description=prog_description)
    parser.add_argument('--taxid2names', help='location of names.dmp file mapping taxid to names')
    parser.add_argument('--taxid2nodes', help='location of nodes.dmp file')
    parser.add_argument('--samples', help='path of the file containing the samples names (one sample per line)')
    parser.add_argument('--batchdir', help='path of the directory containing the output of multiple Pandoras runs')    
    parser.add_argument('--suffixreport', help='suffix string of the report (default: /report_ifilter/report.contig.txt)')
    parser.add_argument('--suffixstats', help='suffix string of the stats report (default: /report_ifilter/report.taxon.txt)')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('-o', '--outputdir', default='aggregate_preprocess', help='the output directory')
    args = parser.parse_args()

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # add key-value pairs to the args dict
    vars(args)['step'] = 'preprocess'
    vars(args)['olog'] = args.outputdir + '/../' + 'log.out'
    vars(args)['elog'] = args.outputdir + '/../' + 'log.err'

    return args

# -------------------------------------

def get_taxid_stuff(args):
    """
    """

    # taxid to (parent taxid, rank)
    id2parentrank = {}
    # taxid to name 
    id2name = {}

    # don't want to loop more than this many iterations through the taxonomic heirarchy
    # in case taxid not found 
    maxiterations = 50 

    # open the nodes.dmp file
    with open(args.taxid2nodes, 'r') as f:
        for line in f:
            # eliminate tabs and split on |
            l = line.replace('\t', '').split('|')
            # map taxid to (parent taxid, rank)
            id2parentrank[l[0]] = (l[1], l[2])

    # open the names.dmp file
    with open(args.taxid2names, 'r') as f:
        for line in f:
            l = line.replace('\t', '').split('|')
            if (l[3] == 'scientific name'):
                id2name[l[0]] = l[1]

    return (id2parentrank, id2name)

# -------------------------------------
  
def taxonomy(taxid, id2parentrank, id2name):

    """
    Function which returns the taxonomic tree traversal as a dictionary for a single taxid
    
    INPUT: taxid string
    OUTPUT: dictionary of string-string pairs with taxonomic levels and names
    """

    # don't want to loop more than this many iterations through the taxonomic heirarchy
    # in case taxid not found
    maxiterations = 50 

    tax_tree = {}

    myiterations = 0
    current_taxid = taxid.strip()
    while (current_taxid != '1' and myiterations < maxiterations):
        
        tax_tree[id2parentrank.get(current_taxid, ['', 'not found'])[1]] = [id2name.get(current_taxid, 'not found'), current_taxid]
        
        # get parent taxid and iterate up the tree of life
        current_taxid = id2parentrank.get(current_taxid, ['1'])[0] 
        
        myiterations += 1

    return tax_tree

# -------------------------------------

# -------------------------------------

# -------------------------------------

# -------------------------------------


def main():
    """Main function"""

    # get arguments
    args = get_arg()

    hp.echostep(args.step)

    # print args
    print(args)
    print

    # mkdir -p
    hp.mkdirp(args.outputdir)

    (id2parentrank, id2name) = get_taxid_stuff(args)

    # example
    print(taxonomy("24", id2parentrank, id2name))

    # end of step
    hp.echostep(args.step, start=0)

# -------------------------------------

if __name__ == '__main__':

    main()
