#!/usr/bin/env python

import argparse
import sys
import os

# This script generates the report

# -------------------------------------

def get_arg():
    """Get Arguments
    :rtype: object
    """

    # parse arguments
    prog_description = 'Make tsv report'
    parser = argparse.ArgumentParser(description=prog_description)
    parser.add_argument('-i', '--input', default='blast/top.concat.txt', help='the contig blast file input')
    parser.add_argument('-i2', '--input2', default='discovery/blast/top.concat.txt', help='the ORF blast file input')
    parser.add_argument('-o', '--outputdir', default='report', help='the output directory')
    parser.add_argument('-l', '--logsdir', help='the logs directory')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('--header', default='blast/header', help='the blast header file')
    parser.add_argument('--id2reads', default='assembly/reads2contigs.stats.txt', help='the output file of samtools idxstats mapping ids to #reads')
    parser.add_argument('--blacklist', help='the file of blacklist taxids')
    parser.add_argument('--noclean', type=int, default=0, help='do not delete temporary intermediate files (default: off)')
    parser.add_argument('--verbose', type=int, default=0, help='verbose mode: echo commands, etc (default: off)')
    parser.add_argument('--id', required=True, help='sample identifier')
    args = parser.parse_args()

    # add key-value pairs to the args dict
    vars(args)['step'] = 'reporting'

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # error checking: exit if previous step produced zero output
    for i in [args.input, args.header, args.id2reads]:
        hp.check_file_exists_and_nonzero(i, step=args.step)

    return args

# -------------------------------------

def makerep(args):
    """Make report"""

    # The goal of this function is to filter blast results based on taxid
    # Don't want human sequences or anything in the taxid blacklist,
    # which contains non-pathogens

    hp.echostep(args.step)

    # print args
    print(args)
    print

    # mkdir -p
    hp.mkdirp(args.outputdir)

    # defaults
    # taxid corresponding to homo sapiens
    humantaxid = '9606'
    # desired header
    desiredfields = ['qseqid', 'sseqid', 'qlen','saccver','staxids','evalue', 'bitscore', 'stitle']

    # taxid blacklist
    filterlist = []
    # header
    header = []
    # id 2 reads dict (output of samtools idxstats)
    idx = {}

    # load blacklist if supplied
    if args.blacklist:
        with open(args.blacklist, 'r') as f:
            # final newline causes empty list elt
            filterlist = f.read().split('\n')[:-1]

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

    # get index of taxid, qseqid:
    taxidindex = myfields.index('staxids')
    qseqidindex = myfields.index('qseqid')

    # list of files
    myfiles = [args.input]
    if hp.check_path_bool(args.input2):
        myfiles.append(args.input2)

    with open(args.outputdir + '/blast.topfilter.unsort.txt', 'w') as f:
        # print header:
        f.write('sampleid\t' + '\t'.join(myfields) + '\tnum_reads\n')
        # loop through inputs
        for myfile in myfiles:
            with open(myfile, 'r') as g:
                for line in g:
                    # don't want lines with predicted (as opposed to real) genes
                    if not 'PREDICTED' in line:
                        # get desired fields
                        fields = [line.split('\t')[i].strip() for i in myindicies]
                        taxid = fields[taxidindex]
                        qseqid = fields[qseqidindex]
                        # get read counts
                        readcounts = '-'
                        if qseqid in idx:
                            readcounts = idx[qseqid]
                        # bypass human taxids
                        if taxid == humantaxid:
                            pass
                        elif taxid not in filterlist:
                            f.write(args.id + '\t' + '\t'.join(fields) + '\t' + readcounts + '\n')

    # sort by staxids then qlen (with bash)
    # careful: you're including the header in the file (make sure it's sorted properly)
    cmd = 'sort -k5,5n -k6,6nr {} > {}'.format(
              args.outputdir + '/blast.topfilter.unsort.txt',
              args.outputdir + '/blast.topfilter.txt',
    )
    hp.run_cmd(cmd, args.verbose, 0)

    os.remove(args.outputdir + '/blast.topfilter.unsort.txt')

    hp.echostep(args.step, start=0)

# -------------------------------------

def main():
    """Main function"""

    # get arguments
    args = get_arg()
    # make report
    makerep(args)

# -------------------------------------

if __name__ == '__main__':

    main()
