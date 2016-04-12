#!/usr/bin/env python
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=4G,time=8::

import argparse
import sys

# This script blasts the entries of a fasta file,
# provided they're over the length threshold

# -------------------------------------

def get_arg():
    """Get Arguments
    :rtype: object
    """
    # parse arguments

    prog_description = 'Blast in parallel'
    parser = argparse.ArgumentParser(description=prog_description)
    parser.add_argument('-i', '--input', default='assembly/contigs_trinity.fasta', help='the input fasta')
    parser.add_argument('-o', '--outputdir', default='blast', help='the output directory')
    parser.add_argument('-l', '--logsdir', default='logs_blast', help='the logs directory')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('--noclean', help='do not delete temporary intermediate files (default: off)')
    parser.add_argument('--verbose', action='store_true', help='verbose mode: echo commands, etc (default: off)')
    parser.add_argument('--threshold', default=0, help='the length threshold')
    parser.add_argument('--db', help='the database prefix')
    parser.add_argument('--whichblast', default='blastn', choices=['blastn', 'blastp'], help='which blast to use (blastn, blastp)')
    parser.add_argument('--nosge', help='no SGE bool')
    parser.add_argument('--id', help='id')
    args = parser.parse_args()

    # blast format string
    fmt="qseqid sseqid saccver staxids pident nident length mismatch gapopen gaps qstart qend qlen qframe qcovs sstart send slen sframe sstrand evalue bitscore stitle"

    # add key-value pairs to the args dict
    vars(args)['fmt'] = fmt

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # print args
    print(args)
    print

    # error checking: exit if previous step produced zero output
    for i in [args.input]:
        hp.check_file_exists_and_nonzero(i)

    return args

# -------------------------------------

def blast(args):
    """Do blast in parallel"""

    # mkdir -p 
    hp.mkdirp(args.outputdir)
    hp.mkdirp(args.logsdir)

    # split fasta file on contigs above threshold length (and return count)
    filecount = hp.fastasplit(args.input, args.outputdir + '/blast', args.threshold)

    print('------------------------------------------------------------------')
    print('BLAST START')
    print('BLAST END')
    print('------------------------------------------------------------------')

# -------------------------------------

def main():
    """Main function"""

    # get arguments
    args = get_arg()
    # blast
    blast(args)

# -------------------------------------

if __name__ == '__main__':

    main()
