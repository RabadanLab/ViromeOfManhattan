#!/usr/bin/env python
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=2G,time=2::

import argparse
import sys

# This script finds ORFs in the contigs which didn't blast

# -------------------------------------

def get_arg():
    """Get Arguments
    :rtype: object
    """
    # parse arguments

    prog_description = 'Find ORFs'
    parser = argparse.ArgumentParser(description=prog_description)
    parser.add_argument('-i', '--input', default='blast/no_blastn.fa', help='the input fasta')
    parser.add_argument('-o', '--outputdir', default='discovery', help='the output directory')
    parser.add_argument('-l', '--logsdir', help='the logs directory')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('--noclean', type=int, default=0, help='do not delete temporary intermediate files (default: off)')
    parser.add_argument('--verbose', type=int, default=0, help='verbose mode: echo commands, etc (default: off)')
    parser.add_argument('--blast', type=int, default=0, help='blast ORFs (default: off)')
    parser.add_argument('--threshold', type=int, default=0, help='the ORF length threshold')
    parser.add_argument('--db', help='the database prefix')
    parser.add_argument('--nosge', type=int, default=0, help='no SGE bool')
    parser.add_argument('--id', help='id')
    args = parser.parse_args()

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

def discovery(args):
    """ORF Discovery"""

    print('------------------------------------------------------------------')
    print('DISCOVERY START')

    # mkdir -p 
    hp.mkdirp(args.outputdir)

    # find ORFs
    hp.getorf(args.input, args.outputdir + '/orf.fa', args.threshold)

    # check if output exists
    hp.check_file_exists_and_nonzero(args.outputdir + '/orf.fa')

    # if blast discovered ORFs
    if args.blast:

        # make directory
        hp.mkdirp(args.outputdir + '/blast')

        # define command: blastp to nr, if blast flag
        cmd = '{}/scripts/blast_wrapper.sh --scripts {} --outputdir {} -i {} --logsdir {} --whichblast {} --threshold {} --db {} --id {} --noclean {}'.format(
                  args.scripts,
                  args.scripts,
                  args.outputdir  + '/blast',
                  args.outputdir + '/orf.fa',
                  'logs_blast2',
                  'blastp',
                  100,
                  args.db,
                  args.id,
                  args.noclean,
        )
        hp.run_cmd(cmd, args.verbose, 0)

    print('DISCOVERY END')
    print('------------------------------------------------------------------')

# -------------------------------------

def main():
    """Main function"""

    # get arguments
    args = get_arg()
    # ORF discovery
    discovery(args)

# -------------------------------------

if __name__ == '__main__':

    main()
