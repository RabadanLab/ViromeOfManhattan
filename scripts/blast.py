#!/usr/bin/env python
#$ -V
#$ -cwd
#$ -l mem=4G,time=4::

# A simple wrapper for blast array job 

import argparse
import sys
import os

# -------------------------------------

def get_arg():
    """Get Arguments
    :rtype: object
    """
    # parse arguments

    prog_description = 'A wrapper for blast'
    parser = argparse.ArgumentParser(description=prog_description)
    # parser.add_argument('-i', '--input', help='the input fasta')
    parser.add_argument('-o', '--outputdir', default='blast', help='the output directory')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('--sgeid', default=os.environ['SGE_TASK_ID'], help='SGE_TASK_ID (set by hand only if qsub is turned off)')
    parser.add_argument('--whichblast', default='blastn', choices=['blastn', 'blastp'], help='which blast to use (blastn, blastp)')
    parser.add_argument('--db', help='the database prefix')
    parser.add_argument('--fmt', help='blast format string')
    parser.add_argument('--noclean', type=int, default=0, help='do not delete temporary intermediate files (default: off)')
    parser.add_argument('--verbose', type=int, default=0, help='verbose mode: echo commands, etc (default: off)')
    args = parser.parse_args()

    # add key-value pairs to the args dict
    vars(args)['input'] = args.outputdir + '/blast_' + args.sgeid + '.fasta'
    vars(args)['step'] = 'blast'

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # error checking: exit if previous step produced zero output
    for i in [args.input]:
        hp.check_file_exists_and_nonzero(i, step=args.step)

    return args

# -------------------------------------

def blastnp(args):
    """Blastn or blastp"""

    hp.echostep(args.step)

    # print(args)
    # print

    flag = ''

    if args.whichblast == 'blastp':
        flag = '-task blastp-fast'

    # do blastn or blastp
    cmd = '{} -outfmt "6 {}" -query {} -db {} {} > {}/blast_{}.result'.format(
        args.whichblast,
        args.fmt,
        args.input,
        args.db,
        flag,
        args.outputdir,
        args.sgeid
    )
    hp.run_cmd(cmd, args.verbose, 0)

    hp.echostep(args.step, start=0)

# -------------------------------------

def main():
    """Main function"""

    # get arguments
    args = get_arg()
    # blast
    blastnp(args)

# -------------------------------------

if __name__ == '__main__':

    main()
