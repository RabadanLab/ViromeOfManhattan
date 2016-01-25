#!/usr/bin/env python

'''
    Pandora 
    ~~~~~~
    Identification and Discovery of Tumor Associated Microbes via RNAseq
'''

__author__ = 'Rabadan Lab'
__version__ = 'Revision: 1.0'
__date__ = 'Date: 11-2015'

import argparse, sys, re, subprocess, os, time
from helpers import helpers
from distutils import spawn

# -------------------------------------

def get_arg():
    '''Get Arguments'''

    prog_description = 'microbial detection from paired-end RNAseq'
    parser = argparse.ArgumentParser(description=prog_description)

    parser.add_argument('-id', '--identifier', required=True, help='5 chars or less sample ID')
    parser.add_argument('-r1', '--mate1', required=True, help='first RNAseq mate')
    parser.add_argument('-r2', '--mate2', required=True, help='second RNAseq mate')
    parser.add_argument('-ct', '--contigthreshold', required=True, help='threshold on contig length for blast')
    parser.add_argument('-sr', '--refstar', required=True, help='STAR host reference')
    parser.add_argument('-br', '--refbowtie', required=True, help='bowtie2 host reference')
    parser.add_argument('-db', '--db', required=True, help='blast (nt) database')
    parser.add_argument('-bl', '--blacklist', help='A text file containing a list of non-pathogen taxids to ignore')
    parser.add_argument('-s', '--steps', default='12345', help='steps to run (default: 12345 - i.e, steps 1 through 5')
    parser.add_argument("--remap", action="store_true", help="create fasta file of pathogen sequences and map reads back onto this reference (default: off)")
    # parser.add_argument("--noerror", action="store_true", help="do not check for errors (default: off)")
    parser.add_argument("--noclean", action="store_true", help="do not delete temporary intermediate files (default: off)")
    parser.add_argument('--verbose', action='store_true', help='verbose mode: echo commands, etc (default: off)')

    args = parser.parse_args()

    # add key-value pairs to the args dict
    # directory where this script resides             
    vars(args)['scripts'] = os.path.dirname(os.path.realpath(__file__))

    # print args
    print(args)
    print

    return args

# -------------------------------------

def getjid(x):
    '''Parse out and return SGE job id from string'''
    # string looks like this: 'Your job 8379811 ("test") has been submitted'
    return x.split('Your job ')[1].split()[0]

# -------------------------------------

def main():
    '''Run all steps'''

    # get arguments dict
    args = get_arg()
    
    # check for errors
    #if not args.noerror: check_error(args)

    # start with job id set to zero
    jid = 0

    # dict which maps each step to 2-tuple, which contains the qsub part of the command,
    # and the shell part of the command
    d = {
             '1': ('qsub -N hsep', '{}/scripts/host_separation.sh {} {} {} {} {}'.format(args.scripts, args.mate1, args.mate2, args.refstar, args.refbowtie, int(args.noclean))),
             '2': ('qsub -N asm', '{}/scripts/assembly.sh {} {}'.format(args.scripts, int(args.noclean), args.scripts)),
             '3': ('qsub -N blst', '{}/scripts/blast_contigs.sh {} {} {} {} {}'.format(args.scripts, args.contigthreshold, args.db, args.identifier, args.scripts, int(args.noclean))),
             '4': ('qsub -N orf', '{}/scripts/orf_discovery.sh'.format(args.scripts)),
             '5': ('qsub -N rep', '{}/scripts/reporting.sh {} {} {} {} {}'.format(args.scripts, args.scripts, int(args.remap), args.db, int(args.noclean), args.blacklist))
    }

    # run steps
    for i in sorted(d):
        # if step requested
        if i in args.steps:
            # define qsub part of command
            cmd = d[i][0] + '_' + args.identifier + ' '
            # if not the first command, hold on previous job id
            if jid > 0: cmd += '-hold_jid ' + jid + ' '
            # define shell (non-qsub) part of command
            cmd += d[i][1]
            # if verbose, print command
            if args.verbose: print(cmd)
            # run command, get job id
            jid = getjid(subprocess.check_output(cmd, shell=True))
            print('Step ' + i + ', jid = ' + jid)

# -------------------------------------

def check_error(args):
    '''Check for errors, check dependencies '''

    # check for required programs
    helpers.check_dependencies(['samtools', 'bam', 'bowtie2', 'STAR', 'blastn', 'Trinity', 'prodigal'])

# -------------------------------------

if __name__ == '__main__':

    main()
