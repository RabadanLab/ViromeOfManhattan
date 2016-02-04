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

def add_common_args(sub):
    '''A function to add common args to subparers, modified from StackOverflow'''

    # http://stackoverflow.com/questions/33463052/how-do-i-specify-two-required-arguments-including-a-subcommand-using-argparse

    # add common args
    sub.add_argument('-id', '--identifier', required=True, help='sample ID (5 chars or less)')
    sub.add_argument("--noclean", action="store_true", help="do not delete temporary intermediate files (default: off)")
    sub.add_argument('--verbose', action='store_true', help='verbose mode: echo commands, etc (default: off)')
    sub.add_argument("--noSGE", action="store_true", help="do not qsub jobs with the Oracle Grid Engine (default: off)")

    return sub

# -------------------------------------

def get_arg():
    '''Get Arguments'''

    prog_description = 'microbial detection from paired-end RNAseq'
    parser = argparse.ArgumentParser(description=prog_description)

    # implement subcommands (Think git: git add, git commit, git push, etc)
    subparsers = parser.add_subparsers(help='sub-command help')

    # create the parser for the 'scan' command
    parser_scan = subparsers.add_parser('scan', help='run the pathogen discovery pipeline')
    parser_scan.add_argument('-r1', '--mate1', required=True, help='first RNAseq mate')
    parser_scan.add_argument('-r2', '--mate2', required=True, help='second RNAseq mate')
    parser_scan.add_argument('-sr', '--refstar', required=True, help='STAR host reference')
    parser_scan.add_argument('-br', '--refbowtie', required=True, help='bowtie2 host reference')
    parser_scan.add_argument('-db', '--blastdb', required=True, help='blast (nt) database')
    parser_scan.add_argument('-ct', '--contigthreshold', default='500', help='threshold on contig length for blast (default: 500)')
    parser_scan.add_argument('-bl', '--blacklist', help='A text file containing a list of non-pathogen taxids to ignore')
    parser.add_argument("--noerror", action="store_true", help="do not check for errors (default: off)")
    # parser_scan.add_argument("--remap", action="store_true", help="create fasta file of pathogen sequences and map reads back onto this reference (default: off)")
    parser_scan.add_argument('-s', '--steps', default='12345', help='steps to run. The steps are as follows: \
        step 1: host separation, step 2: assembly, step 3: blast contigs, step 4: orf discovery, step 5: reporting (default: 12345 - i.e, steps 1 through 5).')
    parser_scan.set_defaults(which='scan')

    # create the parser for the 'remap' command
    parser_remap = subparsers.add_parser('remap', help='map reads back onto contigs')
    parser_remap.set_defaults(which='remap')

    # create the parser for the 'aggregate' command
    parser_agg = subparsers.add_parser('aggregate', help='create report aggregated over multiple sample runs')
    parser_agg.set_defaults(which='aggregate')

    # add common arguments
    for i in [parser_scan, parser_remap, parser_agg]: add_common_args(i)

    args = parser.parse_args()

    # add key-value pairs to the args dict
    # directory where this script resides             
    vars(args)['scripts'] = os.path.dirname(os.path.realpath(__file__))

    # print args
    print(args)
    print

    return args

# -------------------------------------

def main():
    '''Run the appropriate sub-command in the Pandora suite'''

    # dict which maps each subcommand name to its corresponding function (reference)
    d = {
             'scan': scan_main,
             'remap': remap_main,
             'aggregate': agg_main
    }

    # get arguments
    args = get_arg()

    # invoke subcommand
    # print('--> subcommand: ' + args.which)
    d[args.which](args)

# -------------------------------------

def getjid(x):
    '''Parse out and return SGE job id from string'''
    # string looks like this: 'Your job 8379811 ("test") has been submitted'
    return x.split('Your job ')[1].split()[0]

# -------------------------------------

def docmd(mytuple, jid, args):
    '''Run a command on the shell or with SGE qsub'''

    # mytuple - a tuple containing (qsub part, shell part) for a given command
    # jid - job id
    # args - args dict

    # if run in the shell without qsub
    if args.noSGE:
        cmd = mytuple[1]
        # if verbose, print command
        if args.verbose: print(cmd)
        subprocess.check_output(cmd, shell=True)
	return '0'
    # if run command with SGE qsub
    else:
        # define qsub part of command
        cmd = mytuple[0] + '_' + args.identifier + ' '
        # if not the first command, hold on previous job id
        if jid > 0: cmd += '-hold_jid ' + jid + ' '
        # define shell (non-qsub) part of command
        cmd += mytuple[1]
        # if verbose, print command
        if args.verbose: print(cmd)
        # run command, get job id
	return getjid(subprocess.check_output(cmd, shell=True))

# -------------------------------------

def scan_main(args):
    '''Run pathogen discovery steps'''

    # check for errors
    if not args.noerror: check_error(args)

    # start with job id set to zero
    jid = 0

    # dict which maps each step to 2-tuple, which contains the qsub part of the command,
    # and the shell part of the command
    d = {
             '1': ('qsub -N hsep', '{}/scripts/host_separation.sh {} {} {} {} {}'.format(args.scripts, args.mate1, args.mate2, args.refstar, args.refbowtie, int(args.noclean))),
             '2': ('qsub -N asm', '{}/scripts/assembly.sh {} {}'.format(args.scripts, int(args.noclean), args.scripts)),
             '3': ('qsub -N blst', '{}/scripts/blast_contigs.sh {} {} {} {} {} {}'.format(args.scripts, args.contigthreshold, args.blastdb, args.identifier, args.scripts, int(args.noclean), int(args.noSGE))),
             '4': ('qsub -N orf', '{}/scripts/orf_discovery.sh'.format(args.scripts)),
             '5': ('qsub -N rep', '{}/scripts/reporting.sh {} {} {}'.format(args.scripts, args.scripts, args.identifier, args.blacklist))
    }

    # run steps
    for i in args.steps:
        jid = docmd(d[i], jid, args)
        print('Step ' + i + ', jid = ' + jid)

# -------------------------------------

def remap_main(args):
    '''Run remap function'''

    # dict which maps each step to 2-tuple, which contains the qsub part of the command,
    # and the shell part of the command
    d = {
             '1': ('qsub -N rmap', '{}/scripts/remap.sh {}'.format(args.scripts, int(args.noclean)))
    }

    jid = 0
    jid = docmd(d['1'], jid, args)

# -------------------------------------

def agg_main(args):
    '''Run aggregate function'''

    pass

# -------------------------------------

def check_error(args):
    '''Check for errors, check dependencies '''

    # check for required programs
    #helpers.check_dependencies(['samtools', 'bam', 'bowtie2', 'STAR', 'blastn', 'Trinity', 'prodigal'])

    # check for existence of files, if supplied
    for i in [args.mate1, args.mate2, args.blacklist]:
	if i:
            helpers.check_path(i)

# -------------------------------------

if __name__ == '__main__':

    main()
