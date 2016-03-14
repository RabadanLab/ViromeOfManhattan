#!/usr/bin/env python

'''
    Pandora 
    ~~~~~~
    Identification and Discovery of Tumor Associated Microbes via RNAseq
'''

__author__ = 'Rabadan Lab'
__version__ = 'Revision: 1.0'
__date__ = 'Date: 11-2015'

import argparse
import sys
import re
import subprocess
import os
import time
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

    # path of this script
    mycwd = os.path.dirname(os.path.realpath(__file__))

    # implement subcommands (Think git: git add, git commit, git push, etc)
    subparsers = parser.add_subparsers(help='sub-command help')

    # create the parser for the 'scan' command
    parser_scan = subparsers.add_parser('scan', help='run the pathogen discovery pipeline')
    parser_scan.add_argument('-r1', '--mate1', help='first RNAseq mate')
    parser_scan.add_argument('-r2', '--mate2', help='second RNAseq mate')
    parser_scan.add_argument('-sr', '--refstar', help='STAR host reference')
    parser_scan.add_argument('-br', '--refbowtie', help='bowtie2 host reference')
    parser_scan.add_argument('-db', '--blastdb', help='blast (nt) database (contigs are the query set)')
    parser_scan.add_argument('-pdb', '--pblastdb', help='blast protein (nr) database (ORFs are the query set)')
    parser_scan.add_argument('-ct', '--contigthreshold', default='500', help='threshold on contig length for blast (default: 500)')
    parser_scan.add_argument('-ot', '--orfthreshold', default='200', help='threshold on ORF length for protein blast (default: 200)')
    parser_scan.add_argument('-ob', '--orfblast', action='store_true', help='blast the ORFs to protein (nr) database (default: off)')
    parser_scan.add_argument('-bl', '--blacklist', default=mycwd + '/resources/blacklist.txt', help='A text file containing a list of non-pathogen taxids to ignore')
    parser_scan.add_argument('-gz', '--gzip', action='store_true', help='input fastq files are gzipped (default: off)')
    parser_scan.add_argument('--noerror', action='store_true', help='do not check for errors (default: off)')
    parser_scan.add_argument('-s', '--steps', default='12345', help='steps to run. The steps are as follows: \
      step 1: host separation, \
      step 2: assembly, \
      step 3: blast contigs, \
      step 4: orf discovery, \
      step 5: reporting (default: 12345 - i.e, steps 1 through 5).')
    parser_scan.set_defaults(which='scan')

    # create the parser for the 'aggregate' command
    parser_agg = subparsers.add_parser('aggregate', help='create report aggregated over multiple sample runs')
    parser_agg.set_defaults(which='aggregate')

    # add common arguments
    for i in [parser_scan, parser_agg]: add_common_args(i)

    args = parser.parse_args()

    # add key-value pairs to the args dict
    # directory where this script resides             
    vars(args)['scripts'] = mycwd

    # print args
    print(args)
    print

    # error checking
    if '1' in args.steps and ( (not args.mate1) or (not args.mate2) ):
        print('[ERROR] Need --mate1 and --mate2 arguments for Step 1')
	sys.exit(1)
    if '2' in args.steps and ( (not args.refstar) or (not args.refbowtie) ):
        print('[ERROR] Need --refstar and --refbowtie arguments for Step 2')
	sys.exit(1)
    if '3' in args.steps and (not args.blastdb):
        print('[ERROR] Need --blastdb argument for Step 3')
	sys.exit(1)
    if '4' in args.steps and args.orfblast and (not args.pblastdb):
        print('[ERROR] Need --pblastdb argument for Step 4 if blasting ORFs')
	sys.exit(1)
    if '4' in args.steps and (not args.orfblast) and args.pblastdb:
        print('[WARNING] --pblastdb argument supplied but boolean --orfblast is off')

    return args

# -------------------------------------

def main():
    '''Run the appropriate sub-command in the Pandora suite'''

    # dict which maps each subcommand name to its corresponding function (reference)
    d = {
             'scan': scan_main,
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
    # the SGE string must look like this: 'Your job 8379811 ("test") has been submitted'
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
             '1': ('qsub -N hsep', '{}/scripts/host_separation.sh --scripts {} -1 {} -2 {} --refstar {} --refbowtie {} --gzip {} --noclean {}'.format(
                      args.scripts,
                      args.scripts,
                      args.mate1,
                      args.mate2,
                      args.refstar,
                      args.refbowtie,
                      int(args.gzip),
                      int(args.noclean) )
                  ),
             '2': ('qsub -N asm', '{}/scripts/assembly.sh --scripts {} --noclean {}'.format(
                      args.scripts,
                      args.scripts,
                      int(args.noclean) )
                  ),
             '3': ('qsub -N blst', '{}/scripts/blast_wrapper.sh --scripts {} --threshold {} --db {} --id {} --noclean {} --nosge {}'.format(
                      args.scripts,
                      args.scripts,
                      args.contigthreshold,
                      args.blastdb,
                      args.identifier,
                      int(args.noclean),
                      int(args.noSGE) )
                  ),
             '4': ('qsub -N orf', '{}/scripts/orf_discovery.sh --scripts {} --id {} --threshold {} --db {} --blast {} --noclean {}'.format(
                      args.scripts,
                      args.scripts,
                      args.identifier,
                      args.orfthreshold,
                      args.pblastdb,
                      int(args.orfblast),
                      int(args.noclean) )
                  ),
             '5': ('qsub -N rep', '{}/scripts/reporting.sh --scripts {} --id {} --blacklist {}'.format(
                      args.scripts,
                      args.scripts,
                      args.identifier,
                      args.blacklist )
                  )
    }

    # run steps
    for i in args.steps:
        jid = docmd(d[i], jid, args)
        print('Step ' + i + ', jid = ' + jid)

# -------------------------------------

def agg_main(args):
    '''Run aggregate function'''

    # dict which maps each step to 2-tuple, which contains the qsub part of the command,
    # and the shell part of the command
    #d = {
    #         '1': ('qsub -N agg', '{}/scripts/aggregate.sh {}'.format(args.scripts, int(args.noclean)))
    #}

    #jid = 0
    #jid = docmd(d['1'], jid, args)

    pass

# -------------------------------------

def check_error(args):
    '''Check for errors, check dependencies '''

    # check for required programs
    #helpers.check_dependencies(['samtools', 'bam', 'bowtie2', 'STAR', 'blastn', 'Trinity'])

    # check for existence of files, if supplied
    for i in [args.mate1, args.mate2, args.blacklist]:
	if i:
            helpers.check_path(i)

    # check if input files gzipped
    if args.mate1 and args.mate2:
        if args.gzip and not (args.mate1[-3:] == '.gz' and args.mate2[-3:] == '.gz'):
            print('[ERROR] For --gzip option, files must have .gz extension')
    	    sys.exit(1)
        elif (args.mate1[-3:] == '.gz' or args.mate2[-3:] == '.gz') and not args.gzip:
            print('[ERROR] Files have .gz extension: use --gzip option')
    	    sys.exit(1)

# -------------------------------------

if __name__ == '__main__':

    main()
