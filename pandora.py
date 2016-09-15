#!/usr/bin/env python

"""
    Pandora 
    ~~~~~~
    Identification and Discovery of Tumor Associated Microbes via RNAseq
"""

__author__ = 'Rabadan Lab'
__version__ = 'Revision: 1.0'
__date__ = 'Date: 11-2015'

import argparse
import sys
import subprocess
import os
import ConfigParser
from helpers import helpers as hp

def add_common_args(sub):
    """A function to add common args to subparers, modified from StackOverflow"""

    # http://stackoverflow.com/questions/33463052/how-do-i-specify-two-required-arguments-including-a-subcommand-using-argparse

    # add common args
    sub.add_argument('-id', '--identifier', required=True, help='sample ID (5 chars or less)')
    sub.add_argument('-c', '--config', help='config file')
    sub.add_argument('--noclean', action='store_true', help='do not delete temporary intermediate files (default: off)')
    sub.add_argument('--verbose', action='store_true', help='verbose mode: echo commands, etc (default: off)')
    sub.add_argument('--noSGE', action='store_true', help='do not qsub jobs with the Oracle Grid Engine (default: off)')

    return sub

# -------------------------------------

def get_arg():
    """Get Arguments"""

    prog_description = 'microbial detection from paired-end RNAseq'
    parser = argparse.ArgumentParser(description=prog_description)

    # path of this script
    mycwd = os.path.dirname(os.path.realpath(__file__))

    # implement subcommands (Think git: git add, git commit, git push, etc)
    subparsers = parser.add_subparsers(help='sub-command help')

    # create the parser for the 'scan' command
    parser_scan = subparsers.add_parser('scan', help='run the pathogen discovery pipeline')
    parser_scan.add_argument('-r1', '--mate1', default=None, help='RNA-seq mate 1 fastq input')
    parser_scan.add_argument('-r2', '--mate2', default=None, help='RNA-seq mate 2 fastq input')
    parser_scan.add_argument('--bam', default=None, help='bam file input (provide this as an alternative to fastq\'s')
    parser_scan.add_argument('-sr', '--refstar', help='STAR host reference')
    parser_scan.add_argument('-br', '--refbowtie', help='bowtie2 host reference')
    parser_scan.add_argument('-db', '--blastdb', help='blast (nt) database (contigs are the query set)')
    parser_scan.add_argument('-pdb', '--pblastdb', help='blast protein (nr) database (ORFs are the query set)')
    parser_scan.add_argument('-gtf', '--gtf', help='optional host gft for computing gene coverage after host separation')
    parser_scan.add_argument('--contigthreshold', default='500', help='threshold on contig length for blast (default: 500)')
    parser_scan.add_argument('--orfthreshold', default='200', help='threshold on ORF length for protein blast (default: 200)')
    parser_scan.add_argument('--orfblast', action='store_true', help='blast the ORFs to protein (nr) database (default: off)')
    parser_scan.add_argument('--blacklist', default=mycwd + '/resources/blacklist.txt', help='A text file containing a list of non-pathogen taxids to ignore')
    parser_scan.add_argument('--gzip', action='store_true', help='input fastq files are gzipped (default: off)')
    parser_scan.add_argument('--noerror', action='store_true', help='do not check for errors (default: off)')
    parser_scan.add_argument('--steps', default='12345', help='steps to run. The steps are as follows: \
      step 1: host separation, \
      step 2: assembly, \
      step 3: blast contigs, \
      step 4: orf discovery, \
      step 5: reporting (default: 12345 - i.e, steps 1 through 5).')
    parser_scan.add_argument('--trinitymem', default='50', help='max memory for Trinity in gigabytes (default: 50)')
    parser_scan.add_argument('--trinitycores', default='8', help='number of cores for Trinity (default: 8)')
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
    # qsub parameters
    vars(args)['qparam'] = {}

    # if configuration file specified and option not already supplied by flag, read option in from config file
    if args.config:
        Config = ConfigParser.ConfigParser()
        Config.read(args.config)
        for i in [(args.refstar, 'refstar', 'Step1'), 
                  (args.refbowtie, 'refbowtie', 'Step1'), 
                  (args.gtf, 'gtf', 'Step1'),
                  (args.blastdb, 'blastdb', 'Step3'), 
                  (args.pblastdb, 'pblastdb', 'Step4')]:
            if not i[0] and i[1] in hp.config_section_map(Config, i[2]):
                vars(args)[i[1]] = hp.config_section_map(Config, i[2])[i[1]]

        # get qsub params from config file, if specified
        for i in map(str, range(1,6)):
            if 'qparam' in hp.config_section_map(Config, 'Step' + i):
                vars(args)['qparam'][i] = hp.config_section_map(Config, 'Step' + i)['qparam']

    # print args
    print(args)
    print

    # error checking
    if '1' in args.steps and not ((args.mate1 and args.mate2) or args.bam):
        print('[ERROR] Need --mate1 and --mate2 arguments OR --bam argument for Step 1')
        sys.exit(1)
    if '1' in args.steps and ((not args.refstar) or (not args.refbowtie)):
        print('[ERROR] Need --refstar and --refbowtie arguments for Step 1')
        sys.exit(1)
    if '3' in args.steps and (not args.blastdb):
        print('[ERROR] Need --blastdb argument for Step 3')
        sys.exit(1)
    if '4' in args.steps and args.orfblast and (not args.pblastdb):
        print('[ERROR] Need --pblastdb argument for Step 4 if blasting ORFs')
        sys.exit(1)
    if '4' in args.steps and (not args.orfblast) and args.pblastdb:
        print('[WARNING] --pblastdb argument supplied but boolean --orfblast is off')

    # cast as various bools as ints
    for name in 'gzip verbose noclean noSGE orfblast'.split():
        setattr(args, name, int(getattr(args, name, 0)))

    # set bam file to abs path
    if args.bam: 
        setattr(args, 'bam', os.path.abspath(os.path.expanduser(args.bam)))

    return args

# -------------------------------------

def main():
    """Run the appropriate sub-command in the Pandora suite"""

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

def docmd(myqcmd, mycmd, jid, args):
    """Run a command on the shell or with SGE qsub"""

    # myqcmd - qsub portion (or prefix) of the command
    # mycmd - ordinary shell command
    # jid - job id
    # args - args dict

    # if run in the shell without qsub
    if args.noSGE:
        print(hp.run_long_cmd(mycmd, args.verbose, 'log.steps'))
        return '0'
    # if run command with SGE qsub
    else:
        # define qsub part of command
        # "sys.executable contains full path of the currently running Python interpreter"
        cmd = 'qsub -S ' + sys.executable + ' ' + myqcmd + ' '
        # if not the first command, hold on previous job id
        if jid != '0':
            cmd += '-hold_jid ' + jid + ' '
        # define shell (non-qsub) part of command
        cmd += mycmd
        # if verbose, print command
        if args.verbose:
            print(cmd)
        # run command, get job id
	return hp.getjid(subprocess.check_output(cmd, shell=True))

# -------------------------------------

def scan_main(args):
    """Run pathogen discovery steps"""

    # check for errors
    if not args.noerror:
        check_error(args)

    # dict which maps each step to the qsub part of the command
    q = {
             '1': '-N hsep_{args.identifier} -V -cwd -o log.out -e log.err -l mem=16G,time=12:: -pe smp 4 -R y'.format(args=args),
             '2': '-N asm_{args.identifier} -V -cwd -o log.out -e log.err -l mem=12G,time=12:: -pe smp 8 -R y'.format(args=args),
             '3': '-N blst_{args.identifier} -V -cwd -o log.out -e log.err -l mem=4G,time=8::'.format(args=args),
             '4': '-N orf_{args.identifier} -V -cwd -o log.out -e log.err -l mem=2G,time=2::'.format(args=args),
             '5': '-N rep_{args.identifier} -V -cwd -o log.out -e log.err -l mem=1G,time=1::'.format(args=args)
    }

    # dict which maps each step to the shell part of the command
    d = {
             '1': '{args.scripts}/scripts/host_separation.py --scripts {args.scripts} -1 {args.mate1} -2 {args.mate2} --bam {args.bam} --refstar {args.refstar} --refbowtie {args.refbowtie} --gzip {args.gzip} --verbose {args.verbose} --noclean {args.noclean} --gtf {args.gtf}'.format(args=args),

             '2': '{args.scripts}/scripts/assembly.py --scripts {args.scripts} --trinitymem {args.trinitymem} --trinitycores {args.trinitycores} --verbose {args.verbose} --noclean {args.noclean}'.format(args=args),

             '3': '{args.scripts}/scripts/blast_wrapper.py --scripts {args.scripts} --threshold {args.contigthreshold} --db {args.blastdb} --id {args.identifier} --verbose {args.verbose} --noclean {args.noclean} --nosge {args.noSGE}'.format(args=args),

             '4': '{args.scripts}/scripts/orf_discovery.py --scripts {args.scripts} --id {args.identifier} --threshold {args.orfthreshold} --db {args.pblastdb} --blast {args.orfblast} --verbose {args.verbose} --noclean {args.noclean}'.format(args=args),

             '5': '{args.scripts}/scripts/makereport.py --scripts {args.scripts} --id {args.identifier} --verbose {args.verbose} --blacklist {args.blacklist}'.format(args=args)
    }

    # start with job id set to zero string
    jid = '0'

    # run steps
    for i in args.steps:
        # if qsub params specified in config file
        if i in args.qparam:
            jid = docmd(args.qparam[i], d[i], jid, args)
        # if not, use default from dict q
        else:
            jid = docmd(q[i], d[i], jid, args)

        # only print step name if qsub-ing
        if not args.noSGE:
            print('Step ' + i + ', jid = ' + jid)

# -------------------------------------

def agg_main(args):
    """Run aggregate function"""

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
    """Check for errors, check dependencies"""

    # check for required programs
    #hp.check_dependencies(['samtools', 'bam', 'bowtie2', 'STAR', 'blastn', 'Trinity'])

    # check for existence of files, if supplied
    for i in [args.mate1, args.mate2, args.bam, args.blacklist]:
        if i:
            hp.check_path(i)

    # check if input files gzipped
    if args.mate1 and args.mate2:
        if args.gzip and not (args.mate1[-3:] == '.gz' and args.mate2[-3:] == '.gz'):
            print('[ERROR] For --gzip option, files must have .gz extension')
            sys.exit(1)
        elif (args.mate1[-3:] == '.gz' or args.mate2[-3:] == '.gz') and not args.gzip:
            print('[ERROR] Files have .gz extension: use --gzip option')
            sys.exit(1)

    # check if proper extention
    if args.bam:
        if not (args.bam[-4:] == '.bam'):
            print('[ERROR] For --bam option, files must have .bam extension')
            sys.exit(1)

# -------------------------------------

if __name__ == '__main__':

    main()
