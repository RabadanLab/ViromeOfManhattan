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

# -------------------------------------

def get_arg():
	'''Get Arguments'''

	prog_description = 'microbial detection from paired-end RNAseq'
	parser = argparse.ArgumentParser(description=prog_description)

	# directory where this script resides 			
	software = os.path.dirname(os.path.realpath(__file__))

	parser.add_argument('-id', '--identifier', required=True, help='5 chars or less sample ID')
	parser.add_argument('-r1', '--mate1', required=True, help='first RNAseq mate')
	parser.add_argument('-r2', '--mate2', required=True, help='second RNAseq mate')
	parser.add_argument('-ct', '--contigthreshold', required=True, help='threshold on contig length for blast')
	parser.add_argument('-sr', '--refstar', required=True, help='STAR host reference')
	parser.add_argument('-br', '--refbowtie', required=True, help='bowtie2 host reference')
	parser.add_argument('-db', '--db', required=True, help='blast (nt) database')
	parser.add_argument('-s', '--steps', default='12345', help='steps to run (default: 2345 - i.e, steps 2 through 5')
	parser.add_argument("--noclean", action="store_true", help="do not delete temporary intermediate files (default: off)")
	parser.add_argument('--scripts', default=software, help='location of scripts dir (directory where this script resides - use this option only if qsub-ing with the Oracle Grid Engine)')
	parser.add_argument('--verbose', action='store_true', help='verbose mode: echo commands, etc (default: off)')

	args = parser.parse_args()

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
	check_error(args)

	# start with job id set to zero
	jid = 0

	# dict which maps each step to 2-tuple, which contains the qsub part of the command,
	# and the shell part of the command
	d = 	{
			'1': ('qsub -N Pan_tcr', '{}/resources/tcr_repertoire.sh {} {} {}'.format(args.scripts, args.mate1, args.mate2, args.scripts)),
			'2': ('qsub -N Pan_hsep', '{}/resources/host_separation.sh {} {} {} {} {}'.format(args.scripts, args.mate1, args.mate2, args.refstar, args.refbowtie, int(args.noclean))),
			'3': ('qsub -N Pan_asm', '{}/resources/assembly.sh'.format(args.scripts)),
			'4': ('qsub -N Pan_blst', '{}/resources/blast_contigs.sh {} {}'.format(args.scripts, args.contigthreshold, args.db)),
			'5': ('qsub -N Pan_orf', '{}/resources/orf_discovery.sh'.format(args.scripts))
		}

	# run steps
	for i in map(str, range(1,6)): 
		# if step requested
		if i in args.steps:
			# define qsub part of command
			cmd = d[i][0] + '_' + args.identifier + ' '
			# if not the first command, hold on previous job id
			if jid > 0: cmd += '-hold_jid ' + jid + ' '
			# define shell (non-qsub) part of command
			cmd += d[i][1]
			# run command, get job id
			if args.verbose: print(cmd)
			jid = getjid(subprocess.check_output(cmd, shell=True))
			print('Step ' + i + ', jid = ' + jid)

# -------------------------------------

def check_error(args):
	'''Check for errors, check dependencies '''

	# check for required programs
	helpers.check_dependencies(['samtools', 'bam', 'bowtie2', 'STAR', 'findorf', 'blastn', 'Trinity', 'prodigal'])

# -------------------------------------

if __name__ == '__main__':

	main()
