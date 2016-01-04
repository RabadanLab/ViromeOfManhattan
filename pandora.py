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
	# string looks like this: Your job 8379811 ("test") has been submitted 
	return x.split('Your job ')[1].split()[0]

# -------------------------------------

def main():
	'''Run all steps'''

	# get arguments dict
	args = get_arg()
	
	# check for errors
	check_error(args)

	# run steps
	jid_tcr = getjid(subprocess.check_output('qsub -N tcr_{0} {1}/resources/tcr_repertoire.sh {2} {3}'.format(args.identifier, args.scripts, args.mate1, args.mate2), shell=True))
	print 'jid_tcr = {0}'.format(jid_tcr)
	
	jid_hsep = getjid(subprocess.check_output('qsub -hold_jid {0} -N hsep_{1} {2}/resources/host_separation.sh {3} {4} {5} {6}'.format(jid_tcr, args.identifier, args.scripts, args.mate1, args.mate2, args.refstar, args.refbowtie), shell=True))
	print 'jid_hsep = {0}'.format(jid_hsep)
	
	jid_asm = getjid(subprocess.check_output('qsub -hold_jid {0} -N asm_{1} {2}/resources/assembly.sh'.format(jid_hsep, args.identifier, args.scripts), shell=True))
	print 'jid_asm = {0}'.format(jid_asm)
	
	jid_blst = getjid(subprocess.check_output('qsub -hold_jid {0} -N blst_{1} {2}/resources/blast_contigs.sh {3} {4}'.format(jid_asm, args.identifier, args.scripts, args.contigthreshold, args.db), shell=True))
	print 'jid_blst = {0}'.format(jid_blst)
	
	jid_orf = getjid(subprocess.check_output('qsub -hold_jid {0} -N orf_{1} {2}/resources/orf_discovery.sh'.format(jid_blst, args.identifier, args.scripts), shell=True))
	print 'jid_orf = {0}'.format(jid_orf)

# -------------------------------------

def check_error(args):
	'''Check for errors, check dependencies '''

	# check for required programs
	helpers.check_dependencies(['samtools', 'bam', 'bowtie2', 'STAR', 'findorf', 'blastn', 'Trinity', 'prodigal'])

# -------------------------------------

if __name__ == '__main__':

	main()
