#!/usr/bin/env python

import argparse
from Bio.Seq import Seq

# Find ORFs in a fasta file, borrowing from:
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc292

# define amino acid table (see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
# 1. The Standard Code
# 2. The Vertebrate Mitochondrial Code
# 3. The Yeast Mitochondrial Code
# 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
# etc...
table = 1 # standard code

prog_description = 'find ORFs'
parser = argparse.ArgumentParser(description=prog_description)
parser.add_argument('-i', '--input', required=True, help='input path')
parser.add_argument('-t', '--threshold', required=True, help='threshold')
# parser.add_argument('-o', '--output', required=True, help='output path')
args = parser.parse_args()

# important: assume the sequence portion of the fasta is all on one line

with open(args.input, 'r') as f:
    # fasta id
    id = ''
    # ORF counter for a given seq id
    counter = 1
    for line in f:
        # parse fasta
	# grab id
	if line[0] == '>':
            id = line.rstrip()[1:].replace(' ', '_')
	    # reset counter
            counter = 1
        # create biopython seq obj
        myseq = Seq(line.rstrip())
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc292
        for strand, nuc in [(+1, myseq), (-1, myseq.reverse_complement())]:
            #print('strand: ' + str(strand))
            #print('nuc: ' + nuc)
            for frame in range(3):
                # the purpose of this line is get a multiple of three
                # e.g., 3/3 * 3 = 3 and 4/3 * 3 = 3 
                length = 3 * ((len(myseq)-frame) // 3)
	        #print('frame: ' + str(frame))
	        #print('len: ' + str(length))
	        #print('nt: '),
                #print(nuc[frame:frame+length])
	        #print('AA: '),
                #print(nuc[frame:frame+length].translate(table))
	        # split on the stop codon, a '*' character
	        # don't worry about start codons
	        for pro in nuc[frame:frame+length].translate(table).split('*'):
                    if len(pro) >= int(args.threshold):
	                # print('>' + id)
                        # print("%s...%s - length %i, strand %i, frame %i" % (pro[:30], pro[-3:], len(pro), strand, frame))
			# print fasta entry
                        print('>%s_ORF%i_len%i_strand%i_frame%i' % (id, counter, len(pro), strand, frame))
                        print(pro)
			counter += 1
