#!/usr/bin/env python

import sys, os

# convert sam file to fastq
# for this sam file, these reads are unmapped, 
# so no need to worry about complementing or multiple mapping

# argument is output fastq base name
fastqbasename = sys.argv[1]

# recapitulate:
# first mate
# samtools view -f 64 $1 | cut -f1,10,11 | awk 'BEGIN{OFS="\t";}{print "@"$1"/1"; print $2; print "+"; print $3}' > ${2}_1.fastq
# second mate
# samtools view -f 128 $1 | cut -f1,10,11 | awk 'BEGIN{OFS="\t";}{print "@"$1"/2"; print $2; print "+"; print $3}' > ${2}_2.fastq

# key: (https://samtools.github.io/hts-specs/SAMv1.pdf)
# 1 0x1 template having multiple segments in sequencing
# 2 0x2 each segment properly aligned according to the aligner
# 4 0x4 segment unmapped
# 8 0x8 next segment in the template unmapped
# 16 0x10 SEQ being reverse complemented
# 32 0x20 SEQ of the next segment in the template being reverse complemented
# 64 0x40 the first segment in the template
# 128 0x80 the last segment in the template
# 256 0x100 secondary alignment
# 512 0x200 not passing filters, such as platform/vendor quality controls
# 1024 0x400 PCR or optical duplicate
# 2048 0x800 supplementary alignment

fastq1 = fastqbasename + '_1.fastq'
fastq2 = fastqbasename + '_2.fastq'

with open(fastq1, 'w') as f:
    with open(fastq2, 'w') as g:
        # loop thro std:in
        for line in sys.stdin:
            id = line.split()[0]
            flag = line.split()[1]
            myread = line.split()[9]
            qual = line.split()[10]
        
            # here's the idea:
            # get the binary val, reverse it, and look at the 6th bit (zero based counting) (refer to key!)
            # >>> x=77
            # >>> bin(x)
            # '0b1001101'
            # >>> bin(x)[2:]
            # '1001101'
            # >>> bin(x)[2:][::-1]
            # '1011001'
            # >>> bin(77)[2:][::-1][6]
            # '1'
        
            # if mate 1 read
            if (bin(int(flag))[2:][::-1][6] == '1'):
                f.write('@' + id + '/1\n' + myread + '\n' + '+\n' + qual + '\n')
            # else if mate 2 read
            else:
                g.write('@' + id + '/2\n' + myread + '\n' + '+\n' + qual + '\n')
