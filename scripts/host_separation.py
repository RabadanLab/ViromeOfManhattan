#!/usr/bin/env python

# This script removes host (human) reads
# using a pass of STAR mapping, followed by
# a pass of bowtie2 mapping

import argparse
import sys

# -------------------------------------

def get_arg():
    """Get Arguments
    :rtype: object
    """
    # parse arguments

    prog_description = 'Separate host reads'
    parser = argparse.ArgumentParser(description=prog_description)
    parser.add_argument('-1', '--mate1', help='mate1 fastq')
    parser.add_argument('-2', '--mate2', help='mate2 fastq')
    parser.add_argument('--bam', help='bam file')
    parser.add_argument('-o', '--outputdir', default='host_separation', help='the output directory')
    parser.add_argument('-l', '--logsdir', help='the logs directory')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('--numthreads', default='4', help='the number of threads with which to run STAR and bowtie')
    parser.add_argument('--refstar', help='the star reference')
    parser.add_argument('--refbowtie', help='the bowtie reference')
    parser.add_argument('--gtf', help='host gene feature gtf')
    parser.add_argument('--noclean', type=int, default=0, help='do not delete temporary intermediate files (default: off)')
    parser.add_argument('--gzip', type=int, default=0, help='input files are gzipped boolean (default: off)')
    parser.add_argument('--verbose', type=int, default=0, help='verbose mode: echo commands, etc (default: off)')
    args = parser.parse_args()

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # add key-value pairs to the args dict
    vars(args)['step'] = 'host_separation'
    # vars(args)['olog'] = args.outputdir + '/../' + 'log.hostmap.out'
    # vars(args)['elog'] = args.outputdir + '/../' + 'log.hostmap.err'
    vars(args)['olog'] = args.outputdir + '/../' + 'log.out'
    vars(args)['elog'] = args.outputdir + '/../' + 'log.err'

    # error checking: exit if input empty 
    for i in [args.mate1, args.mate2]:
        if i != 'None':
            hp.check_file_exists_and_nonzero(i, step=args.step)

    return args

# -------------------------------------

def hostsep(args):
    """Separate host reads"""

    print('Counting input reads')
    cmd = 'wc -l {args.mate1} > {args.outputdir}/mapping_percent.txt'.format(args=args)
    hp.run_cmd(cmd, args.verbose, 0)

    # flags for STAR
    starflag=''
    # if input files are gzipped
    if args.gzip: 
        starflag='--readFilesCommand zcat'

    print('STAR mapping commenced')
    cmd = 'STAR --runThreadN {args.numthreads} --genomeDir {args.refstar} --readFilesIn {args.mate1} {args.mate2} --outFileNamePrefix {args.outputdir}/ --outSAMtype BAM Unsorted --outSAMunmapped Within {starflag}'.format(args=args, starflag=starflag)
    hp.run_cmd(cmd, args.verbose, 0)
    print('STAR mapping finished')

    print('find unmapped reads')

    cmd = 'samtools flagstat {args.outputdir}/Aligned.out.bam > {args.outputdir}/mapping_stats.txt'.format(args=args)
    hp.run_cmd(cmd, args.verbose, 0)

    # bin(13) = '0b1101', which corresponds to SAM flag bits:
    # read paired; read unmapped; mate unmapped
    cmd = 'samtools view -b -f 13 {args.outputdir}/Aligned.out.bam | samtools sort -n - {args.outputdir}/star_unmapped'.format(args=args)
    hp.run_cmd(cmd, args.verbose, 0)

    cmd = 'samtools view {args.outputdir}/star_unmapped.bam | {args.scripts}/scripts/sam2fastq.py {args.outputdir}/star_unmapped'.format(args=args)
    hp.run_cmd(cmd, args.verbose, 0)

    cmd = 'wc -l {args.outputdir}/star_unmapped_1.fastq >> {args.outputdir}/mapping_percent.txt'.format(args=args)
    hp.run_cmd(cmd, args.verbose, 0)

    print('Bowtie2 mapping commenced')
    cmd = 'bowtie2 -p {args.numthreads} -x {args.refbowtie} -1 {args.outputdir}/star_unmapped_1.fastq -2 {args.outputdir}/star_unmapped_2.fastq -S {args.outputdir}/bwt2.sam'.format(args=args)
    # hp.run_cmd(cmd, args.verbose, 0)
    hp.run_log_cmd(cmd, args.verbose, args.olog, args.elog)
    print('Bowtie2 mapping finished')

    print('find unmapped reads')

    cmd = 'samtools view -S -b -f 13 {args.outputdir}/bwt2.sam | samtools sort -n - {args.outputdir}/bwt2_unmapped'.format(args=args)
    hp.run_cmd(cmd, args.verbose, 0)

    cmd = 'samtools view {args.outputdir}/bwt2_unmapped.bam | {args.scripts}/scripts/sam2fastq.py {args.outputdir}/bwt2_unmapped'.format(args=args)
    hp.run_cmd(cmd, args.verbose, 0)

    # if gtf variable set, get gene coverage
    if args.gtf:
        print('featureCounts commenced')
        cmd = 'featureCounts -a {args.gtf} -o {args.outputdir}/host_gene_counts.txt {args.outputdir}/Aligned.out.bam'.format(args=args)
        # hp.run_cmd(cmd, args.verbose, 0)
        hp.run_log_cmd(cmd, args.verbose, args.olog, args.elog)
        print('featureCounts finished')

    # TO DO: make this code more compact
    if not args.noclean:
        print('clean up')
        cmd = 'rm -rf ' + args.outputdir + '/' + '_STARtmp'
        hp.run_cmd(cmd, args.verbose, 0)
        for i in ['Aligned.out.bam', 'Log.*', 'SJ.out.tab', 'star_unmapped.bam', 'star_unmapped_*.fastq', 'bwt2.sam', 'bwt2_unmapped.bam']:
            cmd = 'rm {args.outputdir}/{i}'.format(args=args, i=i)
            hp.run_cmd(cmd, args.verbose, 0)

    # rename and zip both mates
    for i in ['1', '2']:
        cmd = 'mv {args.outputdir}/bwt2_unmapped_{i}.fastq {args.outputdir}/unmapped_{i}.fastq'.format(args=args, i=i)
        hp.run_cmd(cmd, args.verbose, 0)

        cmd = 'gzip {args.outputdir}/unmapped_{i}.fastq'.format(args=args, i=i)
        hp.run_cmd(cmd, args.verbose, 0)

    cmd = 'wc -l {args.outputdir}/unmapped_1.fastq.gz >> {args.outputdir}/mapping_percent.txt'.format(args=args)
    hp.run_cmd(cmd, args.verbose, 0)

    hp.echostep(args.step, start=0)

# -------------------------------------

def getunmapped(args):
    """Starting with a .bam file, get the unmapped reads"""

    # fix violations of DRY (modify args variable)

    print('find unmapped reads')

    cmd = 'samtools view -b -f 13 {args.bam} | samtools sort -n - {args.outputdir}/unmapped'.format(args=args)
    hp.run_cmd(cmd, args.verbose, 0)

    cmd = 'samtools view {args.outputdir}/unmapped.bam | {args.scripts}/scripts/sam2fastq.py {args.outputdir}/unmapped'.format(args=args)
    hp.run_cmd(cmd, args.verbose, 0)

    # if gtf variable set, get gene coverage
    if args.gtf:
        print('featureCounts commenced')
        cmd = 'featureCounts -a {args.gtf} -o {args.outputdir}/host_gene_counts.txt {args.bam}'.format(args=args)
        hp.run_log_cmd(cmd, args.verbose, args.olog, args.elog)
        print('featureCounts finished')

    # zip both mates
    for i in ['1', '2']:
        cmd = 'gzip {args.outputdir}/unmapped_{i}.fastq'.format(args=args, i=i)
        hp.run_cmd(cmd, args.verbose, 0)

    if not args.noclean:
        print('clean up')
        cmd = 'rm -rf ' + args.outputdir + '/' + 'unmapped.bam'
        hp.run_cmd(cmd, args.verbose, 0)

    hp.echostep(args.step, start=0)

# -------------------------------------

def main():
    """Main function"""

    # get arguments
    args = get_arg()

    hp.echostep(args.step)

    # print args
    print(args)
    print

    # mkdir -p
    hp.mkdirp(args.outputdir)

    if args.bam and args.bam != 'None':
        # get unmapped reads
        getunmapped(args)
    else:
        # host separation
        hostsep(args)

# -------------------------------------

if __name__ == '__main__':

    main()
