#!/usr/bin/env python
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=16G,time=12::
#$ -pe smp 4 -R y

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
    parser.add_argument('-1', '--mate1', help='mate1')
    parser.add_argument('-2', '--mate2', help='mate2')
    parser.add_argument('-o', '--outputdir', default='host_separation', help='the output directory')
    parser.add_argument('-l', '--logsdir', help='the logs directory')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('--numthreads', default='4', help='the number of threads with which to run STAR and bowtie')
    parser.add_argument('--refstar', help='the star reference')
    parser.add_argument('--refbowtie', help='the bowtie reference')
    parser.add_argument('--gtf', help='host gene feature gtf')
    parser.add_argument('--noclean', type=int, default=0, help='do not delete temporary intermediate files (default: off)')
    parser.add_argument('--gzip', type=int, default=0, help='input files are gzipped boolean (default: off)')
    parser.add_argument('--verbose', action='store_true', help='verbose mode: echo commands, etc (default: off)')
    args = parser.parse_args()

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # add key-value pairs to the args dict
    # vars(args)['olog'] = args.outputdir + '/../' + 'log.hostmap.out'
    # vars(args)['elog'] = args.outputdir + '/../' + 'log.hostmap.err'
    vars(args)['olog'] = args.outputdir + '/../' + 'log.out'
    vars(args)['elog'] = args.outputdir + '/../' + 'log.err'

    # print args
    print(args)
    print

    # error checking: exit if input empty 
    for i in [args.mate1, args.mate2]:
        hp.check_file_exists_and_nonzero(i)

    return args

# -------------------------------------

def hostsep(args):
    """Separate host reads"""

    print('------------------------------------------------------------------')
    print('HOST_SEPARATION START')

    # mkdir -p 
    hp.mkdirp(args.outputdir)

    # flags for STAR
    starflag=''
    # if input files are gzipped
    if args.gzip: 
        starflag='--readFilesCommand zcat'

    print('STAR mapping commenced')
    cmd = 'STAR --runThreadN {} --genomeDir {} --readFilesIn {} {} --outFileNamePrefix {} --outSAMtype {} --outSAMunmapped {} {}'.format(
        args.numthreads,
        args.refstar,
        args.mate1,
        args.mate2,
        args.outputdir + '/',
        'BAM Unsorted',
        'Within',
        starflag,
    )
    hp.run_cmd(cmd, args.verbose, 0)
    print('STAR mapping finished')

    print('find unmapped reads')

    cmd = 'samtools flagstat {}/Aligned.out.bam > {}/mapping_stats.txt'.format(
        args.outputdir,
        args.outputdir
    )
    hp.run_cmd(cmd, args.verbose, 0)

    cmd = 'samtools view -b -f 13 {}/Aligned.out.bam | samtools sort -n - {}/star_unmapped'.format(
        args.outputdir,
        args.outputdir
    )
    hp.run_cmd(cmd, args.verbose, 0)

    cmd = 'samtools view {}/star_unmapped.bam | {}/scripts/sam2fastq.py {}/star_unmapped'.format(
        args.outputdir,
        args.scripts,
        args.outputdir
    )
    hp.run_cmd(cmd, args.verbose, 0)

    print('Bowtie2 mapping commenced')
    cmd = 'bowtie2 -p {} -x {} -1 {} -2 {} -S {}'.format(
        args.numthreads,
        args.refbowtie,
        args.outputdir + '/' + 'star_unmapped_1.fastq',
        args.outputdir + '/' + 'star_unmapped_2.fastq',
        args.outputdir + '/' + 'bwt2.sam'
    )
    # hp.run_cmd(cmd, args.verbose, 0)
    hp.run_log_cmd(cmd, args.verbose, args.olog, args.elog)
    print('Bowtie2 mapping finished')

    print('find unmapped reads')

    cmd = 'samtools view -S -b -f 13 {}/bwt2.sam | samtools sort -n - {}/bwt2_unmapped'.format(
        args.outputdir,
        args.outputdir
    )
    hp.run_cmd(cmd, args.verbose, 0)

    cmd = 'samtools view {}/bwt2_unmapped.bam | {}/scripts/sam2fastq.py {}/bwt2_unmapped'.format(
        args.outputdir,
        args.scripts,
        args.outputdir
    )
    hp.run_cmd(cmd, args.verbose, 0)

    # if gtf variable set, get gene coverage
    if args.gtf:
        print('featureCounts commenced')
        cmd = 'featureCounts -a {} -o {}/host_gene_counts.txt {}/Aligned.out.bam'.format(
            args.gtf,
            args.outputdir,
            args.outputdir
        )
        # hp.run_cmd(cmd, args.verbose, 0)
        hp.run_log_cmd(cmd, args.verbose, args.olog, args.elog)
        print('featureCounts finished')

    # TO DO: make this code more compact
    if not args.noclean:
        print('clean up')
        cmd = 'rm -rf ' + args.outputdir + '/' + '_STARtmp'
        hp.run_cmd(cmd, args.verbose, 0)
        for i in ['Aligned.out.bam', 'Log.*', 'SJ.out.tab', 'star_unmapped.bam', 'star_unmapped_*.fastq', 'bwt2.sam', 'bwt2_unmapped.bam']:
            cmd = 'rm ' + args.outputdir + '/' + i
            hp.run_cmd(cmd, args.verbose, 0)

    # rename and zip both mates
    for i in ['1', '2']: 
        cmd = 'mv {}/bwt2_unmapped_{}.fastq {}/unmapped_{}.fastq'.format(
            args.outputdir,
            i,
            args.outputdir,
            i
        )
        hp.run_cmd(cmd, args.verbose, 0)

        cmd = 'gzip {}/unmapped_{}.fastq'.format(
            args.outputdir,
            i
        )
        hp.run_cmd(cmd, args.verbose, 0)

    print('HOST_SEPARATION END')
    print('------------------------------------------------------------------')

# -------------------------------------

def main():
    """Main function"""

    # get arguments
    args = get_arg()
    # host separation
    hostsep(args)

# -------------------------------------

if __name__ == '__main__':

    main()
