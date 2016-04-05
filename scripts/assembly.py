#!/usr/bin/env python
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=12G,time=12::
#$ -pe smp 8 -R y

import argparse
import sys

# This script performs assembly on the reads leftover after host removal

# -------------------------------------

def get_arg():
    """Get Arguments
    :rtype: object
    """
    # parse arguments

    prog_description = 'Trinity assembly'
    parser = argparse.ArgumentParser(description=prog_description)
    parser.add_argument('-1', '--mate1', default='host_separation/unmapped_1.fastq.gz', help='mate1')
    parser.add_argument('-2', '--mate2', default='host_separation/unmapped_2.fastq.gz', help='mate2')
    parser.add_argument('-o', '--outputdir', default='assembly_trinity', help='the output directory')
    parser.add_argument('--trinitymem', default='50G', help='max memory for Trinity')
    parser.add_argument('--trinitycores', default='8', help='number of cores for Trinity')
    parser.add_argument('-l', '--logsdir', help='the logs directory')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('--noclean', help='do not delete temporary intermediate files (default: off)')
    parser.add_argument('--verbose', action='store_true', help='verbose mode: echo commands, etc (default: off)')
    args = parser.parse_args()

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # print args
    print(args)
    print

    # error checking: exit if previous step produced zero output
    for i in [args.mate1, args.mate2]:
        hp.check_file_exists_and_nonzero(i)

    return args

# -------------------------------------

def assembly(args):
    """Do Trinity assembly"""

    # mkdir -p
    hp.mkdirp(args.outputdir)

    # perform Trinity assembly
    print('ASSEMBLY START')
    cmd = 'Trinity --seqType fq --normalize_reads --max_memory {} --CPU {}  --output {} --left {} --right {}'.format(
        args.trinitymem,
        args.trinitycores,
        args.outputdir,
        args.mate1,
        args.mate2
    )
    hp.run_cmd(cmd, args.verbose, 0)
    print('ASSEMBLY END')

    # name of expected output file
    myoutput = args.outputdir + '/Trinity.fasta'

    # exit if no output
    hp.check_file_exists_and_nonzero(myoutput)

    # rename Trinity contigs, join sequence portion of fasta, return number of contigs
    # cat ${outputdir}/Trinity.fasta | awk 'BEGIN{f=0; counter=1}{if ($0~/^>/) {if (f) {printf "\n"; counter++}; print ">contig_"counter; f=1} else printf $0}END{printf "\n"}' > ${output}
    myoutput2 = 'assembly/contigs_trinity.fasta'
    num_contigs = hp.fastajoinlines(myoutput, myoutput2, 'contig')

    # compute simple distribution (TO-DO)
    # cat ${output} | paste - - | awk '{print length($2)}' | sort -nr | ${d}/scripts/tablecount | awk -v tot=${num_contigs} 'BEGIN{x=0}{x+=$2; print $1"\t"$2"\t"x"/"tot"\t"int(100*x/tot)"%"}' > assembly/contigs.distrib.txt

    if not int(args.noclean):
        cmd = 'rm -r assembly_trinity'
        hp.run_cmd(cmd, args.verbose, 0)

    # return the name of assembly file
    return myoutput2

# -------------------------------------

def remap(args, contigs):
    """Map contigs back onto assembly"""

    print('REMAP START')

    hp.mkdirp('assembly/ref_remap')

    refbowtie="assembly/ref_remap/ref"

    cmd = 'bowtie2-build {} {}'.format(contigs, refbowtie)
    hp.run_cmd(cmd, args.verbose, 0)

    cmd = 'bowtie2 -p 4 -x {} -1 {} -2 {} -S {}'.format(refbowtie, args.mate1, args.mate2, 'assembly/reads2contigs.sam')
    hp.run_cmd(cmd, args.verbose, 0)

    # convert to bam
    # BAM index stats
    # mpileup
    cmd = 'samtools view -bS assembly/reads2contigs.sam | samtools sort - assembly/reads2contigs; \
        samtools index assembly/reads2contigs.bam; \
        rm assembly/reads2contigs.sam; \
        samtools idxstats assembly/reads2contigs.bam > assembly/reads2contigs.stats.txt; \
        samtools mpileup -A -B -d 10000 -L 10000 -f ${output} assembly/reads2contigs.bam > assembly/reads2contigs.pileup'
    hp.run_cmd(cmd, args.verbose, 0)

    if not int(args.noclean):
        cmd = 'rm -r assembly/ref_remap'
        hp.run_cmd(cmd, args.verbose, 0)

    print('REMAP END')

# -------------------------------------

def main():
    """Main function"""

    # get arguments
    args = get_arg()
    # assembly
    contigs = assembly(args)
    # remap
    remap(args, contigs)

# -------------------------------------

if __name__ == '__main__':

    main()
