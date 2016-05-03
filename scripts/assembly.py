#!/usr/bin/env python

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
    parser.add_argument('--verbose', type=int, default=0, help='verbose mode: echo commands, etc (default: off)')
    args = parser.parse_args()

    # add key-value pairs to the args dict
    vars(args)['step'] = 'assembly'

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # error checking: exit if previous step produced zero output
    for i in [args.mate1, args.mate2]:
        hp.check_file_exists_and_nonzero(i, step=args.step)

    return args

# -------------------------------------

def assembly(args):
    """Do Trinity assembly"""

    hp.echostep(args.step)

    # print args
    print(args)
    print

    # mkdir -p assembly_trinity
    hp.mkdirp(args.outputdir)

    # perform Trinity assembly
    cmd = 'Trinity --seqType fq --normalize_reads --max_memory {} --CPU {}  --output {} --left {} --right {}'.format(
        args.trinitymem,
        args.trinitycores,
        args.outputdir,
        args.mate1,
        args.mate2
    )
    # use run_long_cmd for programs with verbose output
    hp.run_long_cmd(cmd, args.verbose, 'log.Trinity')

    print('Trinity complete')

    # name of expected output file
    myoutput = args.outputdir + '/Trinity.fasta'

    # exit if no output
    hp.check_file_exists_and_nonzero(myoutput)

    # mkdir -p
    hp.mkdirp('assembly')

    # rename Trinity contigs, join sequence portion of fasta, return number of contigs
    # cat ${outputdir}/Trinity.fasta | awk 'BEGIN{f=0; counter=1}{if ($0~/^>/) {if (f) {printf "\n"; counter++}; print ">contig_"counter; f=1} else printf $0}END{printf "\n"}' > ${output}
    myoutput2 = 'assembly/contigs_trinity.fasta'
    num_contigs = hp.fastajoinlines(myoutput, myoutput2, 'contig')

    # compute simple distribution
    # cat assembly/contigs_trinity.fasta | paste - - | awk '{print length($2)}' | sort -nr | ${d}/scripts/tablecount | awk -v tot=${num_contigs} 'BEGIN{x=0}{x+=$2; print $1"\t"$2"\t"x"/"tot"\t"int(100*x/tot)"%"}' > assembly/contigs.distrib.txt
    computedistrib(myoutput2, 'assembly/contigs.distrib.txt')

    if not int(args.noclean):
        cmd = 'rm -rf assembly_trinity'
        hp.run_cmd(cmd, args.verbose, 0)

    hp.echostep(args.step, start=0)

    # return the name of assembly file
    return myoutput2

# -------------------------------------

def computedistrib(infile, outfile):
    """compute simple distribution of contigs"""

    # a list of contig lengths
    x = []

    with open(infile, 'r') as g:
        for line in g:
            if line[0] != '>': 
                x.append(len(line.rstrip()))

    from collections import Counter

    # count occurences of each length
    xcount = Counter(x)

    # tot length
    tot = sum(dict(xcount).values())

    # running sum
    runningsum = 0

    with open(outfile, 'w') as f:
        for i in sorted(dict(xcount), reverse=True):
            runningsum += xcount[i]
            f.write(str(i) + '\t')
            f.write(str(xcount[i]) + '\t')
            f.write(str(runningsum) + '/' + str(tot) + '\t')
            f.write(str(int(100*runningsum/tot)) + '%')
            f.write('\n')
        
# -------------------------------------

def remap(args, contigs):
    """map contigs back onto assembly"""

    hp.echostep('remap')

    hp.mkdirp('assembly/ref_remap')

    refbowtie="assembly/ref_remap/ref"

    cmd = 'bowtie2-build {} {}'.format(contigs, refbowtie)
    hp.run_cmd(cmd, args.verbose, 0)

    cmd = 'bowtie2 -p 4 -x {} -1 {} -2 {} -S {}'.format(refbowtie, args.mate1, args.mate2, 'assembly/reads2contigs.sam')
    hp.run_cmd(cmd, args.verbose, 0)

    # convert to bam
    cmd = 'samtools view -bS assembly/reads2contigs.sam | samtools sort - assembly/reads2contigs'
    hp.run_cmd(cmd, args.verbose, 0)

    cmd = 'samtools index assembly/reads2contigs.bam'
    hp.run_cmd(cmd, args.verbose, 0)

    cmd = 'rm assembly/reads2contigs.sam'
    hp.run_cmd(cmd, args.verbose, 0)

    # BAM index stats
    cmd = 'samtools idxstats assembly/reads2contigs.bam > assembly/reads2contigs.stats.txt'
    hp.run_cmd(cmd, args.verbose, 0)

    # mpileup
    cmd = 'samtools mpileup -A -B -d 100000 -L 100000 -f assembly/contigs_trinity.fasta assembly/reads2contigs.bam > assembly/reads2contigs.pileup'
    hp.run_cmd(cmd, args.verbose, 0)

    # format pileup file - i.e., add zeros to uncovered positions
    formatpileup('assembly/reads2contigs.pileup', 'assembly/reads2contigs.stats.txt', 'assembly/reads2contigs.format.pileup')

    if not int(args.noclean):
        cmd = 'rm -r assembly/ref_remap'
        hp.run_cmd(cmd, args.verbose, 0)
        cmd = 'rm assembly/reads2contigs.pileup'
        hp.run_cmd(cmd, args.verbose, 0)

    hp.echostep('remap', start=0)

# -------------------------------------

def formatpileup(infile, idxfile, outfile):
    """format the pileup file for computing entropy"""

    # id 2 length dict (output of samtools idxstats)
    idx = {}

    # load idx file
    with open(idxfile, 'r') as f:
        for line in f:
            # map id to length of contig
            idx[line.split()[0].strip()] = line.split()[1].strip()

    myid = ''    # contig id
    pos = ''    # position

    with open(outfile, 'w') as f:
        with open(infile, 'r') as g:
            for line in g:
                # get number reads 
                numrds = line.split()[3]

                # if beginning of a new contig (id != previous id)
                if line.split()[0] != myid:
                    # if change (and not first contig), check if previous contig was covered until the end
                    # if not covered, pad with zeros
                    if myid: 
                        if int(idx[myid]) > int(pos): 
                            for i in range(int(pos) + 1, int(idx[myid]) + 1): 
                                f.write(myid + '\t' + str(i) + '\t0\n')
 
                    # if contig starts at postion > 1, pad with zeros 
                    if int(line.split()[1]) > 1: 
                        for i in range(1, int(line.split()[1])): 
                            f.write(line.split()[0] + '\t' + str(i) + '\t0\n')

                    # set new id 
                    myid = line.split()[0]

                    # write current line
                    f.write(myid + '\t' + line.split()[1] + '\t' + numrds + '\n')

                # if discontinuity (position - previous position > 1), pad with zeros
                elif (int(line.split()[1]) - int(pos)) > 1:
                    for i in range(int(pos) + 1, int(line.split()[1])):
                        f.write(myid + '\t' + str(i) + '\t0\n')

                    f.write(myid + '\t' + line.split()[1] + '\t' + numrds + '\n')

                # otherwise, simply write line
                else:
                    f.write(myid + '\t' + line.split()[1] + '\t' + numrds + '\n')

                # get position (this will become previous position for next iteration)
                pos = line.split()[1]

    # check if last contig covered until the end
    if int(idx[myid]) > int(pos):
        for i in range(int(pos) + 1, int(idx[myid]) + 1):
            f.write(myid + '\t' + str(i) + '\t0\n')

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
