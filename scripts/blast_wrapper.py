#!/usr/bin/env python
#$ -V
#$ -cwd
#$ -o log.out
#$ -e log.err
#$ -l mem=4G,time=8::

import argparse
import sys
import subprocess

# This script blasts the entries of a fasta file,
# provided they're over the length threshold

# -------------------------------------

def get_arg():
    """Get Arguments
    :rtype: object
    """
    # parse arguments

    prog_description = 'Blast in parallel'
    parser = argparse.ArgumentParser(description=prog_description)
    parser.add_argument('-i', '--input', default='assembly/contigs_trinity.fasta', help='the input fasta')
    parser.add_argument('-o', '--outputdir', default='blast', help='the output directory')
    parser.add_argument('-l', '--logsdir', default='logs_blast', help='the logs directory')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('--noclean', type=int, default=0, help='do not delete temporary intermediate files (default: off)')
    parser.add_argument('--verbose', action='store_true', help='verbose mode: echo commands, etc (default: off)')
    parser.add_argument('--threshold', type=int, default=0, help='the length threshold')
    parser.add_argument('--db', help='the database prefix')
    parser.add_argument('--whichblast', default='blastn', choices=['blastn', 'blastp'], help='which blast to use (blastn, blastp)')
    parser.add_argument('--nosge', type=int, default=0, help='no SGE bool')
    parser.add_argument('--id', help='id')
    args = parser.parse_args()

    # blast format string
    fmt="qseqid sseqid saccver staxids pident nident length mismatch gapopen gaps qstart qend qlen qframe qcovs sstart send slen sframe sstrand evalue bitscore stitle"

    # add key-value pairs to the args dict
    vars(args)['fmt'] = fmt

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # print args
    print(args)
    print

    # error checking: exit if previous step produced zero output
    for i in [args.input]:
        hp.check_file_exists_and_nonzero(i)

    return args

# -------------------------------------

def blast(args):
    """Do blast in parallel"""

    print('------------------------------------------------------------------')
    print('BLAST START')

    # mkdir -p 
    hp.mkdirp(args.outputdir)
    hp.mkdirp(args.logsdir)

    # split fasta file on contigs above threshold length (and return count)
    filecount = hp.fastasplit(args.input, args.outputdir + '/blast', args.threshold)

    if filecount == 0:
        print("No contigs above threshold. Exiting")
        sys.exit(1)

    # generate header
    with open(args.outputdir + '/header', 'w') as f:
        f.write(args.fmt.replace(' ', '\t') + '\n')

    # if no qsub
    if args.nosge:
        for i in range(1,filecount+1):
            # define log files
            logs_out = args.logsdir + '/' + 'bc_' + args.id + '.' + str(i) + '.o'
            logs_err = args.logsdir + '/' + 'bc_' + args.id + '.' + str(i) + '.e'
	    # define command: run blast in series (this will be slow!)
            cmd = '{}/scripts/blast.sh --outputdir {} --whichblast {} --db {} --fmt "{}" --sgeid {} > {} 2> {}'.format(
                      args.scripts,
                      args.outputdir,
                      args.whichblast,
                      args.db,
                      args.fmt,
                      str(i),
                      logs_out,
                      logs_err
            )
            hp.run_cmd(cmd, args.verbose, 0)

        # concatenate results
        cmd = '{}/scripts/concat.sh {} {} {}'.format(
                  args.scripts,
                  args.outputdir,
                  args.logsdir,
                  int(args.noclean)
        )
        hp.run_cmd(cmd, args.verbose, 0)
    else:
        # qsub part of command (array job)
        qcmd = 'qsub -N bc_{} -e {} -o {} -t 1-{} '.format(
                  args.id,
                  args.logsdir,
                  args.logsdir,
                  filecount,
        )
        # regular part of command
        cmd = '{}/scripts/blast.sh --outputdir {} --whichblast {} --db {} --fmt "{}"'.format(
                  args.scripts,
                  args.outputdir,
                  args.whichblast,
                  args.db,
                  args.fmt
        )
        # print qsub message
        message = subprocess.check_output(qcmd + cmd, shell=True)
        print(message)
        # get job id
        jid = hp.getjid(message)

        # hold the script up here, until all the blast jobs finish
        # concat top blast hits; concat log files into one, so as not to clutter the file system
        # qsub part of command
        qcmd = 'qsub -N wait_{} -hold_jid {} -sync y '.format(
                  args.id,
                  jid
        )
        # regular part of command
        cmd = '{}/scripts/concat.sh {} {} {}'.format(
                  args.scripts,
                  args.outputdir,
                  args.logsdir,
                  int(args.noclean)
        )
        message = subprocess.check_output(qcmd + cmd, shell=True)
        print(message)

    print('BLAST END')
    print('------------------------------------------------------------------')

# -------------------------------------

def main():
    """Main function"""

    # get arguments
    args = get_arg()
    # blast
    blast(args)

# -------------------------------------

if __name__ == '__main__':

    main()
