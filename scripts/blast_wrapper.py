#!/usr/bin/env python

import argparse
import sys
import subprocess
import os
import glob
import shutil

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
    parser.add_argument('--verbose', type=int, default=0, help='verbose mode: echo commands, etc (default: off)')
    parser.add_argument('--threshold', type=int, default=0, help='the length threshold')
    parser.add_argument('--db', help='the database prefix')
    parser.add_argument('--whichblast', default='blastn', choices=['blastn', 'blastp'], help='which blast to use (blastn, blastp)')
    parser.add_argument('--threads', default='1', help='blast -num_threads option')
    parser.add_argument('--nosge', type=int, default=0, help='no SGE bool')
    parser.add_argument('--id', help='id')
    args = parser.parse_args()

    # blast format string
    fmt="qseqid sseqid saccver staxids pident nident length mismatch gapopen gaps qstart qend qlen qframe qcovs sstart send slen sframe sstrand evalue bitscore stitle"

    # add key-value pairs to the args dict
    vars(args)['fmt'] = fmt
    vars(args)['step'] = 'blast_contigs'

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # error checking: exit if previous step produced zero output
    for i in [args.input]:
        hp.check_file_exists_and_nonzero(i, step=args.step)

    return args

# -------------------------------------

def blast(args):
    """Do blast in parallel"""

    hp.echostep(args.step)

    # print args
    print(args)
    print

    # mkdir -p 
    hp.mkdirp(args.outputdir)

    args.noclean = int(args.noclean)

    # generate header
    with open(args.outputdir + '/header', 'w') as f:
        f.write(args.fmt.replace(' ', '\t') + '\n')

    # if no qsub
    if args.nosge:
        # filter fasta file on contigs above threshold length and hardcode name for blast.py
        # (splitting files doesnt make sense if no cluster)
        filecount = hp.fastafilter(args.input, args.outputdir + '/blast_1.fasta', args.threshold)

        if filecount == 0:
            print("No contigs above threshold. Exiting")
            sys.exit(1)

        # define log files
        logs_out = args.outputdir + '/' + 'blast_log_' + args.id + '.o'
        logs_err = args.outputdir + '/' + 'blast_log_' + args.id + '.e'
        # define command
        cmd = '{args.scripts}/scripts/blast.py --scripts {args.scripts} --outputdir {args.outputdir} --whichblast {args.whichblast} --db {args.db} --threads {args.threads} --fmt "{args.fmt}" --sgeid 1 > {o} 2> {e}'.format(args=args, o=logs_out, e=logs_err)
        hp.run_cmd(cmd, args.verbose, 0)

        # make link
        cmd = 'ln -s blast_1.result {args.outputdir}/concat.txt'.format(args=args)
        hp.run_cmd(cmd, args.verbose, 0)

        # get top hits
        hp.tophitsfilter(args.outputdir + '/blast_1.result', args.outputdir + '/top.concat.txt')
        # get fasta file of entries that didn't blast
        hp.getnohits(args.outputdir + '/top.concat.txt', args.outputdir + '/blast_1.fasta', args.outputdir + '/no_blastn.fa')

    # if qsub
    else:
        # mkdir -p
        hp.mkdirp(args.logsdir)

        # split fasta file on contigs above threshold length (and return count)
        filecount = hp.fastasplit(args.input, args.outputdir + '/blast', args.threshold)

        if filecount == 0:
            print("No contigs above threshold. Exiting")
            sys.exit(1)

        # qsub part of command (array job)
        qcmd = 'qsub -S {mypython} -N bc_{args.id} -e {args.logsdir} -o {args.logsdir} -t 1-{filecount} '.format(mypython=sys.executable, args=args, filecount=filecount)
        # regular part of command
        cmd = '{args.scripts}/scripts/blast.py --scripts {args.scripts} --outputdir {args.outputdir} --whichblast {args.whichblast} --db {args.db} --threads {args.threads} --fmt "{args.fmt}"'.format(args=args)
	if args.verbose:
            print(qcmd + cmd)
        message = subprocess.check_output(qcmd + cmd, shell=True)
        print(message)
        # get job id
        jid = hp.getjid(message)

        # hold the script up here, until all the blast jobs finish
        # concat top blast hits; concat log files into one, so as not to clutter the file system
        # qsub part of command
        qcmd = 'qsub -V -b y -cwd -o log.out -e log.err -l mem=1G,time=1:: -N wait_{args.id} -hold_jid {jid} -sync y echo wait_here'.format(args=args, jid=jid)
        message = subprocess.check_output(qcmd, shell=True)
        print(message)

        # now concatenate blast results
        concat(args)

    hp.echostep(args.step, start=0)

# -------------------------------------

def concat(args):
    """concatenate blast files and logs, so as not to leave many files messily scattered about"""

    print('CONCATENATE START')

    # ids that didn't blast
    noblast = []

    # fasta file of entries that didn't blast
    noblastfile = open(args.outputdir + '/no_blastn.fa', 'w')
    # file of top blast hits
    tophitsfile = open(args.outputdir + '/top.concat.txt', 'w')
    # file of all blast hits
    allhitsfile = open(args.outputdir + '/concat.txt', 'w')
    # all fasta entries
    fastafile = open(args.outputdir + '/above_threshold.fa', 'w')

    # glob blast files
    files = glob.glob(args.outputdir + '/*.result')

    for f in files: 
        # get file length
	cmd = 'cat ' + f + ' | wc -l'
        len = subprocess.check_output(cmd, shell=True).strip()

        # get the name of the fasta file (remove last 6 characters 'result')
        g = f[:-6] + 'fasta'

        with open(f, 'r') as h:
            for line in h:
                allhitsfile.write(line)
        with open(g, 'r') as h:
            for line in h:
                fastafile.write(line)

	# if length of file is zero
        if len == '0':
            noblast.append(os.path.basename(f))
            # cat ${base}.fasta | sed s/X//g
            with open(g, 'r') as h:
                for line in h:
                    noblastfile.write(line)
	# else write tophits file
        else:
            with open(f, 'r') as h:
                ## Filtering function based on entire list of blast hits for every given contig (and tax id and alignment parameters)
                
                # Index reference in blast report file:
                accession_index = 2
                tax_id_index = 3
                contig_len_index = 12
                alignment_length_index = 5
                score_index = 21
                name_index = 22
                e_val_index = 20


                skip_blast_hit = False

                top_line = h.readline()
    
                top_contig = top_line.split('\t')[0].strip()

                top_score = float(top_line.split('\t')[score_index].strip())
                
                # Top number of BLAST hits to parse through in order to determine whether top hit can be trusted as truly non-human
                top_number = 10
    
                if top_line == "":
                    skip_blast_hit = True
                    break
    
                else:
                    count = 0
                    while count <= top_number:
                        next_line = h.readline()
            
                        if next_line == "" or next_line.split('\t')[0].strip() != top_contig:
                            break
                        else:
                            # condition to eliminate contig if accurate mapping to human genome
                            row_values = next_line.split('\t')

                            row_score = float(row_values[score_index].strip())

                            if (row_values[tax_id_index].strip() == '9606') &  ( float(row_values[e_val_index].strip()) < 10**(-4) ) & (row_score > 0.95 * top_score):
                                skip_blast_hit = True
                                break
                
                # Only include the contig in the top concat.txt report of top blast hits if it is NOT a suspected high-alignment score human read
                if not skip_blast_hit:
                    tophitsfile.write(top_line)

        if not args.noclean:
            os.remove(f)
            os.remove(g)

    noblastfile.close()
    tophitsfile.close()
    allhitsfile.close()
    fastafile.close()

    print('No blast hits for: ' + ' '.join(noblast))

    # concat blast logs and remove folder
    print('concatenate blast logs')
    cmd = 'head -100 {args.logsdir}/* > {args.outputdir}/log.blast'.format(args=args)
    hp.run_cmd(cmd, args.verbose, 0)

    if not args.noclean:
        shutil.rmtree(args.logsdir)

    print('CONCATENATE END')

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
