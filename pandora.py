import argparse
import subprocess

parser = argparse.ArgumentParser(description='microbial detection from paired-end RNAseq')
parser.add_argument('-id', '--identifier', required=True,
        help='5 chars or less sample ID')
parser.add_argument('-r1', '--mate1', required=True,
        help='first RNAseq mate')
parser.add_argument('-r2', '--mate2', required=True,
        help='second RNAseq mate')
parser.add_argument('-ct', '--contigthreshold', required=True,
        help='threshold on contig length for blast')
args = parser.parse_args()

jid_tcr = subprocess.check_output('qsub -N tcr_{0} $APPS/Pandora/tcr_repertoire.sh {1} {2}'.format(args.identifier, args.mate1, args.mate2), shell=True).split('Your job ')[1].split()[0]
print 'jid_tcr = {0}'.format(jid_tcr)

jid_hsep = subprocess.check_output('qsub -hold_jid {0} -N hsep_{1} $APPS/Pandora/host_separation.sh {2} {3}'.format(jid_tcr, args.identifier, args.mate1, args.mate2), shell=True).split('Your job ')[1].split()[0]
print 'jid_hsep = {0}'.format(jid_hsep)

jid_asm = subprocess.check_output('qsub -hold_jid {0} -N asm_{1} $APPS/Pandora/assembly.sh'.format(jid_hsep, args.identifier), shell=True).split('Your job ')[1].split()[0]
print 'jid_asm = {0}'.format(jid_asm)

jid_blst = subprocess.check_output('qsub -hold_jid {0} -N blst_{1} $APPS/Pandora/blast_contigs.sh {2}'.format(jid_asm, args.identifier, args.contigthreshold), shell=True).split('Your job ')[1].split()[0]
print 'jid_blst = {0}'.format(jid_blst)

jid_orf = subprocess.check_output('qsub -hold_jid {0} -N orf_{1} $APPS/Pandora/orf_discovery.sh'.format(jid_blst, args.identifier), shell=True).split('Your job ')[1].split()[0]
print 'jid_orf = {0}'.format(jid_orf)
