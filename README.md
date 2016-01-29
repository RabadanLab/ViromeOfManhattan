Pandora
=======

Identification and Discovery of Tumor Associated Microbes via RNAseq

**Introduction**

Pandora is a multi-step pipeline to find pathogen sequences in RNAseq data. 
It includes modules for host separation, assembly, blasting contigs, and orf discovery.

**Dependencies**

The following programs must be in your `PATH`:

- python 2.7.x
- java
- [Samtools](http://www.htslib.org/)
- [bamUtil](https://github.com/statgen/bamUtil)
- [STAR](https://github.com/alexdobin/STAR)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
- [BLAST](http://www.ncbi.nlm.nih.gov/books/NBK279671/)

**Additional Files**

Pandora needs various references: a host genome indexed for STAR; a host genome indexed for bowtie2; and the BLAST nucleotide collection (nt) database.
Optionally, you can provide a text file of "blacklist" of non-pathogen taxids for filtering.
You can find the one we use on the Rabadan Lab homepage [here](https://rabadan.c2b2.columbia.edu/public/pandora_resources/), which currently contains any taxid children of the nodes chordata (Taxonomy ID: 7711) or "other sequences" (Taxonomy ID: 28384).

**Workflow**

To accomplish diverse tasks, Pandora has various subcommands (like, say, the program git).
The primary subcommand is `scan`, which is a pipeline comprising the following steps:

1. Subtraction of reads mapping to human genome
2. De-Novo assembly of remaining reads
3. BLAST of assembled contigs
4. ORF search in contigs of unknown origin
5. Remapping of non-human reads to a microbial meta-genome

The `remap` subcommand maps the reads which did not map (from scan, step1) to the assembly of contigs.

**Usage Examples**

```
pandora.py -id 1 -r1 mate_1.fastq.gz -r2 mate_2.fastq.gz --contigthreshold 500 --refstar /path/ref/STAR --refbowtie /path/ref/bowtie -db /path/ref/blastdb/nt --blacklist /path/blacklist.txt --verbose --remap
```

**Notes**

Currently, Pandora makes use of the [Oracle Grid Engine](https://en.wikipedia.org/wiki/Oracle_Grid_Engine).

Status: Active Development
