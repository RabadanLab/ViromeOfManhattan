Pandora
=======

Identification and Discovery of Tumor Associated Microbes via RNAseq

**Introduction**

Pandora is a multi-step pipeline to find pathogen sequences in RNAseq data. 
It includes modules for host separation, assembly, blasting contigs, and orf discovery.

**Dependencies**

The following programs must be in your `PATH`:

- python (version >= 2.7)
- java
- [Samtools](http://www.htslib.org/)
- [bamUtil](https://github.com/statgen/bamUtil)
- [STAR](https://github.com/alexdobin/STAR)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
- [BLAST](http://www.ncbi.nlm.nih.gov/books/NBK279671/)
- [Prodigal](http://prodigal.ornl.gov/)

The following Python packages are required:

- [BioPython](http://biopython.org/wiki/Main_Page)

**Additional Files**

Pandora needs various references: a host genome indexed for STAR; a host genome indexed for bowtie2; and the BLAST nucleotide collection (nt) database.

**Workflow**

The Pandora pipeline is organized into the following steps:

1. Removal of host (non-pathogen) reads 
2. Assembly of remaining reads
3. Blasting of the assembled contigs
4. ORF discovery
5. Generate summary report

**Usage Examples**

```
pandora.py -id 1 -r1 mate_1.fastq.gz -r2 mate_2.fastq.gz --contigthreshold 500 --refstar /path/ref/STAR --refbowtie /path/ref/bowtie -db /path/ref/blastdb/nt
```

**Notes**

Currently, Pandora makes use of the [Oracle Grid Engine](https://en.wikipedia.org/wiki/Oracle_Grid_Engine).

Status: Active Development
