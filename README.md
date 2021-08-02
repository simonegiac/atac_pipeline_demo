## Overview

This is a pipeline demo specifically put together for AIOnco.

## Assay overview

ATAC-Seq (Assay for Transposase-Accessible Chromatin with high-throughput sequencing) is a method for determining chromatin accessibility across the genome. It utilizes a hyperactive Tn5 transposase to insert sequencing adapters into open chromatin regions. High-throughput sequencing then yields reads that indicate these regions of increased accessibility.

## Pipeline overview

This pipeline is designed for automated and parallel end-to-end quality control and processing of ATAC-seq data.  It can be run end-to-end, starting from raw FASTQ files all the way to peak calling and signal track generation using a single submit command. Illumina paired-end data from human and mouse samples are supported.

This pipeline is designed to run on compute clusters with a job submission engine (e.g SGE), and was originally created for Wynton (a UCSF HPC). Programs are run on Wynton by submitting jobs to SGE via a qsub command.

**Usage**: `qsub -v names=names -v genome=hg19 -v end=pe -o logs -e logs -t 1-10 sumit_atac_pipeline.sh`

`-v var[=value]` passes environment variable `var` to the job:


 - `names` refer to a `names` configuration file containing the names of
   the samples without the `*_R*.fastq` suffix
 - genome can be `hg19` (human) or `mm9` (mouse)
- `end` can be `pe` (paired ended) or `se` (single ended)
- with the `-t` option, this script_path will be run for each task, but with `$SGE_TASK_ID` set to a different value each time.

**Inputs**: Fastq Ilumina reads (paired-ended, unzipped)

**Main Outputs**: Alignments and filtered alignments in BAM format, alignment summary and QC statistics, peaks in BED format, and bigwig signals.

**Requirements**:

 - Paired-end Illumina fastq files, unzipped
 - Genome reference sequences
   in fasta format (and corresponding indexes) for mm9 and hg19
 - `name` configuration file with the names of the samples without the `*_R*.fastq` suffix (one name per line)
 - samtools ==1.9
 - sambamba ==0.6.6
 - samstats ==0.2.1
 - bedtools ==2.29.0

 - deeptools ==3.3.1     picard ==2.20.7     pysam ==0.15.3
 - ucsc-wigtobigwig ==357
 - ucsc-bedgraphtobigwig ==357
 - ucsc-bigwiginfo==357


 - cutadapt ==3.3.4     
 - bwa ==0.7.17     
 - macs2 ==2.2.4     
 - r ==3.5.1
 - java-jdk



 - wget
 - pyopenssl
 - grep
 - tar
 - sed
 - coreutils  # GNU split and shuf
