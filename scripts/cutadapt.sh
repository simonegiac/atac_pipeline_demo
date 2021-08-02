#!/bin/bash
# This script is used to check the quality of the reads and remove sequencer adaptors, if present.

# Variables definided with qsub
homedir=$1
sample=$2

fastq_path=${homedir}/fastq

read1=${fastq_path}/${sample}_R1.fastq
read2=${fastq_path}/${sample}_R2.fastq
echo $read1
echo $read2

# Parameters for removal of Illumina adaptors from both reads. "m"=mimimum read length. "q"=minimum quality  "a","A"=Illumina Universal adaptors
CUTADAPT_ARGS="-m 10 -q 20,20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

# Check quality with Fastqc
mkdir ${homedir}/fastqc
fastqc $read1 -o ${fastqc}
fastqc $read2 -o ${fastqc}

# Remove adaptors with Cutadapt
cutadapt ${CUTADAPT_ARGS} -o ${fastq_path}/${sample}_R1_trimmed.fastq -p ${fastq_path}/${sample}_R2_trimmed.fastq ${read1} ${read2}

# Run Fastqc after trimming adapter sequences
fastqc ${fastq_path}/${sample}_R1_trimmed.fastq -o ${homedir}/fastqc
fastqc ${fastq_path}/${sample}_R2_trimmed.fastq -o ${homedir}/fastqc
