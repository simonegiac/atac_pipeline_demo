#!/bin/bash
# alignment.sh
# Created by Simone Giacometti on 2/2/19

# This script is used to map both paired-end (PE) and single-end (SE) ATAC-Seq data with BWA

dirtools=/netapp/simone.giaco/tools
samtools=${dirtools}/samtools-1.9/samtools
bedtools=${dirtools}/bedtoolread2/bin/bedtools

# variables definided with qsub
homedir=$1
sample=$2
genome=$3
end=$4

echo $homedir
echo $sample

# create folders for outputs
mkdir ${homedir}/bams
bams=${homedir}/bams

mkdir ${homedir}/alignment.summary
summary=${homedir}/alignment.summary

fastq_path=${homedir}/fastq

####################### Alignment ######################
read1=${fastq_path}/${sample}_R1_trimmed.fastq
read2=${fastq_path}/${sample}_R1_trimmed.fastq

echo $read1
echo $read2

# alignment
# paired end
if [ $end == pe ] ; then echo paired end ; ${dirtools}/bwa-0.7.12/bwa mem -t 10  /netapp/simone.giaco/ref/genomes/${genome}.fa $read1 $read2 -M | samtools view -bS - > ${bams}/${sample}.bam ; fi
# single end
if [ $end == se ] ; then echo single end ; ${dirtools}/bwa-0.7.12/bwa mem -t 10 /netapp/simone.giaco/ref/genomes/${genome}.fa $read1 -M | samtools view -bS - > ${bams}/${sample}.bam ; fi

# sort bam files
java -jar ${dirtools}/picard-tools-1.126/picard.jar SortSam INPUT=${bams}/${sample}.bam OUTPUT=${bams}/${sample}.sorted.bam SO=coordinate TMP_DIR=tmp/
mv ${bams}/${sample}.sorted.bam ${bams}/${sample}.bam

samtools flagstat ${bams}/${sample}.bam > ${summary}/${sample}.flagstat

# get bam index
java -jar ${dirtools}/picard-tools-1.126/picard.jar BuildBamIndex INPUT=${bams}/${sample}.bam TMP_DIR=tmp/

# unmapped reads
${samtools} view -bhf 4 ${bams}/${sample}.bam > ${bams}/${sample}.unmapped.bam
samtools flagstat ${bams}/${sample}.bam > ${summary}/${sample}.unmapped.flagstat
java -jar ${dirtools}/picard-tools-1.126/picard.jar BuildBamIndex INPUT=${bams}/${sample}.unmapped.bam TMP_DIR=tmp/
