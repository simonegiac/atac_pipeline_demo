#!/bin/bash
# remove.blr.sh
# Created by Simone Giacometti on 2/2/19

# This script is used to remove remove mitochondrial reads, black listed regions and duplicate reads 


dirtools=/netapp/simone.giaco/tools
samtools=${dirtools}/samtools-1.9/samtools
bedtools=${dirtools}/bedtools2/bin/bedtools

# variables definided with qsub
homedir=$1
sample=$2
genome=$3

echo $homedir
echo $sample
echo $genome

mkdir ${homedir}/bedtools
beds=${homedir}/bedtools

bams=${homedir}/bams
summary=${homedir}/alignment.summary

# get header for mitochondrial DNA
samtools view -H ${bams}/${sample}.bam > ${bams}/${sample}.header

# get mitochondrial and nuclear reads (in human/mouse genome builds, the mitochondrial genome is labeled ‘chrM’)
$samtools view -h ${bams}/${sample}.bam | awk '{if($3 == "chrM"){print $0}}' | cat ${bams}/${sample}.header - |samtools view -Sb - > ${bams}/${sample}.mito.bam
$samtools view -h ${bams}/${sample}.bam | awk '{if($3 != "chrM"){print $0}}' | samtools view -Sb - > ${bams}/${sample}.nuc.bam
rm ${bams}/${sample}.header
samtools flagstat ${bams}/${sample}.mito.bam > ${summary}/${sample}.mito.flagstat
samtools flagstat ${bams}/${sample}.nuc.bam > ${summary}/${sample}.nuc.flagstat

# get nuclear bam index
java -jar ${dirtools}/picard-tools-1.126/picard.jar BuildBamIndex INPUT=${bams}/${sample}.nuc.bam TMP_DIR=tmp/

# remove black listed regions and multiply mapped reads
$samtools view -bF 4 ${bams}/${sample}.nuc.bam | ${bedtools} intersect -v -abam - -b /netapp/simone.giaco/ref/genomes/${genome}-blacklisted.bed > ${bams}/${sample}.un.noblr.bam
samtools flagstat ${bams}/${sample}.un.noblr.bam > ${summary}/${sample}.un.noblr.flagstat

# get black listed regions and multiple mapped reads bam index
java -jar ${dirtools}/picard-tools-1.126/picard.jar BuildBamIndex INPUT=${bams}/${sample}.un.noblr.bam TMP_DIR=tmp/

# get unique reads
java -jar ${dirtools}/picard-tools-1.126/picard.jar MarkDuplicates INPUT=${bams}/${sample}.un.noblr.bam OUTPUT=${bams}/${sample}.un.noblr.nodups.bam REMOVE_DUPLICATES=TRUE METRICS_FILE=${bams}/${sample}.markduplicates.txt TMP_DIR=tmp/

# sort bam file
$samtools sort ${bams}/${sample}.un.noblr.bam -o ${bams}/${sample}.un.noblr.sorted.bam
mv ${bams}/${sample}.un.noblr.sorted.bam ${bams}/${sample}.un.noblr.bam

# get uniquely mapped noblr bam index
java -jar ${dirtools}/picard-tools-1.126/picard.jar BuildBamIndex INPUT=${bams}/${sample}.un.noblr.nodups.bam TMP_DIR=tmp/
samtools flagstat ${bams}/${sample}.un.noblr.nodups.bam > ${summary}/${sample}.un.noblr.nodups.flagstat
