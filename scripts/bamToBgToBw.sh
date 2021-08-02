#!/bin/bash

dirtools=/netapp/simone.giaco/tools
bedtools=${dirtools}/bedtools2/bin/bedtools
bamToBgToBw=${dirtools}/bamToBgToBw/bin/

# variables definided with qsub
homedir=$1
sample=$2
genome=$3
bams=${homedir}/bams

echo $sample
echo $genome

if [ $genome == hg19_ucsc ] ; then g=hg19 ; fi
if [ $genome == mm9 ] ; then g=mm9 ; fi
echo $g

# bam to bedgrah
cd ${bams}
bedtools genomecov -ibam ${sample}.un.noblr.nodups.bam -bg -g ${g} | sort -k1,1 -k2,2n > ${sample}.un.noblr.nodups.bg

# bedgraph to bigwig
${bamToBgToBw}/bedGraphToBigWig ${sample}.un.noblr.nodups.bg ${bamToBgToBw}/${genome}.chrom.sizes ${sample}.un.noblr.nodups.bw
