#!/bin/bash

# variables definided with qsub
homedir=$1
sample=$2
genome=$3
end=$4

echo $sample
echo $genome
echo $end

mkdir $homedir/macs2
macs2=${homedir}/macs2

if [ $end == pe ] ; then end_seq=BAMPE ; fi
if [ $end == se ] ; then end_seq=BAM ; fi
echo $end_seq


###################### Call Peaks ######################

# peak calling with Macs2
export PYTHONPATH=/netapp/simone.giaco/tools/MACS2-2.1.0.20140616/lib64/python2.7/site-packages

scl enable python27 "/netapp/simone.giaco/tools/MACS2-2.1.0.20140616/bin/macs2 callpeak \
  -t ${homedir}/bams/${sample}.un.noblr.nodups.bam -f ${end_seq} -g ${genome} --outdir ${macs2} \
  --nomodel --extsize 200 --shift -100 --bdg --q 0.05 --name ${sample}.un.noblr.nodups.peak"
