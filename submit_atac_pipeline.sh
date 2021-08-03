#!/bin/bash

# submit_atac_pipeline.sh
# Created by Simone Giacometti on 2/2/19

# Instructions to SGE:
#$ -S /bin/bash                     #-- the shell for the job
##$ -o [dir]                        #-- output directory (fill in)
##$ -e [dir]                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=20G                  #-- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               #-- SGE resources (CPU type)
#$ -l netapp=1G,scratch=1G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)
##$ -t 1-10                        #-- remove first '#' to specify the number of tasks

##### Example input: qsub -v names=names -v genome=hg19 -v end=pe -o logs -e logs -t 1-10 sumit_atac_pipeline.sh #####
# -v var[=value] passes environment variable "var" to the job:
  # "names" refer to a "names" configuration file containing the names of the fastq samples without the *_R*.fastq suffix
  # "genome" can be "hg19" (human) or "mm9" (mouse)
  # "end" can be "pe" (paired ended) or "se" (single ended)
# With the -t option, this script_path will be run for each task, but with $SGE_TASK_ID set to a different value each time.

# Requirements:
# samtools ==1.9
# sambamba ==0.6.6
# samstats ==0.2.1
# bedtools ==2.29.0
# deeptools ==3.3.1
# picard ==2.20.7
# pysam ==0.15.3

# ucsc-wigtobigwig ==357
# ucsc-bedgraphtobigwig ==357
# ucsc-bigwiginfo ==357

# cutadapt ==3.3.4
# bwa ==0.7.17
# macs2 ==2.2.4
# r ==3.5.1
# java-jdk

# wget
# pyopenssl
# grep
# tar
# sed
# coreutils  # GNU split and shuf

homedir=$PWD
script_path=${homedir}/scripts

# Location of Fastq files. To be transferred to the server before running the pipeline
fastq_path=${homedir}/fastq

mkdir ${homedir}/logs
export TMPDIR=${homedir}/logs

sample=$(awk 'FNR == i {print}' i=${SGE_TASK_ID} ${homedir}/${names})

echo $names
echo $sample
echo $homedir
echo $genome
echo $end

# Select genome (human or mouse)
if [ $genome == hg19_ucsc ] ; then macs_genome=hs ; fi
if [ $genome == mm9 ] ; then macs_genome=mm ; fi
echo $macs_genome

##### Check read quality and remove Illumina adapters #####
${script_path}/cutadapt.sh ${homedir} ${sample}

##### Align reads with BWA #####
${script_path}/alignment.sh ${homedir} ${sample} ${genome} ${end}

##### Remove mitochondrial reads, black listed regions and duplicate reads #####
${script_path}/remove.blr.sh ${homedir} ${sample} ${genome}

##### Call peaks with Macs2 for ATAC options #####
${script_path}/macs2.sh ${homedir} ${sample} ${macs_genome} ${end}

##### Make bigwig files for visualization #####
${script_path}/bamToBgToBw.sh ${sample} ${genome}

##### Annotate peaks with CHIPseeker and motif discovery with Homer #####
# Not implemented in the current pipeline. To be performed manually depending on the experiment

# End-of-job summary
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
