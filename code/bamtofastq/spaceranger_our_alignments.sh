#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=1G,h_vmem=1G,h_fsize=200G
#$ -pe local 4
#$ -N spaceranger_our_alignments
#$ -o logs/spaceranger_our_alignments.$TASK_ID.txt
#$ -e logs/spaceranger_our_alignments.$TASK_ID.txt
#$ -m e
#$ -t 1-10
#$ -tc 5

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load bamtofastq
module load bamtofastq/1.3.2

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk "NR==${SGE_TASK_ID}" samples.txt)
echo "Processing sample ${SAMPLE}"
date

## List current input file
INPUTBAM="../../raw-data/10x_files/Lieber_Transfer/${SAMPLE}/possorted_genome_bam.bam"
ls -lh ${INPUTBAM}

## Create output directory
mkdir -p ../../raw-data/FASTQ/spaceranger_our_alignments
OUTPUTDIR="../../raw-data/FASTQ/spaceranger_our_alignments/${SAMPLE}/"

## Run bamtofastq
bamtofastq --nthreads=4 --traceback --cr11 ${INPUTBAM} ${OUTPUTDIR}

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
