#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 8
#$ -N spaceranger
#$ -o logs/spaceranger.$TASK_ID.txt
#$ -e logs/spaceranger.$TASK_ID.txt
#$ -m e
#$ -t 1-10
#$ -tc 3

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load SpaceRanger
module load spaceranger/1.3.0

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk "NR==${SGE_TASK_ID}" ../01_bamtofastq/samples.txt)
echo "Processing sample ${SAMPLE}"
date

## Get slide and area
SLIDEPART1=$(echo ${SAMPLE} | cut -c1-6)
SLIDEPART2=$(echo ${SAMPLE} | cut -c7-9)
SLIDE="${SLIDEPART1}-${SLIDEPART2}"
CAPTUREAREA=$(echo ${SAMPLE} | cut -c11-12)    
echo "Slide: ${SLIDE}, capture area: ${CAPTUREAREA}"

## Run SpaceRanger
spaceranger count \
    --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=../../raw-data/FASTQ/spaceranger_our_alignments_nocr11/${SAMPLE}/*/ \
    --image=../../processed-data/Images/VistoSeg/Capture_Areas/${SAMPLE}.tif \
    --slide=${SLIDE} \
    --area=${CAPTUREAREA} \
    --loupe-alignment=../../processed-data/Images/VistoSeg/Capture_Areas/loupe_alignment/${SAMPLE}.json \
    --jobmode=local \
    --localcores=8 \
    --localmem=80

## Move output
echo "Moving results to new location"
date
mkdir -p ../../processed-data/spaceranger/
mv ${SAMPLE} ../../processed-data/spaceranger/

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
