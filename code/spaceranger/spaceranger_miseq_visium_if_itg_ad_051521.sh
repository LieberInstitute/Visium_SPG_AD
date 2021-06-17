#!/bin/bash
#$ -pe local 8
#$ -l mem_free=5G,h_vmem=10G,h_fsize=100G
#$ -V
#$ -cwd

# run on JHPCE cluster
# qsub scripts/spaceranger/spaceranger_NextSeq_dlpfc_post.sh

############################
# Script to run Space Ranger
############################

# locations of files:
# -------------------
# spaceranger reference: /dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A
# fastq (NextSeq): /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/raw-data/FASTQ/FASTQ_10x_2021-05-03/
# images (raw): /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/raw-data/Images/

# summary spreadsheet: /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/sample_info/Visium dlpfc 100520 Master.xlsx
# - contains sample ID, sample name, slide serial number, capture area ID


# load spaceranger module
module use /jhpce/shared/jhpce/modulefiles/libd
module load spaceranger/1.2.2

# run in outputs directory (spaceranger can only save outputs in current working directory)
cwd="/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD"
mkdir -p processed-data/NovaSeq
cd processed-data/NovaSeq


# run spaceranger count for each sample

spaceranger count \
--id=1_Br3874_ITG \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/raw-data/FASTQ/FASTQ_10x_2021-05-03 \
--sample=1 \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/raw-data/Images/1_V10A27-004_A1_Br3874.tiff \
--slide=V10A27-004 \
--area=A1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/raw-data/Images/json/fluorescent/V10A27-004-A1.json \
--jobmode=local \
--localcores=8 \
--localmem=64



# restore working directory
cd $cwd