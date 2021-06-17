#!/bin/bash
#$ -pe local 8
#$ -l mem_free=10G,h_vmem=20G,h_fsize=100G
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
# fastq (NextSeq): /dcl02/lieber/ajaffe/Nina/Kristen/022621_Visium_MiSeq/spaceranger_out/tiny-bcl/outs/fastq_path/CYT55
# images (raw): /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/Images/

# summary spreadsheet: /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/sample_info/Visium dlpfc 100520 Master.xlsx
# - contains sample ID, sample name, slide serial number, capture area ID


# load spaceranger module
module use /jhpce/shared/jhpce/modulefiles/libd
module load spaceranger

# run in outputs directory (spaceranger can only save outputs in current working directory)
cwd=$(pwd)
mkdir -p outputs/MiSeq
cd outputs/MiSeq


# run spaceranger count for each sample

spaceranger count \
--id=Br3880_ITG_montage \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/FASTQ/Br3880_ITG \
--image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/Images/json/montage/Br3880_Montage.tif \
--slide=V10A27-004 \
--area=D1 \
--loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/Images/json/montage/V10A27-004-D1_montage.json \
--jobmode=local \
--localcores=8 \
--localmem=64

# spaceranger count \
# --id=Br3873_ITG_montage \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/FASTQ/Br3873_ITG \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/Images/json/montage/Br3873_Montage.tif \
# --slide=V10A27-004 \
# --area=C1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/Images/json/montage/V10A27-004-C1_montage.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

# spaceranger count \
# --id=Br3854_ITG_montage \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/FASTQ/Br3854_ITG \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/Images/json/montage/Br3854_Montage.tif \
# --slide=V10A27-004 \
# --area=B1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/Images/json/montage/V10A27-004-B1_montage.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

# spaceranger count \
# --id=Br3874_ITG_montage \
# --transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
# --fastqs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/FASTQ/Br3874_ITG \
# --image=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/Images/json/montage/Br3874_Montage.tif \
# --slide=V10A27-004 \
# --area=A1 \
# --loupe-alignment=/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/Visium_IF_AD/Images/json/montage/V10A27-004-A1_montage.json \
# --jobmode=local \
# --localcores=8 \
# --localmem=64

# restore working directory
cd $cwd