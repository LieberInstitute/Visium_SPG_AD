#!/bin/bash
#$ -pe local 8
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -o /dcl01/lieber/ajaffe/Maddy/RNAscope/VisiumIF/logs/spaceranger_VIFAD3A1.txt
#$ -e /dcl01/lieber/ajaffe/Maddy/RNAscope/VisiumIF/logs/spaceranger_VIFAD3A1.txt
#$ -V
#$ -cwd

# run on JHPCE cluster
# qsub scripts/spaceranger/spaceranger_round3_ant.sh

############################
# Script to run Space Ranger
############################

# locations of files:
# -------------------
# spaceranger reference: /dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A
# fastq (NextSeq): /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-09_ASpa061421_01/ and /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-07-06_ASpa061421/
# images (raw): /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3

# summary spreadsheet: /dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC
# - contains sample ID, sample name, slide serial number, capture area ID


# load spaceranger module
module use /jhpce/shared/jhpce/modulefiles/libd
module load spaceranger

# run in outputs directory (spaceranger can only save outputs in current working directory)
# run spaceranger count for each sample

spaceranger count \
--id=VIFAD3_V10T31-036_A1 \
--transcriptome=/dcl02/lieber/ajaffe/SpatialTranscriptomics/refdata-gex-GRCh38-2020-A \
--fastqs=/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/raw-data/FASTQ/spaceranger_our_alignments_nocr11/V10T31036_A1_Br3874/1145930_1160915_0_1_H5GKNDSX2/ \
--darkimage=/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_A1.tif \
--slide=V10T31-036 \
--area=A1 \
--loupe-alignment=/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/loupe_alignment/VIFAD3_V10T31-036_A1.json \
--jobmode=local \
--localcores=8 \
--localmem=64


# restore working directory
cd $cwd