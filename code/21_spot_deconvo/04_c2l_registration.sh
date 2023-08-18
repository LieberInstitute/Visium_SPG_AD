#!/bin/bash

#$ -cwd
#$ -N "c2l_registration"
#$ -o ../../processed-data/21_spot_deconvo/logs/04_c2l_registration.log
#$ -e ../../processed-data/21_spot_deconvo/logs/04_c2l_registration.log
#$ -l caracol,mf=64G,h_vmem=64G,h_fsize=50G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

###############################################################################
#   Dynamically select a GPU based on availability
###############################################################################

USAGE_CUTOFF=10
NUM_GPUS=1

avail_gpus=$(
    nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader |
    cut -d " " -f 1 | awk -v usage="$USAGE_CUTOFF" '$1 < usage {print NR - 1}'
)

#  Simply exit with an error if there are no GPUs left
if [[ -z $avail_gpus ]]; then
    echo "No GPUs are available."
    exit 1
fi

export CUDA_VISIBLE_DEVICES=$(
    echo "$avail_gpus" | head -n $NUM_GPUS | paste -sd ","
)

echo "Chose GPU(s): $CUDA_VISIBLE_DEVICES"

###############################################################################
#   Submit the python script
###############################################################################

module load cell2location/0.1.3
python 04_c2l_registration.py

echo "**** Job ends ****"
date