#!/bin/bash

## Usage:
# sh 01_preprocess_and_harmony.sh

## Create the logs directory
mkdir -p logs

for spetype in "wholegenome" "targeted"; do

    ## Internal script name
    SHORT="preprocess_and_harmony_${spetype}"

    # Construct shell file
    echo "Creating script preprocess_and_harmony_${spetype}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=6G,h_vmem=6G,h_fsize=100G
#$ -pe local 4
#$ -N ${SHORT}
#$ -o logs/${SHORT}.txt
#$ -e logs/${SHORT}.txt
#$ -m e
#$ -hold_jid qc_metrics_and_segmentation

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.2

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_preprocess_and_harmony.R -s ${spetype}

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


EOF

    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
