#!/bin/bash

## Usage:
# sh preprocess_and_harmony.sh

## Create the logs directory
mkdir -p logs

for spefile in spe_postqc.Rdata spe_targeted_postqc.Rdata; do

    ## Internal script name
    SHORT="preprocess_and_harmony_${spefile}"

    # Construct shell file
    echo "Creating script preprocess_and_harmony_${spefile}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -pe local 4
#$ -N ${SHORT}
#$ -o logs/${SHORT}.txt
#$ -e logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript -e "options(width = 120); print('${spefile}'); sessioninfo::session_info()"

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


EOF

    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
