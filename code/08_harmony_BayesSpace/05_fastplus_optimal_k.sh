#!/bin/bash

## Usage:
# sh 05_fastplus_optimal_k.sh

## Create the logs directory
mkdir -p logs

for spetype in wholegenome targeted; do

    ## Internal script name
    SHORT="05_fastplus_optimal_k_${spetype}"

    # Construct shell file
    echo "Creating script 05_fastplus_optimal_k_${spetype}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=80G,h_vmem=80G,h_fsize=100G
#$ -N ${SHORT}
#$ -o logs/${SHORT}.\$TASK_ID.txt
#$ -e logs/${SHORT}.\$TASK_ID.txt
#$ -m e
#$ -t 2-15
#$ -tc 14

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
Rscript -e "options(width = 120); print('${spetype}'); sessioninfo::session_info()"

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


EOF

    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
