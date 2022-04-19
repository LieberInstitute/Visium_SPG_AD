#!/bin/bash

## Usage:
# sh 03_scatter_plots_case_vs_controls.sh

## Create the logs directory
mkdir -p logs

for spetype in wholegenome targeted; do

    ## Internal script name
    SHORT="03_scatter_plots_case_vs_controls_${spetype}"

    # Construct shell file
    echo "Creating script 03_scatter_plots_case_vs_controls_${spetype}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
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
