#!/bin/bash

## Usage:
# sh 02_explore_expr_variability.sh

## Create the logs directory
mkdir -p logs

for spetype in "wholegenome" "targeted"; do #had to include quotes here manually

    ## Internal script name
    SHORT="explore_expr_variability_${spetype}"

    # Construct shell file
    echo "Creating script 02_explore_expr_variability_${spetype}"
    cat > ${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=2G,h_vmem=2G,h_fsize=100G
#$ -N ${SHORT}
#$ -o logs/${SHORT}.txt
#$ -e logs/${SHORT}.txt
#$ -m e
#$ -hold_jid create_pseudobulk_data_${spetype}

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.2.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 02_explore_expr_variability.R -s ${spetype}

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


EOF

    call="qsub ${SHORT}.sh"
    echo $call
    $call
done
