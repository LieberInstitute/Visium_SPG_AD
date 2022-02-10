#!/bin/bash

## Usage:
# sh BayesSpace_k_search.sh

## Create the logs directory
mkdir -p logs

for spefile in spe_harmony_wholegenome spe_harmony_targeted; do

    ## Internal script name
    SHORT="BayesSpace_k_search_${spefile}"

    ## Map to older names in order to use a dynamic hold_jid below
    if [ $spefile == "spe_harmony_wholegenome" ] ; then
    	speprevious="spe_postqc"
    else
    	speprevious="spe_targeted_postqc"
    fi

    # Construct shell file
    echo "Creating script BayesSpace_k_search_${spefile}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N ${SHORT}
#$ -o logs/${SHORT}.\$TASK_ID.txt
#$ -e logs/${SHORT}.\$TASK_ID.txt
#$ -m e
#$ -t 4-15
#$ -tc 20
#$ -hold_jid preprocess_and_harmony_${speprevious}

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.1.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 02_BayesSpace_k_search.R -s ${spefile}

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


EOF

    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
