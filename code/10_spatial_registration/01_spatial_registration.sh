#!/bin/bash
# indicates the file is a bash script file

## Usage:
# sh 01_spatial_registration.sh
#an example of how to run this script

mkdir -p logs

for spetype in wholegenome targeted; do

    SHORT="spatial_registration_${spetype}"
    #create a new var

    ##construct shell file
    echo "Creating script spatial_registration_${spetype}"

    #write to a hidden file using "End of file"
    cat > ${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N ${SHORT}
#$ -o logs/${SHORT}.\$TASK_ID.txt
#$ -e logs/${SHORT}.\$TASK_ID.txt
#$ -m e
#$ -t 2-28
#$ -tc 20


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/devel

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_spatial_registration.R -s ${spetype}

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/

EOF

    call="qsub ${SHORT}.sh"
    echo $call
    $call
done
