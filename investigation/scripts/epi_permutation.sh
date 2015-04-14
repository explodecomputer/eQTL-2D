#!/bin/bash

#$ -N epi_permutation
#$ -cwd
#$ -l h_vmem=500M
#$ -o /ibscratch/wrayvisscher/josephP/job_output/
#$ -e /ibscratch/wrayvisscher/josephP/job_output/
#$ -S /bin/bash
#$ -t 1-501

if [ -n "${1}" ]; then
  echo "${1}"
  SGE_TASK_ID=${1}
fi

i=${SGE_TASK_ID}

R --no-save --args /ibscratch/wrayvisscher/josephP/epi_investigation/investigation_data.RData ${i} /ibscratch/wrayvisscher/josephP/epi_investigation/output_permutation/ < epi_permutation.R



