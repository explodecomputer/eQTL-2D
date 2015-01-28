#!/bin/bash

#$ -N SNP_extract
#$ -cwd
#$ -l h_vmem=10G
#$ -o /clusterdata/josephp/job_reports/
#$ -e /clusterdata/josephp/job_reports/
#$ -S /bin/bash
#$ -t 1-501

if [ -n "${1}" ]; then
  echo "${1}"
  SGE_TASK_ID=${1}
fi

i=${SGE_TASK_ID}

R --no-save --args // ${i} // < epi_empirical.R



