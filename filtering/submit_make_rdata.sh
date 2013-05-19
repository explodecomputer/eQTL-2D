#!/bin/bash

#$ -N makerdata
#$ -cwd
#$ -S /bin/bash
#$ -t 1-7339
#$ -o job_reports/
#$ -e job_reports/
#$ -l h_vmem=10G


set -e

if [ -n "${1}" ]; then
  echo "${1}"
  SGE_TASK_ID=${1}
fi

i=${SGE_TASK_ID}


rootdir="/hpscratch/wrayvisscher/gib/git/eQTL-2D/"
resdir="${rootdir}run/results/result"
hsqdir="${rootdir}run/scratch/resphen"
phendat="${rootdir}data/residuals_all.RData"
output="${resdir}"

R --no-save --args ${i} ${resdir} ${hsqdir} ${phendat} ${output} < make_rdata.R

