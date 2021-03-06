#!/bin/bash

#$ -N nonmar
#$ -cwd
#$ -S /bin/bash
#$ -t 1-5280
#$ -o job_reports/
#$ -e job_reports/
#$ -l h_vmem=4G


set -e

if [ -n "${1}" ]; then
  echo "${1}"
  SGE_TASK_ID=${1}
fi

i=${SGE_TASK_ID}


rootdir="/hpscratch/wrayvisscher/gib/git/eQTL-2D/"
resdir="${rootdir}run/results/result"
hsqdir="${rootdir}run/scratch/resphen"
phendat="${rootdir}data/residuals.RData"
genodat="${rootdir}data/geno.RData"
threshold=13
maxrsq=0.1
output="${rootdir}filtering/filtered_by_nonmarginal/filtered"
minclass=5
snplistfile="${rootdir}filtering/marginal_lists/marginal_list.RData"

R --no-save --args ${i} ${resdir} ${hsqdir} ${phendat} ${genodat} ${threshold} ${maxrsq} ${minclass} ${output} ${snplistfile} < filter_raw_nonmarginal.R

