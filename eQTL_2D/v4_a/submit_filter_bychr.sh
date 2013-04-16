#!/bin/bash

#$ -N epistasis
#$ -cwd
#$ -S /bin/bash
#$ -t 1-1959
#$ -o job_reports/
#$ -e job_reports/

if [ -n "${1}" ]; then
  echo "${1}"
  SGE_TASK_ID=${1}
fi

i=${SGE_TASK_ID}

rootdir="/hpscratch/wrayvisscher/gib/git/wrayvisscher/eQTL_2D/"
resdir="${rootdir}v4/results/result_2_"
hsqdir="${rootdir}v4/scratch/resphen_2_"
phendat="${rootdir}data/residuals2.RData"
genodat="${rootdir}data/geno.RData"
threshold=15.38
maxrsq=0.1
output="${rootdir}v4_a/filtered_by_chr/filtered_2_"
minclass=5

R --no-save --args ${i} ${resdir} ${hsqdir} ${phendat} ${genodat} ${threshold} ${maxrsq} ${minclass} ${output} < filter_raw_bychr.R

