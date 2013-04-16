#!/bin/bash

#$ -N nonadd2
#$ -cwd
#$ -S /bin/bash
#$ -t 1-1959
#$ -o job_reports/
#$ -e job_reports/


set -e

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
output="${rootdir}v4_a/filtered_by_nonadditive/filtered_2_"
minclass=5
snplistfile="Full_Merlin_out_NLP_8.txt"

R --no-save --args ${i} ${resdir} ${hsqdir} ${phendat} ${genodat} ${threshold} ${maxrsq} ${minclass} ${output} ${snplistfile} < filter_raw_nonadditive.R

