#!/bin/bash

#$ -N nonmarsig
#$ -cwd
#$ -S /bin/bash
#$ -t 1-7339
#$ -o job_reports/
#$ -e job_reports/
#$ -l h_vmem=5G


set -e

if [ -n "${1}" ]; then
  echo "${1}"
  SGE_TASK_ID=${1}
fi

i=${SGE_TASK_ID}


rootdir="/hpscratch/wrayvisscher/gib/git/eQTL-2D/"
resdir="${rootdir}run/results/result"
phendat="${rootdir}data/residuals_all.RData"
genodat="${rootdir}data/geno.RData"
threshold=13
threshold2=10.89
maxrsq=0.1
output="${rootdir}filtering/filtered_by_nonmarginal_sig/filtered"
minclass=5
snplistfile="${rootdir}filtering/marginal_lists/marginal_list.RData"

R --no-save --args ${i} ${resdir} ${phendat} ${genodat} ${threshold} ${maxrsq} ${minclass} ${output} ${snplistfile} ${threshold2} < filter_rdata_nonmarginal_sig.R

