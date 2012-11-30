#!/bin/bash

#$ -N nonadd
#$ -cwd
#$ -S /bin/bash
#$ -t 1-5280
#$ -o job_reports/
#$ -e job_reports/

# SGE_TASK_ID=${1}
i=${SGE_TASK_ID}

rootdir="/hpscratch/wrayvisscher/gib/git/wrayvisscher/eQTL_2D/"
resdir="${rootdir}v4/results/result"
hsqdir="${rootdir}v4/scratch/resphen"
phendat="${rootdir}data/residuals.RData"
genodat="${rootdir}data/geno.RData"
threshold=15.38
maxrsq=0.1
output="${rootdir}v4_a/filtered_by_nonadditive/filtered"
minclass=5
snplistfile="Full_Merlin_out_NLP_8.txt"

R --no-save --args ${i} ${resdir} ${hsqdir} ${phendat} ${genodat} ${threshold} ${maxrsq} ${minclass} ${output} ${snplistfile} < filter_raw_nonadditive.R

