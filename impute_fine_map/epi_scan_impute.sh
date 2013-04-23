#!/bin/bash

#$ -N impute_epi_scan
#$ -cwd
#$ -o /clusterdata/josephp/job_reports/
#$ -e /clusterdata/josephp/job_reports/
#$ -S /bin/bash
#$ -t 1-385

if [ -n "${1}" ]; then
  echo "${1}"
  SGE_TASK_ID=${1}
fi
n=${SGE_TASK_ID}


trunkdir="/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/"
outdir="${trunkdir}epi_scan_imput_outputs/"
blockdir="${trunkdir}snp_region_data/"
pheno="${trunkdir}residuals.RData"
sets="${trunkdir}set3.RData"



R --no-save --args ${n} ${pheno} ${sets} ${blockdir} ${outdir} < epi_scan_impute.R


