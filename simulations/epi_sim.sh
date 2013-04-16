#!/bin/bash

#$ -N epi_sim_cis_cis
#$ -cwd
#$ -o /clusterdata/josephp/job_reports/
#$ -e /clusterdata/josephp/job_reports/
#$ -S /bin/bash
#$ -t 1-1000

if [ -n "${1}" ]; then
  echo "${1}"
  SGE_TASK_ID=${1}
fi
n=${SGE_TASK_ID}


trunkdir="/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/"
outdir="${trunkdir}sim_output_files/"
blockdir="${trunkdir}sim_impute_regions/"
pheno="${trunkdir}residuals.RData"
geno="${trunkdir}geno.RData"
g_info="${trunkdir}GWAS.map"
i_info="${trunkdir}bsgs_imputed_R2_80_cleaned_stage2_chr_all_SNP_info.txt"
type="cis_cis" 

R --no-save --args ${n} ${pheno} ${geno} ${outdir} ${blockdir} ${g_info} ${i_info} ${type} < epi_scan_impute.R


