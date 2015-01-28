#!/bin/bash

#$ -N SNP_extract
#$ -cwd
#$ -l h_vmem=4G
#$ -o /clusterdata/josephp/job_reports/
#$ -e /clusterdata/josephp/job_reports/
#$ -S /bin/bash
#$ -t 1-48

if [ -n "${1}" ]; then
  echo "${1}"
  SGE_TASK_ID=${1}
fi

trunkdir="/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Var_eQTL/data/"
resultsdir="${trunkdir}Output_Z/"
datadir="${trunkdir}data_files/"
SnpInfo="${datadir}bsgs_imputed_R2_80_cleaned_stage2_chr_all_SNP_info.txt"
outdir="${trunkdir}SNP_extract/"
snplist="${trunkdir}phil_snps.txt"


i=${SGE_TASK_ID}
snp=`awk -v i=$i '{if(NR == i) {print $1; exit}}' ${snplist}`
chr=`grep ${snp} ${SnpInfo} | awk '{print $1}'`
out="${outdir}${snp}.txt"

if [ "${chr}" -lt "10" ]; then
	grep ${snp} ${resultsdir}ILMN_*/ILMN_*_Z-fastassoc-chr0${chr}.tbl > ${out}
fi

if [ "${chr}" -gt "9" ]; then
	grep ${snp} ${resultsdir}ILMN_*/ILMN_*_Z-fastassoc-chr${chr}.tbl > ${out}
fi


