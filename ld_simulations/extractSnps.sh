#!/bin/bash

#$ -N extract
#$ -cwd
#$ -S /bin/bash
#$ -t 1-22
#$ -o job_reports/
#$ -e job_reports/
#$ -l h_vmem=50G

if [ -n "${1}" ]; then
  echo "${1}"
  SGE_TASK_ID=${1}
fi

chr=${SGE_TASK_ID}

plinkrt="/ibscratch/wrayvisscher/imputation/arichu/data/imputed/chr${chr}/arichu_1kg_p1v3_${chr}"
snplist="/clusterdata/uqgheman/repo/eQTL-2D/ld_simulations/allsnps.txt"
output="/clusterdata/uqgheman/repo/eQTL-2D/ld_simulations/data/arichu_${chr}"

plink --noweb --bfile ${plinkrt} --extract ${snplist} --make-bed --out ${output}


# merge into one file
