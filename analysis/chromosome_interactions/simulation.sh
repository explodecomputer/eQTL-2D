#!/bin/bash

#$ -N ci
#$ -cwd
#$ -l h_vmem=1G
#$ -t 1-1000
#$ -S /bin/bash
#$ -o job_reports/
#$ -e job_reports/

set -e

if [ -n "${1}" ]; then
  echo "${1}"
  SGE_TASK_ID=${1}
fi

i=${SGE_TASK_ID}

cifile="~/repo/eQTL-2D/data/supFile3_K562_interactingLoci_clusters.csv"
bimfile="~/repo/eQTL-2D/data/clean_geno_final.bim"
n=549
win=25000
outroot="~/repo/eQTL-2D/analysis/chromosome_interactions/out"

R --no-save --args ${i} ${n} ${win} ${cifile} ${bimfile} ${outroot} < simulation.R
