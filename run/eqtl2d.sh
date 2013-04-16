#!/bin/bash

#PBS -A Desc006
#PBS -l nodes=1:gpus=1
#PBS -l pmem=4000MB
#PBS -N epistasis
#PBS -t 1501-1959
#PBS -l walltime=12:00:00
#PBS -o job_reports/
#PBS -e job_reports/

set -e

rootdir="/home/projects/Desc006/repo/eqtl2d/"
epigpu="/home/projects/Desc006/repo/epiGPU/exe/epiGPU"

cd ${rootdir}run/

module add R
module add cuda
# module add cuda-sdk

## PBS_ARRAY_INDEX=${1}

if [ -n "${1}" ]; then
  PBS_ARRAYID=${1}
fi
id=${PBS_ARRAYID}

if [ ! -d "scratch/" ]; then
  mkdir scratch/
fi

if [ ! -d "results/" ]; then
  mkdir results/
fi

# Not really required anymore
bplink="${rootdir}data/clean_geno_final"
bimfile="${bplink}.bim"
bedfile="${bplink}.bed"
famfile="${bplink}.fam"
phenfile="../data/probe_pheno_plink_final.txt"
covfile="../data/covariates_plink_final.txt"

# These are the files that will be used in the analysis
egufile="${bplink}.egu"
objfile="../data/eqtl2d_objects.RData"
prbfile="../data/probes5381-7339.RData"

# Outputs
resphen="scratch/resphen_2_${id}"
outfile="results/result_2_${id}.txt"

R="R"

# Correct phenotype for polygenic effects

${R} --no-save --args ${objfile} ${prbfile} ${id} ${resphen} < polygenic.R

# Run the analysis

${epigpu} -A ${egufile} ${outfile} -i 2048 -t i -F 9.5 -I 16 -P ${resphen}

gzip ${outfile}

