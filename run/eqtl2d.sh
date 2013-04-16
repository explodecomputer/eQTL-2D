#!/bin/bash

#PBS -W group_list=director553
#PBS -q workq
#PBS -l select=1:ncpus=12:ngpus=1:mem=4000mb
#PBS -l place=excl
#PBS -N epistasis
#PBS -J 184-185
#PBS -l walltime=12:00:00
#PBS -o job_reports/
#PBS -e job_reports/

set -e

rootdir="/home/ghemani/repo/eqtl2d/"
epigpu="/home/ghemani/repo/epiGPU/exe/linux/epiGPU"

cd ${rootdir}run/

module add R
module add cuda
# module add cuda-sdk

if [ -n "${1}" ]; then
  PBS_ARRAY_INDEX=${1}
fi
id=${PBS_ARRAY_INDEX}

if [ ! -d "scratch/" ]; then
  mkdir scratch/
fi

if [ ! -d "results/" ]; then
  mkdir results/
fi


# These are the files that will be used in the analysis
bplink="${rootdir}data/clean_geno_final"
egufile="${bplink}.egu"
objfile="${rootdir}data/eqtl2d_objects.RData"
prbfile="${rootdir}data/probes5381-7339.RData"

# Outputs
resphen="scratch/resphen_2_${id}"
outfile="results/result_2_${id}.txt"

R="R"

# Correct phenotype for polygenic effects

${R} --no-save --args ${objfile} ${prbfile} ${id} ${resphen} < polygenic.R

# Run the analysis

${epigpu} -A ${egufile} ${outfile} -i 2048 -t i -F 9.5 -I 16 -P ${resphen}

gzip ${outfile}

