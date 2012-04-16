#!/bin/bash


#PBS -N eqtl2d
#PBS -J 1-7000
#PBS -q gpu_queue
#PBS -l walltime=12:00:00,nodes=1

# read in pedigree and phenotype information

id=${PBS_ARRAY_INDEX}

bimfile="../data/.bim"
bedfile="../data/.bed"
famfile="../data/.fam"
phefile="../data/.pheno"

epigpu="../exe/epiGPU"



# Correct phenotype for polygenic effects
# This will make a file called ${famfile}${id}

R --no-save --args ${phefile} ${famfile} ${id} < generate_famfile.R



# Run the analysis

${epigpu} -A ${bedfile} ${bimfile} ${famfile}${id} ${outfile} -t i -F 8 -I 12




