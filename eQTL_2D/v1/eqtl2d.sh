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
objfile="../data/.RData"
amatrix="../data/"
residuals="../data/"
outfile="result${id}.txt"

epigpu="../exe/epiGPU"
gcta="../exe/gcta64_test_new"


# Correct phenotype for polygenic effects
# This will make a file called gcta.phen.${id}

R --no-save --args ${objfile} ${id} < setup_gcta.R
${gcta} --reml --grm ${amatrix} --pheno gcta.phen.${id} --residuals --out ${residuals}${id}


# Create the .fam file

R --no-save --args ${famfile} ${residuals}${id} ${id} < setup_epigpu.R

# Run the analysis

${epigpu} -A ${bedfile} ${bimfile} ${famfile}${id} ${outfile} -t i -F 8 -I 12

