#!/bin/bash

#PBS -N ci2
#PBS -o /home/tghemani/repo/ibd_ibs/grm/scripts/job_reports/
#PBS -e /home/tghemani/repo/ibd_ibs/grm/scripts/job_reports/
#PBS -J 1-10000


set -e

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAY_INDEX=${1}
fi

i=${PBS_ARRAY_INDEX}

cifile="~/repo/eQTL-2D/data/supFile3_K562_interactingLoci_clusters.csv"
bimfile="~/repo/eQTL-2D/data/hg19locations.bim"
n=549
win=10000
outroot="~/repo/eQTL-2D/analysis/chromosome_interactions/results/out"

R --no-save --args ${i} ${n} ${win} ${cifile} ${bimfile} ${outroot} < ~/repo/eQTL-2D/analysis/chromosome_interactions/simulation.R
