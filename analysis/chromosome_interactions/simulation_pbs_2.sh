#!/bin/bash

#PBS -N ci2
#PBS -o /home/tghemani/repo/ibd_ibs/grm/scripts/job_reports/
#PBS -e /home/tghemani/repo/ibd_ibs/grm/scripts/job_reports/
#PBS -J 1-1000


set -e

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAY_INDEX=${1}
fi

i=${PBS_ARRAY_INDEX}

cifile="~/repo/eQTL-2D/data/supFile3_K562_interactingLoci_clusters.csv"
sigfile="~/repo/eQTL-2D/analysis/interaction_list_replication_summary.RData"
n=549
win=5000000
outroot="~/repo/eQTL-2D/analysis/chromosome_interactions/results/out"

R --no-save --args ${i} ${n} ${win} ${cifile} ${sigfile} ${outroot} < ~/repo/eQTL-2D/analysis/chromosome_interactions/simulation_2.R
