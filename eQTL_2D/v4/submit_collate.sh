#!/bin/bash

#$ -N collate
#$ -cwd -V
#$ -S /bin/bash
#$ -t 1-200

#SGE_TASK_ID=${1}
i=${SGE_TASK_ID}
n=10
start=$(((i-1)*n+1))
end=$((n*i))

/clusterdata/apps/R-2.14/bin/R --no-save --args results/result ../data/residuals.RData ../data/clean_geno_final.RData scratch/resphen ${start} ${end} 12 ../data/eqtl2d_objects.RData res${start}-${end}.RData < collate_results.R



