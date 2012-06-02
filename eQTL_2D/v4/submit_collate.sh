#!/bin/bash

#$ -N collate
#$ -cwd -V
#$ -S /bin/bash
#$ -t 1-100

#SGE_TASK_ID=${1}
i=${SGE_TASK_ID}
n=20
start=$(((i-1)*n+1))
end=$((n*i))

/clusterdata/apps/R-2.14/bin/R --no-save --args results/result scratch/resphen scratch/resphen ${start} ${end} 12 ../data/eqtl2d_objects.RData res${start}-${end}.RData < collate_results.R



