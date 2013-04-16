#!/bin/bash

#$ -N colall
#$ -cwd -V
#$ -S /bin/bash

start=301
end=400

/clusterdata/apps/R-2.14/bin/R --no-save --args cols ${start} ${end} chunks < collate_chunks.R

