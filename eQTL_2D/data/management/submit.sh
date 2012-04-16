#!/bin/bash

#$ -N datasort
#$ -cwd
#$ -t 1-20

./paste_genome2.sh ${SGE_TASK_ID}

