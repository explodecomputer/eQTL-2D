#!/bin/bash

#$ -N dataploink
#$ -cwd

../plink --merge-list merge.txt --make-bed --out eqtl2dall

