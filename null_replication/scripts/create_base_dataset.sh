#!/bin/bash

# The original dataset is here on bc3 here:
# /panfs/panasas01/shared/alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/hrc/released/2017-05-04/data/plink/combined/combined
# Copied across to bc4 because no mirror there yet

datafile="../data/combined"
snplist="../../data/hg19locations.bim"
idlist="../data/inclusion_list.txt"

disc="../data/disc"
rep="../data/rep"

ndisc=846
nrep=2131
total=$((ndisc + nrep))
echo $total


# Get individual list

head -n $ndisc $idlist > ../data/disclist.txt
head -n $total $idlist | tail -n $nrep > ../data/replist.txt
wc -l ../data/disclist.txt
wc -l ../data/replist.txt

# Get SNP list

cat <(cut -d " " -f 2 $snplist) | sort -u > ../data/base_snplist.txt

# Create discovery and replication datasets

plink --bfile $datafile --extract ../data/base_snplist.txt --keep ../data/disclist.txt --make-bed --out ../data/disc

plink --bfile $datafile --extract ../data/base_snplist.txt --keep ../data/replist.txt --make-bed --out ../data/rep
