#!/bin/bash

# Get list of SNPs
cat <(cut -d " " -f 2 ../data/finemap.txt) <(cut -d " " -f 4 ../data/finemap.txt) | sort -u > ../data/snplist.txt

wc -l ../data/snplist.txt

# Get LD between all SNPs
plink --bfile ../data/combined --extract ../data/snplist.txt --make-bed --out ../data/snplist
plink --bfile ../data/snplist --r2 square gz --out ../data/snplist
awk '{ print $2 }' ../data/snplist.bim > ../data/snplist.snps
wc -l ../data/snplist.snps
head ../data/snplist.snps
head ../data/snplist.txt

