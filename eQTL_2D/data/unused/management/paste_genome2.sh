#!/bin/bash


# first run makemapfiles.R

id=${1}
filename="GWAS_Adolescent_chr${id}"
out="chr${id}"

cp ../../${filename}.ped.gz .
gunzip ${filename}.ped.gz
mv ${filename}.ped ${out}.ped
mv GWAS_all.map${id} ${out}.map

../plink --file ${out} --make-bed --out ${out}

