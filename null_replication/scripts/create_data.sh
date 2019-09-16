#!/bin/bash

# The original dataset is here on bc3 here:
# /panfs/panasas01/shared/alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/hrc/released/2017-05-04/data/plink/combined/combined
# Copied across to bc4 because no mirror there yet

datafile="../data/combined"
sensnp="rs67903230"
cissnp="rs13069559"
varexp="0.5"

disc="../data/disc"
rep="../data/rep"

cischr=`grep -w $cissnp $disc.bim | cut -f 1`
echo $cischr

sd="../data/${sensnp}"
mkdir -p ${sd}

# Create epigpu format
# cis SNP = chr1
# trans SNP = chr2
plink --bfile ${disc} --snps ${cissnp} --make-bed --out ${sd}/${cissnp}_disc
plink --bfile ${rep} --snps ${cissnp} --make-bed --out ${sd}/${cissnp}_rep

plink --bfile ${disc} --not-chr ${cischr} --make-bed --out ${sd}/${cissnp}_exclude_disc
plink --bfile ${rep} --not-chr ${cischr} --make-bed --out ${sd}/${cissnp}_exclude_rep


awk -v OFS='\t' '{ print "1", $2, $3, $4, $5, $6 }' ${sd}/${cissnp}_disc.bim > ${sd}/${cissnp}_disc.bim2 && mv ${sd}/${cissnp}_disc.bim2 ${sd}/${cissnp}_disc.bim
awk -v OFS='\t' '{ print "1", $2, $3, $4, $5, $6 }' ${sd}/${cissnp}_rep.bim > ${sd}/${cissnp}_rep.bim2 && mv ${sd}/${cissnp}_rep.bim2 ${sd}/${cissnp}_rep.bim

awk -v OFS='\t' '{ print "2", $2, $3, $4, $5, $6 }' ${sd}/${cissnp}_exclude_disc.bim > ${sd}/${cissnp}_exclude_disc.bim2 && mv ${sd}/${cissnp}_exclude_disc.bim2 ${sd}/${cissnp}_exclude_disc.bim
awk -v OFS='\t' '{ print "2", $2, $3, $4, $5, $6 }' ${sd}/${cissnp}_exclude_rep.bim > ${sd}/${cissnp}_exclude_rep.bim2 && mv ${sd}/${cissnp}_exclude_rep.bim2 ${sd}/${cissnp}_exclude_rep.bim


# Create discovery dataset
plink --bfile ${sd}/${cissnp}_disc --bmerge ${sd}/${cissnp}_exclude_disc --recode --out ${sd}/${cissnp}_merge_disc
episcan -Drm ${sd}/${cissnp}_merge_disc.ped ${sd}/${cissnp}_merge_disc.map ${sd}/${cissnp}_merge_disc.egu


# Create phenotypes using sentinel SNP
plink --bfile ${datafile} --snps ${sensnp} --keep ${disc}list.txt --make-bed --out ${sd}/${sensnp}_disc
plink --bfile ${datafile} --snps ${sensnp} --keep ${rep}list.txt --make-bed --out ${sd}/${sensnp}_rep

plink --bfile ${sd}/${sensnp}_disc --recode A --out ${sd}/${sensnp}_disc
Rscript makephen.r ${sd}/${sensnp}_disc.raw ${sd}/${sensnp}_disc.fam ${varexp}

plink --bfile ${sd}/${sensnp}_rep --recode A --out ${sd}/${sensnp}_rep
Rscript makephen.r ${sd}/${sensnp}_rep.raw ${sd}/${sensnp}_rep.fam ${varexp}


# Perform scan
episcan -A ${sd}/${cissnp}_merge_disc.egu ${sd}/${cissnp}_disc_out -t i -F 10000 -I 11 -1 1 -2 2 -f ${sd}/${sensnp}_disc.fam

cat ${sd}/${cissnp}_disc_out[0-9]* > ${sd}/${cissnp}_disc_out
rm ${sd}/${cissnp}_disc_out[0-9]*

Rscript formatres.r ${sd}/${cissnp}_disc_out ${sd}/${cissnp}_merge_disc.map ${sd}/${cissnp}_tophits.txt 


# Perform replication
plink --bfile ${sd}/${cissnp}_exclude_rep --extract ${sd}/${cissnp}_tophits.txt --make-bed --out ${sd}/${cissnp}_exclude_rep
plink --bfile ${sd}/${cissnp}_rep --bmerge ${sd}/${cissnp}_exclude_rep --recode --out ${sd}/${cissnp}_merge_rep
episcan -Drm ${sd}/${cissnp}_merge_rep.ped ${sd}/${cissnp}_merge_rep.map ${sd}/${cissnp}_merge_rep.egu
episcan -A ${sd}/${cissnp}_merge_rep.egu ${sd}/${cissnp}_rep_out -t i -F 10000 -I 0.0001 -1 1 -2 2 -f ${sd}/${sensnp}_rep.fam
cat ${sd}/${cissnp}_rep_out[0-9]* > ${sd}/${cissnp}_rep_out





