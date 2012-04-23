#!bin/sh

~/../apps/plink/plink-1.07-x86_64/plink --bfile /hpscratch/wrayvisscher/GWAS/GWAS/GeneralRelease/Release6/GWAS --nonfounders --keep ids_plink_final.txt --maf 0.01 --recode --allele1234 --make-bed --out clean_geno_final


