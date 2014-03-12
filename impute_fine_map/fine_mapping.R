#=======================================================#
#-------------------------------------------------------#
#														#
#	fine_mapping.R	 									#
#														#
#	run the analysis to fine map the interactions for 	#
#	the list of probes and snps given in the wood et al	#
#	list.												#
#														#
#	joseph.powell@uq.edu.au								#
#	V1. March.2014										#
#														#
#-------------------------------------------------------#
#=======================================================#

#=======================================================#
#		READ IN THE FILTERED DATASETS (FROM GIB)		#
#=======================================================#


#load("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/agg_filtered.RData")

info <- read.csv("/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/Frayling_et_al/data_files/inc_info_bsgs.csv", header=T)
load("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/residuals.RData")


#=======================================================#
#		MAKE THE SNP BLOCKS FROM WHOLE DATASET 			#
#=======================================================#


make_snp_region_files.fun(info)

#=======================================================#
#		RUN THE EPISCAN ON THE IMPUTED DATA 			#
#=======================================================#


for


blockdir <- "/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/"

ped1 <- read.table(paste(blockdir, as.character(info$SNP1[i]), ".ped", sep=""), header=F)
ped2 <- read.table(paste(blockdir, as.character(info$SNP2[i]), ".ped", sep=""), header=F)
 
map1 <- read.table(paste(blockdir, as.character(info$SNP1[i]), ".map", sep=""), header=F)
map2 <- read.table(paste(blockdir, as.character(info$SNP2[i]), ".map", sep=""), header=F)

block1 <- plink_to_012.fun(ped1, map1)
block2 <- plink_to_012.fun(ped2, map2)

pheno <- resphen[,which(colnames(resphen)==as.character(info$Probe[i]))]

epi_scan_out <- epi_scan.fun(block1, block2, pheno)



epi_scan_out[which.max(epi_scan_out$intP),]

#


snp_info <- read.table("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Var_eQTL/bsgs_imputed_R2_80_cleaned_stage2_chr_all_SNP_info.txt", header=T)
snp_info_geno <- read.table("/ibscratch/wrayvisscher/josephP/BSGS/Genotype/GWAS.map", header=T)
# pheno
load("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/residuals.RData")


