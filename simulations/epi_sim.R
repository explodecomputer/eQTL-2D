#=======================================================#
#-------------------------------------------------------#
#														#
#	epi_sim.R 											#
#														#
#	Run the analysis using the functioms in epi_sim		#
#	_functions.R 										#
#														#
#	joseph.powell@uq.edu.au								#
#	V1. March.2013										#
#														#
#-------------------------------------------------------#
#=======================================================#


#=======================================================#
#				READ IN THE FILES				 		#
#=======================================================#


n <- as.numeric(commandArgs(T)[1])
pheno <-  commandArgs(T)[2]
geno <- commandArgs(T)[3]
outdir <- commandArgs(T)[4]
blockdir <- commandArgs(T)[5]
g_info <- commandArgs(T)[6]
i_info <- commandArgs(T)[7]
type <- commandArgs(T)[8]
fun <- commandArgs(T)[9]

load(pheno)
load(geno)

geno_info <- read.table(g_info, header=T)
impute_info <- read.table(i_info, header=T)

source(fun)

#=======================================================#
#		OBTAIN CLEAN SETS OF SNP AND PROBE DATA 		#
#=======================================================#


probe_info <- probeinfo[which(probeinfo$CHROMOSOME_NEW >=1 & probeinfo$CHROMOSOME_NEW < 23),]
snps <- colnames(geno)
index <- which(snps %in% impute_info$rs_id)
snps <- snps[index]
index <- match(snps, geno_info$rs_id)
geno_info <- geno_info[index,]


#=======================================================#
#			SELECT THE SNP PAIR AND PROBE 		  		#
#=======================================================#

pick_out <- snp_pick.fun(type, probe_info, geno_info)


#=======================================================#
#		OBTAIN CLEAN SETS OF SNP AND PROBE DATA 		#
#=======================================================#

# Data surrounding SNP 1
snp1 <- as.character(pick_out[[1]]$rs_id[1])	
system(paste("/clusterdata/apps/plink/plink-1.07-x86_64/plink --bfile /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Clean_R2_Imputed_BSGS_data/bsgs_imputed_R2_80_cleaned_stage2_chr_all --snp ", snp1, " --window 100 --recode12 --out /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/sim_impute_regions/", snp1,  sep=""))	

# Data surrounding SNP 2
snp2 <- as.character(pick_out[[1]]$rs_id[2])	
system(paste("/clusterdata/apps/plink/plink-1.07-x86_64/plink --bfile /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Clean_R2_Imputed_BSGS_data/bsgs_imputed_R2_80_cleaned_stage2_chr_all --snp ", snp2, " --window 100 --recode12 --out /ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/sim_impute_regions/", snp2,  sep=""))	



#=======================================================#
#			CONVERT THE IMPUTED GENOTYPES	  			#
#=======================================================#

ped1 <- read.table(paste(blockdir, snp1, ".ped", sep=""), header=F)
ped2 <- read.table(paste(blockdir, snp2, ".ped", sep=""), header=F)
 
map1 <- read.table(paste(blockdir, snp1, ".map", sep=""), header=F)
map2 <- read.table(paste(blockdir, snp2, ".map", sep=""), header=F)

block1 <- plink_to_012.fun(ped1, map1)
block2 <- plink_to_012.fun(ped2, map2)

pheno <- resphen[,which(colnames(resphen)==as.character(pick_out[[2]]$PROBE_ID))]


#=======================================================#
#	RUN THE 4 AND 8DF SCAN ACROSS THE IMPUTED REGION	#
#=======================================================#

epi_scan_out <- epi_scan.fun(block1, block2, pheno)

#=======================================================#
#		CHOOSE THE RELEVENT (IE SMALLEST) PVALUES		#
#=======================================================#

# (Genoed) SNP
i1 <- which(epi_scan_out$snp1==snp1 & epi_scan_out$snp2==snp2)
imputed_geno_snps <- epi_scan_out[i1,]


# Max p-value
filter <- which(epi_scan_out$nclass==9 & as.numeric(as.matrix(epi_scan_out$minclass))>4)
epi_scan_out <- epi_scan_out[filter,]
epi_scan_maxp <- epi_scan_out[which.max(epi_scan_out$intP),]


#=======================================================#
#		PERFOM THE 4DF AND 8DF ASSOCS ON GENO SNPS		#
#=======================================================#

geno_snp1 <- geno[ ,which(colnames(geno)==snp1)]
geno_snp2 <- geno[ ,which(colnames(geno)==snp2)]

epi_geno_out <- epi_geno.fun(geno_snp1, geno_snp2, pheno, snp1, snp2)


#=======================================================#
#	CALCULATE THE CORRELATION BETWEEN GENO AND IMPUTE 	#
#=======================================================#

impute_snp1 <- block1[ ,which(colnames(block1)==snp1)]
impute_snp2 <- block2[ ,which(colnames(block2)==snp2)]

geno_impute_snp1_cor <- abs(cor(impute_snp1, geno_snp1))
geno_impute_snp2_cor <- abs(cor(impute_snp2, geno_snp2))


#=======================================================#
#		COMBINE ALL THE RESULTS AND WRITE OUT 			#
#=======================================================#

# Results
save(pick_out, geno_impute_snp1_cor, geno_impute_snp2_cor, epi_geno_out, imputed_geno_snps, epi_scan_maxp, file=paste(outdir, "results_", type, "_sim_", n, ".RData", sep=""))

# relevent data
save(pick_out, block1, block2, pheno, epi_scan_out, file=paste(outdir, "data_", type, "_sim_", n, ".RData", sep=""))




