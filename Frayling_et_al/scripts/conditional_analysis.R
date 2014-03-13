#=======================================================#
#-------------------------------------------------------#
#														#
#	conditional_analysis.R								#
#														#
#	conditional analysis using GCTA					 	#
#	after the snps in westra are fitted					#
#	list.												#
#														#
#	joseph.powell@uq.edu.au								#
#	V1. March.2014										#
#														#
#-------------------------------------------------------#
#=======================================================#

source("/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/Frayling_et_al/scripts/fine_mapping_functions.R")

#=======================================================#
#				READ IN THE DATA 						#
#=======================================================#


info <- read.csv("/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/Frayling_et_al/data_files/inc_info_bsgs.csv", header=T)
load("/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/data/bsgs_egcut_fehr_data.RData")
bsgs <- phenlist[[1]]

ped <- read.table("/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/Frayling_et_al/test/snp_list.ped", header=F)
map <- read.table("/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/Frayling_et_al/test/snp_list.map", header=F)
block <- plink_to_012.fun(ped, map)



# file mapped by Kostya
westra <- 

# read in 





three_snp_test.fun <- function(info, block, bsgs){

	out <- array(0, c(nrow(info), 4)) 

	for(i in 1:nrow(info)) {
		
		if(i == 6) {

			out[i,] <- NA
		}

		else {
			snp1  <- block[,which(colnames(block)==as.character(info$SNP1[i]))] 
			snp2  <- block[,which(colnames(block)==as.character(info$SNP2[i]))] 
			inc_snp <- block[,which(colnames(block)==as.character(info$rs_id[i]))] 

			pheno <- bsgs[,which(colnames(bsgs)==as.character(info$Probe[i]))]
			pheno_adj <- summary(lm(pheno ~ inc_snp))$residuals


			fullmod <- lm(pheno ~ as.factor(inc_snp) + as.factor(snp1) + as.factor(snp2) + as.factor(snp1):as.factor(snp2):as.factor(inc_snp))
			redmod <- lm(pheno ~ as.factor(snp1) + as.factor(snp2) + as.factor(inc_snp))
			# This is the interaction terms on their own (nested test)
			intmod <- anova(redmod, fullmod)	
			# Extract statistics	
			tmp <- summary(fullmod)$fstatistic
			out[i,1] <- round(-log10(pf(tmp[1], tmp[2], tmp[3], low=F)),2)		
			out[i,2] <- (summary(fullmod))$fstatistic[2]	
			out[i,3] <- round(-log10(intmod$Pr[2]), 2)
			out[i,4] <- 
		}
	}

	colnames(out) <- c("")






