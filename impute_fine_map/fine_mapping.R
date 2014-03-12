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
#load("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/residuals.RData")
load("/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/data/bsgs_egcut_fehr_data.RData")
bsgs <- phenlist[[1]]


#=======================================================#
#		MAKE THE SNP BLOCKS FROM WHOLE DATASET 			#
#=======================================================#


make_snp_region_files.fun(info)

#=======================================================#
#		RUN THE EPISCAN ON THE IMPUTED DATA 			#
#=======================================================#

blockdir <- "/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/"
outdir <- "/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D//impute_fine_map/frayling_results/"

for(i in 1:nrow(info)) {
	ped1 <- read.table(paste(blockdir, as.character(info$SNP1[i]), ".ped", sep=""), header=F)
	ped2 <- read.table(paste(blockdir, as.character(info$SNP2[i]), ".ped", sep=""), header=F)
	 
	map1 <- read.table(paste(blockdir, as.character(info$SNP1[i]), ".map", sep=""), header=F)
	map2 <- read.table(paste(blockdir, as.character(info$SNP2[i]), ".map", sep=""), header=F)

	block1 <- plink_to_012.fun(ped1, map1)
	block2 <- plink_to_012.fun(ped2, map2)

	pheno <- bsgs[,which(colnames(bsgs)==as.character(info$Probe[i]))]

	epi_scan_out <- epi_scan.fun(block1, block2, pheno)
	write.csv(epi_scan_out, paste(outdir, as.character(info$Probe[i]), "_", as.character(info$SNP1[i]), "_", as.character(info$SNP2[i]), ".csv", sep=""), quote=F, row.names=F)
}





#=======================================================#
#		ANALYSIS OF THE IMPUTED REGION EPI SCAN 		#
#=======================================================#

load("/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/analysis/interaction_list_meta_analysis.RData")

out <- array(0, c(nrow(info), 18))

for(i in 1:nrow(info)) {

	m <- which(meta$probename==as.character(info$Probe[i]) & meta$snp1==as.character(info$SNP1[i]) & meta$snp2==as.character(info$SNP2[i]))
	out[i,1:6] <- as.matrix(meta[m,c(1,2, 8, 9, 11, 12)])

	tmp <- read.csv(paste(outdir, as.character(info$Probe[i]), "_", as.character(info$SNP1[i]), "_", as.character(info$SNP2[i]), ".csv", sep=""), header=T)
	tmp$probe <- as.character(info$Probe[i])
	filter <- which(tmp$rsq < 0.2 & tmp$nclass==9 & tmp$minclass > 5)
	tmp2 <- tmp[filter,]
	hit <- tmp2[which.max(tmp2$intP),]
	out[i,8:17] <- as.matrix(hit)

}



out <- as.data.frame(out)
names(out)[8:17] <- names(hit)
names(out)[18] <- "R2_new"
names(out)[7] <- "R2_original"
names(out)[1:6] <- c("snp1", "snp2", "probe", "gene", "8df_p_original", "4df_p_original")
write.csv(out, "/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/Frayling_et_al/data/fine_mapping.csv", quote=F, row.names=T)

lf <- list.files()













#=======================================================#
#		ANALYSIS ADJUSTING ON THE INC SNP ADDITIVE 		#
#=======================================================#


blockdir <- "/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/snp_region_data/"
outdir <- "/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/impute_fine_map/frayling_results/"

ped3 <- read.table("/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/Frayling_et_al/test/snp_list.ped", header=F)
map3 <- read.table("/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/Frayling_et_al/test/snp_list.map", header=F)
block3 <- plink_to_012.fun(ped3, map3)

for(i in 15:nrow(info)) {
	ped1 <- read.table(paste(blockdir, as.character(info$SNP1[i]), ".ped", sep=""), header=F)
	ped2 <- read.table(paste(blockdir, as.character(info$SNP2[i]), ".ped", sep=""), header=F)
	
	map1 <- read.table(paste(blockdir, as.character(info$SNP1[i]), ".map", sep=""), header=F)
	map2 <- read.table(paste(blockdir, as.character(info$SNP2[i]), ".map", sep=""), header=F)

	block1 <- plink_to_012.fun(ped1, map1)
	block2 <- plink_to_012.fun(ped2, map2)

	pheno <- bsgs[,which(colnames(bsgs)==as.character(info$Probe[i]))]
	
	inc_snp <- which(colnames(block3)==as.character(info$rs_id[i]))
	
	if(length(inc_snp)==0) {
		print("snp missing ", as.character(info$rs_id[i]))
	}	

	if(length(inc_snp)!=0) {
		adj_pheno <- summary(lm(pheno~block3[,inc_snp]))$residuals
		epi_scan_out <- epi_scan.fun(block1, block2, adj_pheno)
		write.csv(epi_scan_out, paste(outdir, as.character(info$Probe[i]), "_", as.character(info$SNP1[i]), "_", as.character(info$SNP2[i]), "adj.csv", sep=""), quote=F, row.names=F)
	}
	print(i)
}





