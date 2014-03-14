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


# Clean Westra data (kostya)

westra <- read.table(pipe("gunzip -c 2012-12-21-CisAssociationsProbeLevelFDR0.5.zip"),
                                                             header = TRUE, sep = "\t")
ill_orig <- read.delim(pipe("gunzip -c humanht-12_v3_0_r3_11283641_a_txt.zip"),
                                                        skip = 8, header = TRUE)
#------------------------------------------------------------------------------
# Educated guess that Westra et al used Array_Address_Id as probe id instead of ILMN_... ids
ill_ids <- ill_orig$Probe_Id
names(ill_ids) <- ill_orig$Array_Address_Id
westra$ill_id <- ill_ids[as.character(westra$ProbeName)]
#------------------------------------------------------------------------------
c_as.num <- plyr::colwise(as.numeric)
z_val <- stringr::str_split_fixed(westra$DatasetsZScores, ",", n = 9)
z_val <- c_as.num(as.data.frame(z_val))
colnames(z_val) <- c("EGCUT","SHIP_TREND","Groningen-HT12",
                     "Groningen-H8v2","Rotterdam","DILGOM",
                     "INCHIANTI","HVH-HT12v3", "HVH-HT12v4")
westra <- cbind(westra, z_val)



# Three snp haplotype analysis


three_snp_analysis <- three_snp_test.fun(info, block, bsgs)
write.csv(three_snp_analysis, "/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/Frayling_et_al/three_snp_analysis.csv", quote=F, row.names=F)


#=======================================================#
#			PREDICTING GENOTYPES 						#
#=======================================================#

ld_info <- read.csv



out <- array(0, c(nrow(info), 3))
for(i in 1:nrow(info)) {
	if(i==6) {
		out[i,] <- NA

	}
	else{
		snp1  <- block[,which(colnames(block)==as.character(info$SNP1[i]))] 
		snp2  <- block[,which(colnames(block)==as.character(info$SNP2[i]))] 
		inc_snp <- block[,which(colnames(block)==as.character(info$rs_id[i]))] 
		
		out[i,1] <- round(summary(lm(inc_snp ~ snp1 + snp2 + snp1:snp2))$adj.r.squared,2)
		out[i,2] <- round(cor(inc_snp, snp1)^2,2)
		out[i,3] <- round(cor(inc_snp, snp2)^2,2)
	}

}
colnames(out) <- c("r2_full", "r2_snp1", "r2_snp2")
out <- cbind(info, out)

write.csv(out, "Fraylingsnp_prediction.csv")














