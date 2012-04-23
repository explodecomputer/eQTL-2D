#///=============================================================\\\#
#===================================================================#
#																	#
#	FINAL_STAGE														#
#                                                               	#
# 	Name: Make_plink_data_final.R		           					#
#                                                               	#
# 	Read in normalised data and sample information from				#
# 	"Data_files/Final_stage_data_files" and produce plink 	    	#
# 	input files	for an eQTL analysis								#
#                                                               	#
#                                                               	#
# 	Joseph Powell     17/04/2012                                	#
#                                                               	#
#===================================================================#
#\\\=============================================================///#

setwd("/Volumes/group_wrayvisscher/josephP/Projects/Genomics/Expression_main/Data/Data_analysis/Data_files/Final_stage_data_files")

#/--- Read in data - N100 ---\#

sample_info <- read.csv("sample_info_clean_S2.csv", header=T)
probe_signal_N100 <- read.csv("probe_signal_N100.csv", header=T)
probe_info_N100 <- read.csv("probe_info_N100.csv", header=T)

#///=============================================================\\\#
# Two individuals with expression data do not have genotype data 	#
# These are individuals 8624503 and 8667102. Remove from files 		#
# They correspond to rows / cols 829 and 860						#
#\\\=============================================================///#

index <- as.numeric(c(which(sample_info[,1]=="8624503"), which(sample_info[,1]=="8667102"))) 
 

sample_info <- sample_info[-index,]
probe_signal_N100 <- probe_signal_N100[,-index]


#/--- Turn missing NA's to -9 ---\#
probe_signal_N100[is.na(probe_signal_N100)] <- -9


#///=============================================================\\\#
# 			GENERATE AND WRITE OUT THE PLINK FILES					#
#\\\=============================================================///#

#/---\# individual ids

indi_ids_plink <- as.data.frame(cbind(substr(sample_info[,1], 1, 5), sample_info[,1]))


#/---\# Covariance file

cov_file_plink <- as.data.frame(cbind(indi_ids_plink, sample_info$sex, sample_info$Generation))

#/---\# Phenotype file

probe_signal_N100_plink <- cbind(indi_ids_plink, t(probe_signal_N100))
probe_signal_N100_plink <- as.data.frame(probe_signal_N100_plink)
names(probe_signal_N100_plink)[1:2] <- c("FID","IID")
names(probe_signal_N100_plink)[3:ncol(probe_signal_N100_plink)] <- as.character(probe_info_N100[,3])


#/---\# Test phenotype file - single gene expression 

test_probe_pheno_N100_plink <- probe_signal_N100_plink[,1:12]


#/---------------------\#
#-----------------------#
#	Write out files	#
#-----------------------#
#\---------------------/#
setwd("/Volumes/group_wrayvisscher/josephP/Projects/Genomics/Expression_main/Main - eQTL project/Final_data_plink_eQTL")

write.table(indi_ids_plink, "ids_plink_final.txt", quote=F, row.names=F, col.names=F, sep=" ")

write.table(cov_file_plink, "covariates_plink_final.txt", quote=F, row.names=F, col.names=F, sep=" ")

write.table(probe_signal_N100_plink, "probe_pheno_plink_final.txt", quote=F, row.names=F, col.names=T, sep=" ")


#/--- Test files ---\#

write.table(test_probe_pheno_N100_plink, "probe_pheno_plink_final_test.txt", quote=F, row.names=F, col.names=T, sep=" ")




