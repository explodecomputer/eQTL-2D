# tabulation of ld and maf info


info <- read.csv("../data_files/inc_info_bsgs.csv", header=T)
maf <- read.table("maf_info.txt", header=T)
ld <- read.table("ld_info.ld", header=F)
snps <- read.table("snp_list.map", header=F)

row.names(ld) <- snps$V2
names(ld) <- snps$V2



out <- array(0, c(nrow(info), 6))
for(i in 1:nrow(info)) {

	out[i,1] <- maf$MAF[which(maf$SNP == as.character(info$SNP1[i]))]
	out[i,2] <- maf$MAF[which(maf$SNP == as.character(info$SNP2[i]))]
	tmp <- maf$MAF[which(maf$SNP == as.character(info$rs_id[i]))]
	if(length(tmp)==0) {
		out[i,3] <- NA
	}
	if(length(tmp)!=0) {
		out[i,3] <- tmp
	}

	s1 <- which(row.names(ld)==as.character(info$SNP1[i]))
	s2 <- which(row.names(ld)==as.character(info$SNP2[i]))
	s3 <- which(row.names(ld)==as.character(info$rs_id[i])) 

	out[i,4] <- ld[s1,s2]

	if(length(s3)==0) {
		out[i,5] <- NA
		out[i,6] <- NA
	}

	if(length(s3)!=0) {
		out[i,5] <- ld[s1,s3]
		out[i,6] <- ld[s2,s3]
	}

}


out <- round(out,3)
out <- as.data.frame(out) 
names(out) <- c("snp1_maf", "snp2_maf", "inc_snp_maf", "r2_snp1_snp2", "r2_snp1_inc_snp", "rs_snp2_inc_snp") 

out <- cbind(info, out)
head(out)
write.csv(out, "ld_maf_bsgs_info.csv", quote=F, row.names=F)


getwd()
