# imputed analysis




info <- read.csv("/fileserver/group/wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/Frayling_et_al/data_files/inc_info_bsgs.csv", header=T)

lf <- list.files("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/epi_scan_imput_outputs_2/")


out <- array(0, c(nrow(info)))
for(i in 1:nrow(info)) {

	a <- which(lf==paste("epi_impute_scan_", as.character(info$Probe[i]), "_", as.character(info$SNP1[i]), "_", as.character(info$SNP2[i]), ".txt", sep=""))
	if(length(a)>0) {
		out[i] <- a
	} 

}

info2 <- info[which(out>0),]

a1 <- NULL
a2 <- NULL
a3 <- NULL
for(i in 1:nrow(info2)) {

	tmp <- read.table(paste("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Epistasis/epi_scan_imput_outputs_2/epi_impute_scan_", as.character(info2$Probe[i]), "_", as.character(info2$SNP1[i]), "_", as.character(info2$SNP2[i]), ".txt", sep=""), header=T)

	a <- tmp[which(tmp$snp1==as.character(info2$SNP1[i]) & tmp$snp2==as.character(info2$SNP2[i])),]
	b <- tmp[which(tmp$snp1==as.character(info2$rs_id[i]) & tmp$snp2==as.character(info2$SNP2[i])),]
	c <- tmp[which(tmp$snp1==as.character(info2$SNP1[i]) & tmp$snp2==as.character(info2$rs_id[i])),]
	a1 <- rbind(a1, a)
	a2 <- rbind(a2, b)
	a3 <- rbind(a3, c)
}



