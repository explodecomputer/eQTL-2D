# incseq analysis


inc <- read.table("/fileserver/group/wrayvisscher/josephP/incseq.txt", header=T)
snp <- read.table("/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Clean_R2_Imputed_BSGS_data/bsgs_imputed_R2_80_cleaned_stage2_chr_all_SNP_info.txt", header=T)

inc_split <- strsplit(as.character(inc$incseq), ":")


out <- array(0, c(length(inc_split)))
for(i in 1:length(inc_split)) {
	a <- which(snp$chr==inc_split[[i]][1] & snp$position==inc_split[[i]][2])
	if(length(a)!=0) {
		out[i] <- a
	}

	print(i)
}

snp_info <- snp[out,]
tmp <- cbind(snp_info, inc)
write.csv(tmp, "inc_bsgs_info.csv", quote=F, row.names=F)


# extract bsgs additive effect info 
inc_info <- read.csv("/fileserver/group/wrayvisscher/josephP/inc_info.csv", header=T)

# read in results




# pvalues
pval <- read.csv("/Users/jpowell/repo/eQTL-2D/Frayling_et_al/data_files/pval.csv", header=F)

# chisq
out <- array(0, c(nrow(pval), ncol(pval)))
out[,1] <- qchisq(pval[,1], df=1, lower.tail=F)
out[,2] <- qchisq(pval[,2], df=1, lower.tail=F)
out[,3] <- qchisq(pval[,3], df=1, lower.tail=F)
out[,4] <- qchisq(pval[,4], df=8, lower.tail=F)

# var explained
out2 <- array(0, c(nrow(pval), ncol(pval)))
out2[,1] <- round((out[,1]/500)/(1+(out[,1]/500))*100, 2)
out2[,2] <- round((out[,2]/500)/(1+(out[,2]/500))*100, 2)
out2[,3] <- round((out[,3]/500)/(1+(out[,3]/500))*100, 2)
out2[,4] <- round((8*out[,4]/500)/(1+8*(out[,4]/500))*100, 2)

out_all <- cbind(out, out2)
write.csv(out_all, "/Users/jpowell/repo/eQTL-2D/Frayling_et_al/data_files/chisq_var.csv", quote=F, row.names=F)