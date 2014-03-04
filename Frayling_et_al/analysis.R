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
n <- 200
out2 <- array(0, c(nrow(pval), ncol(pval)))
out2[,1] <- round((out[,1]/n)/(1+(out[,1]/n))*100, 2)
out2[,2] <- round((out[,2]/n)/(1+(out[,2]/n))*100, 2)
out2[,3] <- round((out[,3]/n)/(1+(out[,3]/n))*100, 2)
out2[,4] <- round((8*out[,4]/n)/(1+8*(out[,4]/n))*100, 2)

out_all <- cbind(out, out2)
write.csv(out_all, "/Users/jpowell/repo/eQTL-2D/Frayling_et_al/data_files/chisq_var.csv", quote=F, row.names=F)




# get additive effects of their inc snps
info <- read.csv("Frayling_et_al/data_files/inc_info_bsgs.csv", header=T)


dir <- "/ibscratch/wrayvisscher/josephP/BSGS/Imputed/Var_eQTL/data/Output_Z/"

s <- array(0, c(nrow(info), 11))
for(i in 1:nrow(info)) {

	if(!is.na(info$chr[i])) {

		if(nchar(info$chr[i])==1) {
			tmp <- read.table(paste(dir, as.character(info$Probe[i]), "/", as.character(info$Probe[i]), "_Z-fastassoc-chr0", as.character(info$chr[i]), ".tbl", sep=""), header=T)
			out <- tmp[which(tmp$SNP==as.character(info$rs_id[i])),]
			s[i,] <- as.matrix(out)	
			
		}

		if(nchar(info$chr[i])==2) {
			tmp <- read.table(paste(dir, as.character(info$Probe[i]), "/", as.character(info$Probe[i]), "_Z-fastassoc-chr", as.character(info$chr[i]), ".tbl", sep=""), header=T)
			out <- tmp[which(tmp$SNP==as.character(info$rs_id[i])),]
			s[i,] <- as.matrix(out)	
		
		}
		print(i)
	}

	else {

		

	}	
}

s <- as.data.frame(s)
names(s) <- names(tmp)
write.csv(s, "Frayling_et_al/data_files/inc_snps_in_bsgs.csv", quote=F, row.names=F)




# extract sig results

index <- array(0, c(nrow(info)))
for(i in 1:nrow(info)) {
	tmp <- which(sig$probename==as.character(info$Probe[i]) & sig$snp1==as.character(info$SNP1[i]) & sig$snp2==as.character(info$SNP2[i]))
	if(length(tmp)>0) {

		index[i] <- tmp
	}


}










