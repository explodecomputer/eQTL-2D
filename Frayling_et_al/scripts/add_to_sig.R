


out <- array(0, c(nrow(info)*2, ncol(sig)))
head(info)


k <- 0
for(i in seq(1, nrow(out), 2)) {
	k <- k+1	

	# probe name add
	out[c(i, i+1),13] <- as.character(info$Probe[k])

	# Snp 1 vs inc snp add
	out[i, 5] <- as.character(info$SNP1[k])
	out[i, 6] <- as.character(info$rs_id_inc_in_bsgs[k])

	# Snp 2 vs inc snp add
	out[i+1, 5] <- as.character(info$SNP2[k])
	out[i+1, 6] <- as.character(info$rs_id_inc_in_bsgs[k])

	a <- which(snp_info$rs_id==as.character(info$SNP1[k]))
	if(length(a)==1) {
		out[i, 1] <- snp_info$chr[a]	
		out[i, 3] <- snp_info$position[a]
	}

	b <- which(snp_info$rs_id==as.character(info$SNP2[k]))
	if(length(b)==1) {
		out[i+1, 1] <- snp_info$chr[b]	
		out[i+1, 3] <- snp_info$position[b]
	}

	c <- which(snp_info$rs_id==as.character(info$rs_id_inc_in_bsgs[k]))
	if(length(c)==1) {
		out[i, 2] <- snp_info$chr[c]	
		out[i, 4] <- snp_info$position[c]

		out[i+1, 2] <- snp_info$chr[c]	
		out[i+1, 4] <- snp_info$position[c]

	}

	# add gene info
	out[c(i,i+1), 15] <- as.character(info$GENE[k])

}

head(sig)
head(out, 10)
head(info)