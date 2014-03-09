

# read in geno and pheno files
map <- read.table("snp_list.map", header=F)
ped <- read.table("snp_list.ped", header=F)
info <- read.csv("../data_files/inc_info_bsgs.csv", header=T)
load("../../data/bsgs_egcut_fehr_data.RData")
bsgs <- phenlist[[1]]


#=======================================================#
#		CONVERT THE PED FORMT TO A 0, 1, 2, FORMAT 		#
#=======================================================#


plink_to_012.fun <- function(
	ped, 		# ped file
	map) {		# map file

	ids <- ped[, 1:6]
	nid <- nrow(ids)
	ped <- ped[, -c(1:6)]
	index <- seq(1, ncol(ped), 2)
	geno <- matrix(0, nid, length(index))

	# Convert to 0, 1, 2 format
	for(i in 1:length(index)) {
		snp <- ped[,c(index[i], index[i]+1)]
		x <- array(NA, nid)
		snp[snp == "0"] <- NA

		i0 <- snp[,1] == 1 & snp[,2] == 1
		i2 <- snp[,1] == 2 & snp[,2] == 2
		i1 <- (snp[,1] == 1 & snp[,2] == 2) | (snp[,1] == 2 & snp[,2] == 1)
		x[i0] <- 0
		x[i1] <- 1
		x[i2] <- 2
		geno[, i] <- x
	}

	colnames(geno) <- map$V2
	rownames(geno) <- ids$V2
	return(geno)
}

block <- plink_to_012.fun(ped, map)



out <- array(0, c(nrow(info), 10))
for(i in 1:nrow(info)) {
	
	pheno <- bsgs[,which(colnames(bsgs)==as.character(info$Probe[i]))]
	snp1 <- block[,which(colnames(block)==as.character(info$SNP1[i]))]
	snp2 <- block[,which(colnames(block)==as.character(info$SNP2[i]))]
	tmp  <- which(colnames(block)==as.character(info$rs_id[i]))

	if(length(tmp)==0) {
		out[i,] <- NA	

	}
	if(length(tmp!=0)) {

		inc_snp <- block[,tmp]

		Osnp1 <- summary(lm(pheno ~ snp1))
		Osnp2 <- summary(lm(pheno ~ snp2))
		Oinc <- summary(lm(pheno ~ inc_snp))

		apheno <- summary(lm(pheno ~ inc_snp))$residuals

		Asnp1 <- summary(lm(apheno ~ snp1))
		Asnp2 <- summary(lm(apheno ~ snp2))


		out[i,1] <- Osnp1$r.squared*100
		out[i,2] <- -log10(Osnp1$coefficients[2,4])
		out[i,3] <- Osnp2$r.squared*100
		out[i,4] <- -log10(Osnp2$coefficients[2,4])
		out[i,5] <- Oinc$r.squared*100
		out[i,6] <- -log10(Oinc$coefficients[2,4])
		out[i,7] <- Asnp1$r.squared*100
		out[i,8] <- -log10(Asnp1$coefficients[2,4])
		out[i,9] <- Asnp2$r.squared*100
		out[i,10] <- -log10(Asnp2$coefficients[2,4])

	}
}

out <- round(out, 3)
out <- as.data.frame(out)
out <- cbind(out)







pheno <- resphen[,which(colnames(resphen)=="ILMN_1710752")]

adj <- lm(pheno~snp3)
pheno_adj <- adj$residuals



fullmod1 <- lm(pheno~ as.factor(snp1)+as.factor(snp2)+as.factor(snp1):as.factor(snp2))
redmod1 <- lm(pheno~ as.factor(snp1)+as.factor(snp2))
intmod1 <- anova(fullmod1, redmod1)


fullmod_adj <- lm(pheno_adj~ as.factor(snp1)+as.factor(snp2)+as.factor(snp1):as.factor(snp2))
redmod_adj <- lm(pheno_adj~ as.factor(snp1)+as.factor(snp2))
intmod_adj <- anova(fullmod_adj, redmod_adj)


fullmod2 <- lm(pheno~ as.factor(snp1)+as.factor(snp3)+as.factor(snp1):as.factor(snp3))
redmod2 <- lm(pheno~ as.factor(snp1)+as.factor(snp3))
intmod2 <- anova(fullmod2, redmod2)







