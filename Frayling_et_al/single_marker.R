
# single inc snp analysis

load("repo/eQTL-2D/data/bsgs_egcut_fehr_data.RData"

pheno <- phenlist[[1]]
map <- read.table("repo/eQTL-2D/Frayling_et_al/test/snp_list.map", header=F)
ped <- read.table("repo/eQTL-2D/Frayling_et_al/test/snp_list.ped", header=F)
info <- read.csv("repo/eQTL-2D/Frayling_et_al/data_files/inc_info_bsgs.csv", header=T)



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


R2 <- array(0, c(nrow(info)))
for(i in 1:nrow(info)) {
	a1 <- which(colnames(pheno)==as.character(info$Probe[i]))
	a2 <- which(colnames(block)==as.character(info$rs_id[i]))

	if(length(a2)==0){
		R2[i] <- NA
	}

	if(length(a2)!=0){
		pheno1 <- pheno[,a1]
		snp1 <- block[,a2]
		R2[i] <- summary(lm(pheno1~snp1))$r.squared
	}

}

R2 <- round(R2, 3)

out <- cbind(info, R2)

head(in)

