#=======================================================#
#-------------------------------------------------------#
#														#
#	vc_effects.R	 									#
#														#
#	estimate the additive and non-additive vc before 	#
#	and after fitting the snps given in the wood et al	#
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



#=======================================================#
#				
#=======================================================#

i <- 1
snp1  <- block[,which(colnames(block)==as.character(info$SNP1[i]))] 
snp2  <- block[,which(colnames(block)==as.character(info$SNP2[i]))] 
inc_snp <- block[,which(colnames(block)==as.character(info$rs_id[i]))] 

pheno <- bsgs[,which(colnames(bsgs)==as.character(info$Probe[i]))]
pheno_adj <- summary(lm(pheno ~ inc_snp))$residuals



VCTests.fun <- function(snp1, snp2, inc_snp, pheno)
{
	require(noia)
	# Extract data

	effects_all <- NULL
	effects <- array(0, c(1,36))
	se_all <- NULL
	se <- array(0, c(1,36))


	for(i in 1:nrow(info)) {
	
		if(i==6) {
			effects[1:36] <- NA

		}

		else {
			snp1  <- block[,which(colnames(block)==as.character(info$SNP1[i]))] 
			snp2  <- block[,which(colnames(block)==as.character(info$SNP2[i]))] 
			inc_snp <- block[,which(colnames(block)==as.character(info$rs_id[i]))] 

			pheno <- bsgs[,which(colnames(bsgs)==as.character(info$Probe[i]))]
			pheno_adj <- summary(lm(pheno ~ inc_snp))$residuals


			# Summary statistics original
			tab1 <- table(snp1 + 3*snp2)
			gcm1 <- tapply(pheno, list(snp1, snp2), function(x) { mean(x, na.rm=T)})
			mod1 <- linearRegression(pheno, cbind(snp1, snp2)+1)

			# Summary statistics adjusted
			tab2 <- table(snp1 + 3*snp2)
			gcm2 <- tapply(pheno, list(snp1, snp2), function(x) { mean(x, na.rm=T)})
			mod2 <- linearRegression(pheno_adj, cbind(snp1, snp2)+1)

			# Summary statistics inc snp1
			tab3 <- table(snp1 + 3*inc_snp)
			gcm3 <- tapply(pheno, list(snp1, inc_snp), function(x) { mean(x, na.rm=T)})
			mod3 <- linearRegression(pheno, cbind(snp1, inc_snp)+1)

			# Summary statistics inc snp1
			tab4 <- table(inc_snp + 3*snp2)
			gcm4 <- tapply(pheno, list(inc_snp, snp2), function(x) { mean(x, na.rm=T)})
			mod4 <- linearRegression(pheno, cbind(inc_snp, snp2)+1)

			effects[1:9] <- mod1$E
			effects[10:18] <- mod2$E
			effects[19:27] <- mod3$E
			effects[28:36] <- mod4$E

			se[1:9] <- mod1$std.err
			se[10:18] <- mod2$std.err
			se[19:27] <- mod3$std.err
			se[28:36] <- mod4$std.err


		}	

		effects_all <- rbind(effects_all, effects)
		se_all <- rbind(se_all, se)
	}	
	
	effects_all <- round(effects_all, 3)	
	effects_all <- as.data.frame(effects_all)
	names(effects_all)[1:9] <- paste(names(mod1$E), "_ori", sep="")
	names(effects_all)[10:18] <- paste(names(mod1$E), "_adj", sep="")
	names(effects_all)[19:27] <- paste(names(mod1$E), "_s1inc", sep="")
	names(effects_all)[28:36] <- paste(names(mod1$E), "_s2inc", sep="")

	se_all <- round(se_all, 3)	
	se_all <- as.data.frame(se_all)
	names(se_all)[1:9] <- paste(names(mod1$E), "_ori", sep="")
	names(se_all)[10:18] <- paste(names(mod1$E), "_adj", sep="")
	names(se_all)[19:27] <- paste(names(mod1$E), "_s1inc", sep="")
	names(se_all)[28:36] <- paste(names(mod1$E), "_s2inc", sep="")


	return(l)
}


#' Run replication analysis
#'
#' Test all SNP pairs in interaction list in replication dataset.
#'
#' @param sig Output from \link{LoadIntList}
#' @param checked Output from \link{DataChecks}
#'
#' @return Returns \code{data.frame} with new columns for results from replication data
#' @export
RunReplication <- function(sig, checked)
{
	sig$replication_pfull <- NA
	sig$replication_pnest <- NA
	sig$replication_p1 <- NA
	sig$replication_p2 <- NA
	sig$replication_r <- NA
	sig$replication_nclass <- NA
	sig$replication_minclass <- NA
	sig$replication_nid <- NA

	geno <- checked$geno
	probes <- checked$probes

	gcm <- list()
	gcs <- list()
	mod <- list()

	for(i in 1:nrow(sig))
	{
		cat(i, "of", nrow(sig), "\n")
		out <- ReplicationTests(geno, probes, sig, i)
		sig <- out$sig
		gcm[[i]] <- out$gcm
		gcs[[i]] <- out$gcs
		mod[[i]] <- out$mod
	}

	l <- list()
	l$sig <- sig
	l$gcm <- gcm
	l$gcs <- gcs
	l$mod <- mod

	return(l)
}



plot(pheno, pheno_adj)
