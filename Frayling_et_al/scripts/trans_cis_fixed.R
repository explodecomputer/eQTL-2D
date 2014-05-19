# Trans mapping fixed for cis-snp
# Q: do we see an inflation in the p-vaues in a genome-wide analysis between 
# epistatsis snps when there are additive effects


# Run functions at the bottom of the script

# Read in data 
load("/Users/jpowell/repo/eQTL-2D/data/bsgs_egcut_fehr_data.RData")			### Path dir to corrected expression data
pheno <- phenlist[[1]]

# PLINK geno change function
load("/Users/jpowell/repo/eQTL-2D/data/geno.RData")




#=======================================================#
#		RUN ANALYSIS USING FUNCTIONS LISTED BELOW		#
#=======================================================#

# TMEM149 - ILMN_1786426
# CIS SNPS - rs8106959 

TMEM149_out <- trans_cis_epi.fun(geno, "rs8106959", pheno, "ILMN_1786426")
write.csv(TMEM149_out, "/Users/jpowell/repo/eQTL-2D/Frayling_et_al/data_files/TMEM149_out.csv", quote=F, row.names=F)

# MBLN1 - ILMN_2313158
# CIS SNPS - rs13069559 

MBLN1_out <- trans_cis_epi.fun(geno, "rs13069559", pheno, "ILMN_2313158")
write.csv(MBLN1_out, "/Users/jpowell/repo/eQTL-2D/Frayling_et_al/data_files/MBLN1_out.csv", quote=F, row.names=F)



#=======================================================#
#				ANALYSIS OF THE RESULTS					#
#=======================================================#

TMEM149_out <- cbind(bim[which(which(bim$V2=="rs8106959"))], TMEM149_out)
MBLN1_out <- cbind(bim[-which(bim$V2=="rs13069559"),], MBLN1_out)


# TMEM149 CHR
index <- which(TMEM149_out$V1==19 | TMEM149_out$V1==6 | TMEM149_out$V1==1 | TMEM149_out$V1==4 | TMEM149_out$V1==2 | TMEM149_out$V1==8 | TMEM149_out$V1==13)
pdf(file="/Users/jpowell/repo/eQTL-2D/Frayling_et_al/data_files/TMEM149_QQ.pdf")
qqplot.fun(10^-TMEM149_out$P[-index])
dev.off()

lambda.fun(10^-TMEM149_out$P)


# MBLN1 CHR
index <- which(MBLN1_out$V1==3 | MBLN1_out$V1==6 | MBLN1_out$V1==14 | MBLN1_out$V1==17 | MBLN1_out$V1==7)
pdf(file="/Users/jpowell/repo/eQTL-2D/Frayling_et_al/data_files/MBLN1_QQ.pdf")
qqplot.fun(10^-MBLN1_out$P[-index])
dev.off()

lambda.fun(10^-MBLN1_out$P)



#=======================================================#
#		*******			FUNCTIONS 		********		#
#=======================================================#


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



#=======================================================#
#		ANALYSIS OF EPI SCAN BY FIXED CIS SNP 			#
#=======================================================#

trans_cis_epi.fun <- function(geno, snp1_id, probe, probe_id) {


	snp1 <- geno[,which(colnames(geno)==snp1_id)]
	g <- geno[,-which(colnames(geno)==snp1_id)]				
	probe <- pheno[,which(colnames(pheno)==probe_id)]						

	nsnps <- ncol(g)	# ncol of geno data

	out <- array(0, c(nsnps, 5))				# blank array of results
	telliter <- 1000

	for(i in 1:nsnps) {
	
		snp2 <- g[,i]			# this is to vary by i (1:nsnps)
	

		fullmod <- lm(probe ~ as.factor(snp1) + as.factor(snp2) + as.factor(snp1):as.factor(snp2))
		redmod <- lm(probe ~ as.factor(snp1) + as.factor(snp2))
		intmod <- anova(fullmod, redmod)

		out[i,1] <- round(intmod$F[2],2)		# store F-statistics
		out[i,2] <- -log10(intmod$Pr[2])		# store P-values

		out[i,3] <- length(table(snp1 + 3*snp2))   	# nclass size
		out[i,4] <- min(table(snp1 + 3*snp2))		# min class size

		out[i,5] <- round(cor(snp1, snp2, use="pairwise.complete.obs"),2)^2		# LD / correlation between 2 snps 


		# print iteraction
		if(i %% telliter==0) {
			cat(paste("iteraction ", i, " complete\n"))
		}


	}

	out <- as.data.frame(out)
	names(out) <- c("F", "P", "nclass", "minclass", "LD")
	return(out)

}



#=======================================================#
#				MANHATTEN PLOT FUNCTION					#
#=======================================================#


Manhat_plot.fun <- function(data, x, name) {
	# data: Data
	# x: Position of the "red" line 
	# name: name of the figure to be written out

	#**# Define the colour system
	use=c("blue", "steelblue1")

	#**# Name for the output file
	name <- as.character(name)

	#==# order data by chromosome and bp
	# Data provided in correct order
	# Add cumulative bp
	data <- cbind(data, c(1:nrow(data)))
	
	#==# Determine cumulative bp per chromosome
	chr_ends_bp <- array(0, c(length(unique(data[,1]))))
	for(i in unique(data[,1])) {
		chr_ends_bp[i] <- max(data[data[,1]==i, 5])
	}

	for(i in 2:length(chr_ends_bp)) {
		chr_ends_bp[i] <- chr_ends_bp[i-1]+chr_ends_bp[i]
	}

	#==# Make cumulative bp position
	data_c <- cbind(data, data[,5])
	for(i in 2:length(chr_ends_bp)) {
		data_c[data_c[,4]==i, 6] <- data_c[data_c[,4]==i,5]+chr_ends_bp[i-1]		
	}

	#==# the new data file with cumulative position is called "data_c"
	#==# Make the graphic identifiers

	#==# Make colors for chromosomes
	col=use[data[,1]%%length(use)+1]

	#==# Make the x-axis position identifier 
	chr_mid <- array(0, c(length(unique(data[,1]))))

	out <- array(0, c(length(unique(data[,1]))))
	for(i in unique(data[,1])) {
		tmp <- length(which(data[,1]==unique(data[,1])[i]))
		out[i] <- tmp
	}

	for(i in unique(data[,1])) {
		tmp <- length(which(data[,1]==unique(data[,1])[i]))
		
		if(i==1) {
			chr_mid[i] <- tmp/2
		}
		else {
			chr_mid[i] <- sum(out[1:i-1])+(tmp/2) 
		}	
	}

	#==# Plot the "whole genome" plot
	jpeg(paste(name, ".jpg", sep=""), height=480, width=2200)
	plot(-log10(data[,4]), ylim=c(0,7), col=col, xlab="Chromosome", ylab="-log10 p values", axes=F, pch=19, cex=0.7)
	axis(1, at=chr_mid, labels=c(1:22), las=2, cex.axis=1, las=1)
	axis(2, at=seq(0,max(-log10(data[,4]))+1,1), labels=seq(0,max(-log10(data[,4]))+1,1), las=1, cex.axis=1)
	abline(x, 0, col="red")
	dev.off()
}


#=======================================================#
#					QQ-PLOT FUNCTION					#
#=======================================================#


qqplot.fun <- function(pvector, main=NULL, ...) {
    o = -log10(sort(pvector,decreasing=F))
    e = -log10( 1:length(o)/length(o) )
    plot(e,o,pch=19,cex=1, main=main, ...,
        xlab=expression(Expected~~-log[10](italic(p))),
        ylab=expression(Observed~~-log[10](italic(p))),
        xlim=c(0,max(e)), ylim=c(0,max(o)))
    lines(e,e,col="red")
}
 

#=======================================================#
#			ESTIMATING LAMBDA FUNCTION					#
#=======================================================#



lambda.fun <- function (data, plot = FALSE, proportion = 1, method = "regression", 
    filter = TRUE, df = 1, ...) 
	{
	    data <- data[which(!is.na(data))]
	    if (proportion > 1 || proportion <= 0) 
	        stop("proportion argument should be greater then zero and less than or equal to one")
	    ntp <- round(proportion * length(data))
	    if (ntp < 1) 
	        stop("no valid measurements")
	    if (ntp == 1) {
	        warning(paste("One measurement, lambda = 1 returned"))
	        return(list(estimate = 1, se = 999.99))
	    }
	    if (ntp < 10) 
	        warning(paste("number of points is too small:", ntp))
	    if (min(data) < 0) 
	        stop("data argument has values <0")
	    if (max(data) <= 1) {
	        data <- qchisq(data, 1, lower.tail = FALSE)
	    }
	    if (filter) {
	        data[which(abs(data) < 1e-08)] <- NA
	    }
	    data <- sort(data)
	    ppoi <- ppoints(data)
	    ppoi <- sort(qchisq(ppoi, df = df, lower.tail = FALSE))
	    data <- data[1:ntp]
	    ppoi <- ppoi[1:ntp]
	    out <- list()
	    if (method == "regression") {
	        s <- summary(lm(data ~ 0 + ppoi))$coeff
	        out$estimate <- s[1, 1]
	        out$se <- s[1, 2]
	    }
	    else if (method == "median") {
	        out$estimate <- median(data, na.rm = TRUE)/qchisq(0.5, 
	            df)
	        out$se <- NA
	    }
	    else if (method == "KS") {
	        limits <- c(0.5, 100)
	        out$estimate <- estLambdaKS(data, limits = limits, df = df)
	        if (abs(out$estimate - limits[1]) < 1e-04 || abs(out$estimate - 
	            limits[2]) < 1e-04) 
	            warning("using method='KS' lambda too close to limits, use other method")
	        out$se <- NA
	    }
	    else {
	        stop("'method' should be either 'regression' or 'median'!")
	    }
	    if (plot) {
	        lim <- c(0, max(data, ppoi, na.rm = TRUE))
	        oldmargins <- par()$mar
	        par(mar = oldmargins + 0.2)
	        plot(ppoi, data, xlab = expression("Expected " ~ chi^2), 
	            ylab = expression("Observed " ~ chi^2), ...)
	        abline(a = 0, b = 1)
	        abline(a = 0, b = out$estimate, col = "red")
	        par(mar = oldmargins)
	    }
    return(out)
}



