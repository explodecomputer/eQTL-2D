F(4,1000)= (chi(4)/4) / (chi(1000)/1000)) ~ chi(4)/4

This is assuming that your test with 4 df has a large number of df in the
denominator.

So 4*F ~ chi(4)

This implies that T = 4*F(cohort1) + 4*F(cohort2) is distributed as chi(8)
under the null of no epistasis, and you can get a p-value from a chisquare
with 8 df. You could also compare T for the null SNPs and see if it
deviates from its expectation of 8.0 under the null.

You could also do a MA on just using the p-values from the 2 replication
studies. Under the null, -2*ln(p) follows a chi(2) distribution, so T =
-2ln(p from cohort1) - 2ln(p from cohort2) follows a chi(4) distribution.


library(ggplot2)

nestedTestVars <- function(mod)
{
	phen <- mod$phen
	x <- mod$gen
	mod1 <- lm(phen ~ as.factor(x[,1]) * as.factor(x[,2]))
	mod2 <- lm(phen ~ as.factor(x[,1]) + as.factor(x[,2]))
	nest <- anova(mod2, mod1)
	return(nest)
}


compareLargestVc <- function(s, mod1, mod2)
{
	vcepi <- cbind(c(4,5,7,8), s$vc[[1]][c(4, 5, 7, 8)])

	i <- vcepi[which.max(vcepi[,2]), 1] + 1
	p1 <- mod1[[s$code]]$pvalues[i]
	p2 <- mod2[[s$code]]$pvalues[i]
	return(-log10(pchisq(-2*log(p1) - 2*log(p2), 4, lower.tail=FALSE)))
}	


metaTestChisq <- function(nest1, nest2)
{
	F1 <- nest1$F[2]
	F2 <- nest2$F[2]
	return(-log10(pchisq(4*F1 + 4*F2, 8, lower.tail=FALSE)))
}


metaTestPval <- function(nest1, nest2)
{
	p1 <- nest1$P[2]
	p2 <- nest2$P[2]
	return(-log10(pchisq(-2*log(p1) - 2*log(p2), 4, lower.tail=FALSE)))
}


performMeta <- function(sig, mod1, mod2)
{
	sig$pnest_meta <- NA
	sig$pnest_meta2 <- NA
	sig$pnest_meta_vc <- NA
	sig$code <- with(sig, paste(probename, snp1, snp2))

	for(i in 1:nrow(sig))
	{
		cat(i, "\n")
		if(is.na(sig$pnest_egcut[i]) | is.na(sig$pnest_fehr[i]) | length(sig$gcs_egcut[[i]]) != 9 | length(sig$gcs_fehr[[i]]) != 9)
		{
			cat("skipped\n")
			next
		}
		b1 <- nestedTestVars(mod1[[sig$code[i]]])
		b2 <- nestedTestVars(mod2[[sig$code[i]]])
		sig$pnest_meta[i] <- metaTestChisq(b1, b2)
		sig$pnest_meta2[i] <- metaTestPval(b1, b2)
		sig$pnest_meta_vc[i] <- compareLargestVc(sig[i, ], mod1, mod2)
	}
	return(sig)
}


confInt <- function(n, alpha)
{
	k <- c(1:n)
	upper <- -log10(qbeta(alpha/2,k,n+1-k))
	lower <- -log10(qbeta((1-alpha/2),k,n+1-k))
	expect <- -log10((k-0.5)/n)
	return(data.frame(expect, lower, upper))
}


qqDat <- function(sig, alpha, col)
{
	a <- subset(sig, filter==3)
	a <- a[order(a[[col]], decreasing=T), ]
	ci <- confInt(nrow(a), alpha)
	a <- data.frame(a, ci)

	b <- subset(sig, filter!=3)
	b <- b[order(b[[col]], decreasing=T), ]
	ci <- confInt(nrow(b), alpha)
	b <- data.frame(b, ci)

	return(rbind(b, a))
}

makeQqDat <- function(sig, alpha, col)
{
	meta <- subset(sig, !is.na(pnest_meta))
	meta <- qqDat(meta, 0.05, "pnest_meta")
	meta$ex <- meta$pnest_meta > meta$upper
	meta$fake <- "Null"
	meta$fake[meta$filter != 3] <- "Empirical"
	meta$observed <- meta$pnest_meta
	return(meta)
}


qqPlot <- function(dat, lim)
{
	p <- qqplot2 <- ggplot(dat) +
		geom_ribbon(aes(x=expect, ymin=lower, ymax=upper), colour="white", alpha=0.5) +
		geom_abline(intercept=0, slope=1) +
		geom_point(aes(x=expect, y=observed, colour=ex), size=1.5) +
		facet_grid(. ~ fake) +
		labs(colour="Above FDR 5% CI?", y="Observed", x="Expected") +
		scale_colour_brewer(type="qual", palette=3) +
		theme(legend.position="none") + ylim(c(0, lim))
	return(p)
}

qqPlot2 <- function(dat, thresh)
{
	dat <- subset(dat, filter != 3)
	dat1 <- subset(dat, pnest_meta < thresh)
	dat$code <- "All"
	dat1$code <- "Below Bonferroni threshold"
	dat <- rbind(dat, dat1)
	p <- qqplot2 <- ggplot(dat) +
		geom_ribbon(aes(x=expect, ymin=lower, ymax=upper), colour="white", alpha=0.5) +
		geom_abline(intercept=0, slope=1) +
		geom_point(aes(x=expect, y=observed, colour=ex), size=1.5) +
		facet_grid(code ~ ., scales = "free_y") +
		geom_hline(yintercept=thresh, linetype="dotted") +
		labs(colour="Above FDR 5% CI?", y="Observed", x="Expected") +
		scale_colour_brewer(type="qual", palette=3) +
		theme(legend.position="none")
	return(p)
}

qqPlot3 <- function(dat, thresh)
{
	dat <- subset(dat, filter != 3)
	dat1 <- subset(dat, pnest_meta < thresh)
	dat$code <- "All"
	dat1$code <- "Below Bonferroni threshold"
	dat <- rbind(dat, dat1)
	p <- qqplot2 <- ggplot(dat) +
		geom_ribbon(aes(x=expect, ymin=lower, ymax=upper), colour="white", alpha=0.5) +
		geom_abline(intercept=0, slope=1) +
		geom_point(aes(x=expect, y=observed, colour=ex), size=1.5) +
		facet_grid(code ~ ., scales = "free_y") +
		geom_hline(yintercept=thresh, linetype="dotted") +
		labs(colour="Above FDR 5% CI?", y="Observed", x="Expected") +
		scale_colour_manual(values=c("#556270", "#556270")) +
		theme(legend.position="none")
	return(p)
}


alleleFreq <- function(x)
{
	if(is.na(x[1]))
	{
		return(c(NA, NA))
	}
	a1 <- rev(colSums(x))
	a2 <- rev(rowSums(x))

	a1 <- 1 - (a1[1] + a1[2]/2)/sum(a1)
	a2 <- 1 - (a2[1] + a2[2]/2)/sum(a2)

	return(c(a1, a2))
}


combineReplicationData <- function(mod1, mod2)
{
	x <- rbind(mod1$gen, mod2$gen)
	y <- c(mod1$phen, mod2$phen)
	return(list(phen = y, gen = x))
}


createAverageGp <- function(gcm, gcs)
{
	if(dim(gcm[[1]]) != c(3,3))
	{
		return(NA)
	}

	nrep <- length(gcm)

	gcsT <- matrix(0, 3,3)
	gcmT <- matrix(0, 3,3)

	for(i in 1:nrep)
	{
		gcsT[1,1] <- gcsT[1,1] + gcm[[i]][1,1]
		gcsT[1,2] <- gcsT[1,2] + gcm[[i]][1,2]
		gcsT[1,3] <- gcsT[1,3] + gcm[[i]][1,3]
		gcsT[2,1] <- gcsT[2,1] + gcm[[i]][2,1]
		gcsT[2,2] <- gcsT[2,2] + gcm[[i]][2,2]
		gcsT[2,3] <- gcsT[2,3] + gcm[[i]][2,3]
		gcsT[3,1] <- gcsT[3,1] + gcm[[i]][3,1]
		gcsT[3,2] <- gcsT[3,2] + gcm[[i]][3,2]
		gcsT[3,3] <- gcsT[3,3] + gcm[[i]][3,3]
	}

	for(i in 1:nrep)
	{
		gcs[[i]] <- gcs[[i]] / gcsT
		
		gcmT[1,1] <- gcmT[1,1] + gcm[[i]][1,1] * gcs[[i]][1,1]
		gcmT[1,2] <- gcmT[1,2] + gcm[[i]][1,2] * gcs[[i]][1,2]
		gcmT[1,3] <- gcmT[1,3] + gcm[[i]][1,3] * gcs[[i]][1,3]
		gcmT[2,1] <- gcmT[2,1] + gcm[[i]][2,1] * gcs[[i]][2,1]
		gcmT[2,2] <- gcmT[2,2] + gcm[[i]][2,2] * gcs[[i]][2,2]
		gcmT[2,3] <- gcmT[2,3] + gcm[[i]][2,3] * gcs[[i]][2,3]
		gcmT[3,1] <- gcmT[3,1] + gcm[[i]][3,1] * gcs[[i]][3,1]
		gcmT[3,2] <- gcmT[3,2] + gcm[[i]][3,2] * gcs[[i]][3,2]
		gcmT[3,3] <- gcmT[3,3] + gcm[[i]][3,3] * gcs[[i]][3,3]
	}

	return(list(gcm = gcmT, gcs = gcsT))
}



load(file="~/repo/eQTL-2D/analysis/interaction_list_meta_analysis.RData")
load(file="~/repo/eQTL-2D/data/replication_null.RData")

meta2 <- subset(meta, filter != 3)
meta3 <- subset(meta, filter != 3)
meta3$filter <- 3
meta3$fake <- "Null"
meta3$pnest_meta <- meta3$observed <- sort(-log10(pnull$meta), decreasing=T)
meta3$ex <- meta3$observed > meta3$upper
meta1 <- rbind(meta2, meta3)


thresh <- -log10(0.05 / 434)
with(meta1, table(filter, pnest_meta > thresh))
with(meta1, table(filter, pnest_meta2 > thresh))
with(meta1, table(filter, pnest_meta_vc > thresh))
with(meta1, table(filter, is.na(pnest_meta)))


#=============================================================#
#=============================================================#

# Make Q-Q plots for meta analysis


qqPlot(meta1, thresh)
ggsave(file="~/repo/eQTL-2D/analysis/images/qqMetaNonsig.pdf", width=10, height=5)
qqPlot(meta, max(meta$pnest_meta))
ggsave(file="~/repo/eQTL-2D/analysis/images/qqMetaAll.pdf", width=5, height=10)
qqPlot2(meta, thresh)
ggsave(file="~/repo/eQTL-2D/analysis/images/qqMeta.pdf", width=5, height=10)

with(meta, table(filter == 3, pnest_meta > upper))
with(meta, table(filter == 3, pnest_fehr > upper_fehr))
with(meta, table(filter == 3, pnest_egcut > upper_egcut))

# Make Q-Q plot for poster (Yellow)
qqPlot3(meta, thresh)
ggsave(file="~/repo/eQTL-2D/analysis/images/qqMeta_yellow.pdf", width=5, height=10)




#=============================================================#
#=============================================================#

# Calculate lambda values for meta analysis


lambdaVal <- function(pval)
{
	qchisq(1-mean(10^-pval, na.rm=T), 1) / qchisq(0.5, 1)
}


lambdaVal2 <- function(pval)
{
	mean(qchisq(10^-pval, 1, lower.tail=FALSE))
}


p1 <- subset(meta, !is.na(pnest_meta) & filter != 3)$pnest_meta
p2 <- subset(meta, !is.na(pnest_meta) & filter == 3)$pnest_meta
p3 <- subset(meta, !is.na(pnest_meta) & filter == 3)$expect
p4 <- subset(meta, pnest_meta < thresh & filter != 3)$pnest_meta
mean(qchisq(p2/2.3, 2), na.rm=T)

p1 <- p1[order(p1, decreasing=T)]
p1a <- p1[-c(1:20)]

lambdaVal2(p1)
lambdaVal2(p2)
lambdaVal2(p3)
lambdaVal2(p4)

lambdaVal2(sig$pnest_fehr)
lambdaVal2(sig$pnest_egcut)
lambdaVal2(subset(sig_all, filter == 3)$pnest_egcut)
lambdaVal2(subset(sig_all, filter == 3)$pnest_fehr)

lambdaVal2(subset(sig_all, filter != 3)$pnest_egcut)
lambdaVal2(subset(sig_all, filter != 3)$pnest_fehr)

lambdaVal2(subset(meta, filter != 3)$pnest_meta)
mean(c(lambdaVal2(subset(sig_all, filter != 3)$pnest_egcut), lambdaVal2(subset(sig_all, filter != 3)$pnest_fehr)))


#=============================================================#
#=============================================================#

# 

af_bsgs <- 1 - do.call(rbind, lapply(sig_all$gcs, alleleFreq))
af_egcut <- do.call(rbind, lapply(sig_all$gcs_egcut, alleleFreq))
af_fehr <- do.call(rbind, lapply(sig_all$gcs_fehr, alleleFreq))

pdf("~/repo/eQTL-2D/analysis/images/alleleFreq.pdf")
pairs(data.frame(BSGS = af_bsgs[,1], EGCUT = af_egcut[,1], Fehrmann = af_fehr[,1]))
dev.off()


pairs(data.frame(BSGS = c(af_bsgs), EGCUT = c(af_egcut), Fehrmann = c(af_fehr)))


x_fehr <- cbind(abs(af_bsgs[,1] - af_fehr[,1]), abs(af_bsgs[,2] - af_fehr[,2]))
y_fehr <- abs(sig_all$pnest - sig_all$pnest_fehr)
anova(lm(y_fehr ~ x_fehr))



#=============================================================#
#=============================================================#



sig_all$pnest_combined <- NA
for(i in 1:nrow(sig_all))
{
	cat(i, "\n")
	if(is.na(sig_all$pnest_egcut[i]) | is.na(sig_all$pnest_fehr[i]) | length(sig_all$gcs_egcut[[i]]) != 9 | length(sig_all$gcs_fehr[[i]]) != 9)
	{
		cat("skipped\n")
		next
	}
	m <- combineReplicationData(egcut[[sig_all$code[i]]], fehr[[sig_all$code[i]]])
	sig_all$pnest_combined[i] <- -log10(nestedTestVars(m)$P[2])
}

thresh <- -log10(0.05 / 473)
with(sig_all, table(filter, pnest_combined > thresh))
with(sig_all, table(filter, is.na(pnest_combined)))

combined <- subset(sig_all, !is.na(pnest_combined))
combined <- qqDat(combined, 0.05, "pnest_combined")
combined$ex <- combined$pnest_combined > combined$upper
combined$fake <- combined$filter == 3
combined$observed <- combined$pnest_combined
qqPlot(combined)

with(combined, table(filter, ex))





sig_all$pnest_combined <- NA
for(i in 1:nrow(sig_all))
{
	cat(i, "\n")
	if(is.na(sig_all$pnest_egcut[i]) | is.na(sig_all$pnest_fehr[i]) | length(sig_all$gcs_egcut[[i]]) != 9 | length(sig_all$gcs_fehr[[i]]) != 9)
	{
		cat("skipped\n")
		next
	}
	m <- combineReplicationData(egcut[[sig_all$code[i]]], fehr[[sig_all$code[i]]])
	sig_all$pnest_combined[i] <- -log10(nestedTestVars(m)$P[2])
}

thresh <- -log10(0.05 / 473)
with(sig_all, table(filter, pnest_combined > thresh))
with(sig_all, table(filter, is.na(pnest_combined)))

combined <- subset(sig_all, !is.na(pnest_combined))
combined <- qqDat(combined, 0.05, "pnest_combined")
combined$ex <- combined$pnest_combined > combined$upper
combined$fake <- combined$filter == 3
combined$observed <- combined$pnest_combined
qqPlot(combined)

with(combined, table(filter, ex))



#=============================================================#
#=============================================================#


