nestedTestVars <- function(mod)
{
	phen <- mod$phen
	x <- mod$gen
	mod1 <- lm(phen ~ as.factor(x[,1]) * as.factor(x[,2]))
	mod2 <- lm(phen ~ as.factor(x[,1]) + as.factor(x[,2]))
	nest <- anova(mod2, mod1)
	return(nest)
}


metaNestedTest <- function(nest1, nest2)
{
	nT <- nest1$Res[2] + 9 + nest2$Res[2] + 9
	RSST <- nest1$RSS + nest2$RSS
	Fval <- ((RSST[1] - RSST[2]) / 4) / (RSST[2] / (nT - 9))
	p <- pf(Fval, 4, nT - 4, lower.tail=FALSE)
	return(-log10(p))
}


metaTestChisq <- function(nest1, nest2)
{
	F1 <- nest1$F[2]
	F2 <- nest2$F[2]
	return(-log10(pchisq(F1 + F2, 8, lower.tail=TRUE)))
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


combineReplicationData <- function(mod1, mod2)
{
	x <- rbind(mod1$gen, mod2$gen)
	y <- c(mod1$phen, mod2$phen)
	return(list(phen = y, gen = x))
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


qqPlot <- function(dat)
{
	p <- qqplot2 <- ggplot(dat) +
		geom_ribbon(aes(x=expect, ymin=lower, ymax=upper), colour="white", alpha=0.5) +
		geom_abline(intercept=0, slope=1) +
		geom_point(aes(x=expect, y=observed, colour=ex), size=1.5) +
		facet_grid(. ~ fake) +
		labs(colour="Above FDR 5% CI?", y="Observed", x="Expected") +
		scale_colour_brewer(type="qual", palette=3) +
		theme(legend.position="none") + ylim(c(0, 7))
	return(p)
}


#=============================================================#
#=============================================================#


load("~/repo/eQTL-2D/analysis/interaction_list_replication_summary.RData")

load("~/repo/eQTL-2D/replication/results/replication_GrngHT12v3.RData")
newsig$code <- with(newsig, paste(probename, snp1, snp2))
fehr <- mod
names(fehr) <- newsig$code

load("~/repo/eQTL-2D/replication/results/replication_EGCUT.RData")
newsig$code <- with(newsig, paste(probename, snp1, snp2))
egcut <- mod
names(egcut) <- newsig$code



#=============================================================#
#=============================================================#


sig_all$pnest_meta <- NA
sig_all$pnest_meta2 <- NA
sig_all$code <- with(sig_all, paste(probename, snp1, snp2))

for(i in 1:nrow(sig_all))
{
	cat(i, "\n")
	if(is.na(sig_all$pnest_egcut[i]) | is.na(sig_all$pnest_fehr[i]) | length(sig_all$gcs_egcut[[i]]) != 9 | length(sig_all$gcs_fehr[[i]]) != 9)
	{
		cat("skipped\n")
		next
	}
	b1 <- nestedTestVars(egcut[[sig_all$code[i]]])
	b2 <- nestedTestVars(fehr[[sig_all$code[i]]])
	sig_all$pnest_meta[i] <- metaNestedTest(b1, b2)
	sig_all$pnest_meta2[i] <- metaTestChisq(b1, b2)
}

thresh <- -log10(0.05 / 473)
with(sig_all, table(filter, pnest_meta > thresh))
with(sig_all, table(filter, pnest_meta2 > thresh))
with(sig_all, table(filter, is.na(pnest_meta)))

meta <- subset(sig_all, !is.na(pnest_meta))
meta <- qqDat(meta, 0.05, "pnest_meta")
meta$ex <- meta$pnest_meta > meta$upper
meta$fake <- meta$filter == 3
meta$observed <- meta$pnest_meta
qqPlot(meta)



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


alleleFreq(sig$gcs_egcut[[1]])

af_bsgs <- do.call(rbind, lapply(sig_all$gcs, alleleFreq))
af_egcut <- do.call(rbind, lapply(sig_all$gcs_egcut, alleleFreq))
af_fehr <- do.call(rbind, lapply(sig_all$gcs_fehr, alleleFreq))

pairs(cbind(af_bsgs[,1], af_egcut[,1], af_fehr[,1]))

sig_all$gcmRep <- NA
sig_all$gcsRep <- NA



