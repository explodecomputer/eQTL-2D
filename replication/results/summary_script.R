# Summary of the replication results
# Results sent by Harm-Jan Westra	


# Read in RData files

load("~/repo/eQTL-2D/replication/results/EGCUT_replication.RData")
newsig$id <- with(newsig, paste(probename, snp1, snp2))
EGCUT <- newsig

load("~/repo/eQTL-2D/replication/results/FehrmannHT12v3_replication.RData")
newsig$id <- with(newsig, paste(probename, snp1, snp2))
Feh <- subset(newsig, select=-c(chr1, chr2, pos1, pos2, snp1, snp2, complete, probeid, probename, probechr, probegene, probehsq, minclasssize, snpcor, propG, propA, pnest, pfull, pint))
rm(newsig)

rep <- merge(EGCUT, Feh, by="id", suff=c("_EGCUT", "_Feh"), all=T)
dim(rep)

# Pairwise correlations between 
cor(subset(rep, select=c(pnest, replication_pnest_EGCUT, replication_pnest_Feh)), use="pair")


panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use="pair"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)

    test <- cor.test(x,y)
    text(0.5, 0.5, txt, cex = cex)
}



temp <- subset(rep, select=c(pnest, replication_pnest_EGCUT, replication_pnest_Feh))
names(temp) <- c("BSGS", "EGCUT", "Fehrmann")
pdf("pnest_rep_pvals.pdf")
pairs(temp, lower.panel=panel.smooth, upper.panel=panel.cor, main="Interaction p-values in 3 datasets (4df)")
dev.off()

temp <- subset(rep, select=c(pfull, replication_pfull_EGCUT, replication_pfull_Feh))
names(temp) <- c("BSGS", "EGCUT", "Fehrmann")
pdf("pfull_rep_pvals.pdf")
pairs(temp, lower.panel=panel.smooth, upper.panel=panel.cor, main="Full genetic p-values in 3 datasets (8df)")
dev.off()


# QQ plots




# Overlap in replication results

match_rep_results.fun <- function(EGCUT, Feh) {


	out <- NULL

	for(i in 1:nrow(Feh)) {

		m1 <- as.character(Feh$probename[i])
		m2 <- as.character(Feh$snp1[i])
		m3 <- as.character(Feh$snp2[i])	

		index <- which(EGCUT$probename==m1 & EGCUT$snp1==m2 & EGCUT$snp2==m3)
		if(length(index)==1) {

			tmp <- EGCUT[index, c(22:29)]
			names(tmp) <- paste(names(tmp), "_EGCUT", sep="")
			tmp <- cbind(Feh[i,], tmp)
			out <- rbind(out, tmp)

		}

		if(length(index)!=1) {

			print(length(index))
			print(i)
			print("Look here!!")
		}
	}

	return(out)
}


out <- match_rep_results.fun(EGCUT, Feh)
write.csv(out, "/Volumes/group_wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/replication/results/matched_replication_results.csv", quote=F, row.names=F)

png("/Volumes/group_wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/docs/manuscript/images/EGCUT_vs_Fehrmann_all.png", width=600, height=600)
plot(out$replication_pnest, out$replication_pnest_EGCUT, xlab="-log10 pvalues from Fehrmann", ylab="-log10 pvalues from EGCUT")
dev.off()


index <- which(out$replication_pnest_EGCUT < 10)
png("/Volumes/group_wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/docs/manuscript/images/EGCUT_vs_Fehrmann_pval10.png", width=600, height=600)
plot(out$replication_pnest[index], out$replication_pnest_EGCUT[index], xlab="-log10 pvalues from Fehrmann", ylab="-log10 pvalues from EGCUT")
dev.off()











