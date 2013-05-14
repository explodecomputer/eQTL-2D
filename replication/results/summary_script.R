# Summary of the replication results
# Results sent by Harm-Jan Westra	


# Read in RData files

load("/Volumes/group_wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/replication/results/EGCUT_replication.RData")
EGCUT <- newsig
load("/Volumes/group_wrayvisscher/josephP/Projects/Genomics/Expression_main/Epistasis/eQTL-2D/replication/results/FehrmannHT12v3_replication.RData")
Feh <- newsig
rm(newsig)


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











