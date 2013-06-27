summariseRep <- function(filename)
{
	load(filename)
	newsig$gcm <- gcm
	newsig$gcs <- gcs

	v <- matrix(0, nrow(newsig), 8)
	for(i in 1:nrow(newsig))
	{
		v[i,] <- mod[[i]]$variances[-1]
		vp <- var(mod[[i]]$phen, na.rm=TRUE)
		v[i,] <- v[i,] / vp
	}
	v <- as.data.frame(v)
	nom <- names(mod[[1]]$variances[-1])
	names(v) <- nom

	sig <- cbind(newsig, v)

	return(sig)
}

fehr <- summariseRep("~/repo/eQTL-2D/replication/results/replication_GrngHT12v3.RData")
egcut <- summariseRep("~/repo/eQTL-2D/replication/results/replication_EGCUT.RData")

save(fehr, egcut, file="~/repo/eQTL-2D/replication/results/replication2_summarised.RData")



load("~/repo/eQTL-2D/replication/results/replication_GrngHT12v3.RData")



