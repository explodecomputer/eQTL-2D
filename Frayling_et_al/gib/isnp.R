library(plyr)
library(lattice)
library(latticeExtra)
library(gridExtra)


load("~/repo/eQTL-2D/data/bsgs_egcut_fehr_data.RData")
bsgsphen <- phenlist[[1]]
bsgsgeno <- genolist[[1]]

bim <- read.table("isnpall.bim")
names(bim) <- c("chr", "snp", "gd", "pd", "a1", "a2")

tab <- read.table("isnptable.txt", he=T)
tab$ipos <- as.numeric(do.call(rbind, strsplit(as.character(tab$snpi), split=":"))[,2])
tab <- merge(tab, bim, by.x="ipos", by.y="pd")


readXmat <- function(filename)
{
	a <- read.table(filename)
	snames <- as.character(unlist(a[1,-c(1:2)]))
	ids <- a[-c(1:2), 1:2]
	isnp <- matrix(as.numeric(as.character(unlist(a[-c(1:2), -c(1:2)]))), nrow(a)-2)
	colnames(isnp) <- snames
	return(isnp)
}

isnp <- readXmat("isnpall.xmat.gz")



l <- list()
m1 <- list()
m2 <- list()
for(i in 1:nrow(tab))
{
	cat(i, "\n")
	phen <- bsgsphen[,which(colnames(bsgsphen) == tab$probe[i])[1]]
	snp1 <- bsgsgeno[,which(colnames(bsgsgeno) == tab$snp1[i])[1]]
	snp2 <- bsgsgeno[,which(colnames(bsgsgeno) == tab$snp2[i])[1]]
	snpi <- isnp[,which(colnames(isnp) == tab$snp[i])[1]]
	# snpi[1:100] <- sample(snpi, 100, rep=T)
	res <- residuals(lm(phen ~ snpi))
	l[[i]] <- data.frame(phen = phen, snp1 = snp1, snp2 = snp2, snpi = snpi, res = res)
	m1[[i]] <- tapply(phen, list(snp1, snp2), function(x) mean(x, na.rm=T))
	m2[[i]] <- tapply(res, list(snp1, snp2), function(x) mean(x, na.rm=T))
}


pvals <- ldply(l, function(x)
{
	mod1 <- lm(phen ~ as.factor(snp1) * as.factor(snp2), x)
	mod2 <- lm(phen ~ as.factor(snp1) + as.factor(snp2), x)
	p1 <- anova(mod1, mod2)$P[2]

	mod1 <- lm(res ~ as.factor(snp1) * as.factor(snp2), x)
	mod2 <- lm(res ~ as.factor(snp1) + as.factor(snp2), x)
	p2 <- anova(mod1, mod2)$P[2]

	return(data.frame(p1 = p1, p2 = p2))
})

tab$p1 <- pvals$p1
tab$p2 <- pvals$p2
tab

subset(tab, cistrans == "cis")
subset(tab, cistrans == "trans")
subset(tab, p2 <= 0.05)

table(tab$p2 < 0.05, tab$cistrans)


plot3dGp <- function(gp, title="", snp1="SNP1", snp2="SNP2", z=-45)
{
	p <- cloud(
		gp, # This has to be a table!?
		panel.3d.cloud=panel.3dbars, 
		col="black", 
		col.facet=c("#e5f5e0", "#A1D99B", "#31A354"), 
		# col.facet=c("#31A354"), 
		xbase=0.4, 
		ybase=0.4,
		xlab=list(
			label=snp1, 
			cex=0.8), 
		ylab=list(
			label=snp2, 
			cex=0.8), 
		scales=list(
			arrows=F, 
			z = list(cex = 0), 
			y = list(cex = 0.5), 
			x=list(cex=0.5)),
		zlab="",
		cex.title = 0.5,
		screen = list(
			z = z, 
			x = -60, 
			y = 3),
		main = list(
			label=title, 
			cex=0.8)
	)
	return(p)
}


save(m1, m2, tab, l, file="isnp_analysis.RData")

m1 <- lapply(m1, function(x) x - min(x))
m2 <- lapply(m2, function(x) x - min(x))

ga <- function(i, m1=m1, m2=m2)
{
	res <- grid.arrange(plot3dGp(m1[[i]], title=paste("Original", i), plot3dGp(m2[[i]], title=paste("Residuals", i)))
	return(res)
}

i <- 21
res <- grid.arrange(plot3dGp(m1[[i]], title=paste("Original", i), plot3dGp(m2[[i]], title=paste("Residuals", i))))


pdf(file="isnp_gpmaps.pdf")
i <- 4
res <- grid.arrange(plot3dGp(m1[[i]], title=paste("Original", i), z=-35), plot3dGp(m2[[i]], title=paste("Residuals", i), z=-35))
i <- 6
res <- grid.arrange(plot3dGp(m1[[i]], title=paste("Original", i), z=-35), plot3dGp(m2[[i]], title=paste("Residuals", i), z=-35))
i <- 7
res <- grid.arrange(plot3dGp(m1[[i]], title=paste("Original", i), z=-35), plot3dGp(m2[[i]], title=paste("Residuals", i), z=-35))
i <- 8
res <- grid.arrange(plot3dGp(m1[[i]], title=paste("Original", i), z=-35), plot3dGp(m2[[i]], title=paste("Residuals", i), z=-35))
i <- 9
res <- grid.arrange(plot3dGp(m1[[i]], title=paste("Original", i), z=35), plot3dGp(m2[[i]], title=paste("Residuals", i), z=35))
i <- 22
res <- grid.arrange(plot3dGp(m1[[i]], title=paste("Original", i), z=35), plot3dGp(m2[[i]], title=paste("Residuals", i), z=35))
dev.off()


l20 <- l[[20]]

l20$res2 <- residuals(lm(phen ~ snp2, l20))
head(l20)

summary(lm(res2 ~ snpi, l20))
summary(lm(res ~ snp1, l20))
summary(lm(phen ~ snpi, l20))
summary(lm(phen ~ snp1, l20))
summary(lm(phen ~ snp2, l20))


write.table(subset(tab, select=c(cistrans, probe, snp1, snp2, snpi, p1, p2)), file="isnp_analysis.txt", row=T, col=T, qu=F)
