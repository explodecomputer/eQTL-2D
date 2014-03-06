library(lattice)
library(latticeExtra)
library(gridExtra)

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



snp1 <- bsgs[,which(colnames(bsgs) == "rs8106959")]
snp2 <- bsgs[,which(colnames(bsgs) == "rs6718480")]

1 - p[which(colnames(bsgs) == "rs8106959")]
1 - p[which(colnames(bsgs) == "rs6718480")]


(1 - 0.82)^2 * 0.42^2 * 450


subset(sig, snp2 == "rs6718480" & snp1 == "rs8106959")$gcs
subset(sig, snp2 == "rs6718480" & snp1 == "rs8106959")$gcm


all(rownames(bsgs) == a$V2) 


phen <- phenlist[[1]][,which(colnames(phenlist[[1]]) == "ILMN_1786426")]

res <- residuals(lm(phen ~ snp))
res2 <- residuals(lm(phen ~ snp1))

gp1 <- tapply(phen, list(snp1, snp2), mean)
gp1 <- gp1 - min(gp1)
gp2 <- tapply(res, list(snp1, snp2), mean)
gp2 <- gp2 - min(gp2)
gp3 <- tapply(res2, list(snp1, snp2), mean)
gp3 <- gp3 - min(gp3)

plot3dGp(gp1, z=30)
plot3dGp(gp2, z=30)
plot3dGp(gp3, z=30)


pdf("tmem_19x6.pdf")
plot3dGp(gp1, z=30)
plot3dGp(gp2, z=30)
dev.off()


plot3dGp(gp3, z=30)


summary(lm(res ~ as.factor(snp1) * as.factor(snp2)))
summary(lm(res2 ~ as.factor(snp1) * as.factor(snp2)))
summary(lm(phen ~ as.factor(snp1) * as.factor(snp2)))
cor(snp1, snp2)

mod1 <- lm(res2 ~ as.factor(snp1) * as.factor(snp2))
mod2 <- lm(res2 ~ as.factor(snp1) + as.factor(snp2))
anova(mod2, mod1)

mod1 <- lm(res2 ~ as.factor(snp) * as.factor(snp2))
mod2 <- lm(res2 ~ as.factor(snp) + as.factor(snp2))
anova(mod2, mod1)

mod1 <- lm(phen ~ as.factor(snp1) * as.factor(snp2))
mod2 <- lm(phen ~ as.factor(snp1) + as.factor(snp2))
anova(mod1, mod2)

summary(lm(res ~ as.factor(snp1) : as.factor(snp2)))


library(noia)

linearRegression(rep(phen, 3), cbind(rep(snp1,3), rep(snp2,3))+1)
linearRegression(rep(res, 5), cbind(rep(snp1,5), rep(snp2,5))+1)
linearRegression(res, cbind(snp1, snp2)+1, ref="P2")


