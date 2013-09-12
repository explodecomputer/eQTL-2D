# gcta64 --bfile data/arichu_selected_snps --keep arichu_discovery.ids --ld qtls.txt --ld-wind 1000 --ld-sig 0.05 --out arichu

a <- read.table("~/repo/eQTL-2D/ld_simulations/arichu.rsq.ld", header=T)
head(a)

hist(a$max_rsq)

table(a$max_rsq > 0.95)

b <- subset(a, max_rsq > 0.9)
dim(b)
head(b)

write.table(b, "arichu.rsq.ld.9.discovery", row=F, col=T, qu=F)

snplist <- with(b, unique(c(as.character(target_SNP), as.character(max_rsq_snp))))
length(snplist)
dim(b)

write.table(snplist, "arichu.rsq.ld.9.discovery.snps", row=F, col=F, qu=F)

plink --bfile data/arichu_selected_snps --keep arichu_replication.ids --extract arichu.rsq.ld.9.discovery.snps --make-bed --out data/arichu_replication_9

library(snpStats)
rawdata <- read.plink("data/arichu_replication_9")

xi <- match(b$target_SNP, rawdata$map$snp.name)
yi <- match(b$max_rsq_snp, rawdata$map$snp.name)

length(xi)
length(yi)

x <- rawdata$genotypes[,xi]
y <- rawdata$genotypes[,yi]

win <- 500
breaks <- floor(length(xi)/win)
res <- array(0, length(xi))
for(i in 1:breaks)
{
	cat(i, "of", breaks, "\n")
	index <- 1:win + win * (i-1)
	X <- x[,index]
	Y <- y[,index]
	mat <- ld(X, Y, stats="R.squared")
	res[index] <- diag(mat)
}
index <- c((breaks*win+1) : length(xi))
X <- x[,index]
Y <- y[,index]
mat <- ld(X, Y, stats="R.squared")
res[index] <- diag(mat)

hist(res)
hist(b$max_rsq)
plot(res ~ b$max_rsq, ylim=c(0.65, 1), xlim=c(0.65,1))
abline(lm(res~b$max_rsq))
summary(lm(res ~ b$max_rsq))

b$replication_rsq <- res

mean(subset(b, max_rsq > 0.99)$max_rsq)
mean(subset(b, max_rsq > 0.99)$replication_rsq)

mean(subset(b, max_rsq > 0.9)$max_rsq)
mean(subset(b, max_rsq > 0.9)$replication_rsq)





# gcta64 --bfile data/1kg_selected_snps --keep 1kg_discovery.ids --ld 1kg_qtls.txt --ld-wind 1000 --ld-sig 0.05 --out 1kg_discovery

a <- read.table("~/repo/eQTL-2D/ld_simulations/1kg_discovery.rsq.ld", header=T)
head(a)

hist(a$max_rsq)

table(a$max_rsq > 0.95)

b <- subset(a, max_rsq > 0.9)
dim(b)
head(b)

write.table(b, "1kg_discovery.rsq.ld.9.discovery", row=F, col=T, qu=F)

snplist <- with(b, unique(c(as.character(target_SNP), as.character(max_rsq_snp))))
length(snplist)
dim(b)

write.table(snplist, "1kg_discovery.rsq.ld.9.discovery.snps", row=F, col=F, qu=F)

# plink --bfile data/1kg_selected_snps --keep 1kg_replication.ids --extract 1kg_discovery.rsq.ld.9.discovery.snps --make-bed --out data/1kg_replication_9

library(snpStats)
rawdata <- read.plink("data/1kg_replication_9")

xi <- match(b$target_SNP, rawdata$map$snp.name)
yi <- match(b$max_rsq_snp, rawdata$map$snp.name)

length(xi)
length(yi)

x <- rawdata$genotypes[,xi]
y <- rawdata$genotypes[,yi]

win <- 500
breaks <- floor(length(xi)/win)
res <- array(0, length(xi))
for(i in 1:breaks)
{
	cat(i, "of", breaks, "\n")
	index <- 1:win + win * (i-1)
	X <- x[,index]
	Y <- y[,index]
	mat <- ld(X, Y, stats="R.squared")
	res[index] <- diag(mat)
}
index <- c((breaks*win+1) : length(xi))
X <- x[,index]
Y <- y[,index]
mat <- ld(X, Y, stats="R.squared")
res[index] <- diag(mat)
res[res > 0.999] <- NA

hist(res)
hist(b$max_rsq)
plot(res ~ b$max_rsq)
abline(lm(res~b$max_rsq))
summary(lm(res ~ b$max_rsq))

b$replication_rsq <- res

save(b, file="ld_simulations.RData")

mean(subset(b, max_rsq > 0.9)$max_rsq)
mean(subset(b, max_rsq > 0.9)$replication_rsq, na.rm=T)

mean(subset(b, max_rsq > 0.9)$max_rsq)^8
mean(subset(b, max_rsq > 0.9)$replication_rsq, na.rm=T)^8


lm(replication_rsq ~ max_rsq, data=subset(b, max_rsq > 0.95))



plot(replication_rsq ~ max_rsq, data=subset(b, max_rsq > 0.9))
abline(lm(replication_rsq ~ max_rsq, data=subset(b, max_rsq > 0.95)), col="red")
abline(a=0,b=1, col="blue")

dim(b)

b2 <- subset(b, max_rsq > 0.9)
b2$pow <- 2
b2$thresh <- b2$max_rsq > 0.95
b4 <- b2
b4$pow <- 4
b4$max_rsq <- b4$max_rsq^2
b4$replication_rsq <- b4$replication_rsq^2
b6 <- b2
b6$max_rsq <- b6$max_rsq^3
b6$replication_rsq <- b6$replication_rsq^3
b6$pow <- 6
b8 <- b2
b8$max_rsq <- b8$max_rsq^4
b8$replication_rsq <- b8$replication_rsq^4
b8$pow <- 8

b_all <- rbind(b2, b4, b6, b8)


library(ggplot2)


ggplot(b_all, aes(x=max_rsq, y=replication_rsq)) +
geom_point() +
geom_abline(intercept=0, slope=1, colour="blue") +
geom_smooth(colour="red") +
facet_wrap(~ pow) +
labs(y = "Replication", x = "Discovery")
ggsave("~/repo/eQTL-2D/analysis/images/ld_ascertained.png")


library(plyr)

a <- ddply(b_all, .(pow), function(x)
{
	x <- mutate(x)
	m <- mean(x$max_rsq - x$replication_rsq, na.rm=T)
	return(m)
})


a <- ddply(b_all, .(pow), function(x)
{
	x <- mutate(x)
	m <- mean((x$replication_rsq - x$max_rsq) / x$max_rsq, na.rm=T)

	return(m)
})
a

b_all$rdiff <- b_all$replication_rsq - b_all$max_rsq
b_all$prdiff <- (b_all$replication_rsq - b_all$max_rsq) / b_all$max_rsq
ggplot(b_all, aes(x=factor(pow), y=prdiff)) + 
geom_boxplot() +
labs(y = expression( (r[R]^x - r[D]^x)/r[D]^x) , x = "Power term")
ggsave("~/repo/eQTL-2D/analysis/images/ld_reduction.pdf")



