---
title: Checking lambda estimates for 4df test
author: Gibran Hemani
---

Here are a few different ways that could be used to estimate lambda, that should all give the same answer:


```{r}
# convert p-value to 1df chisq:
# Get the z score from the p-value, get the median, square it
estlambda1 <- function(x)
{
	da <- qnorm(1 - x/2)
	median(da, na.rm=TRUE)^2/qchisq(0.5, 1)
}

# convert p-value to 1df chisq:
# get the chisq value directly, get the median
estlambda2 <- function(x)
{
	da <- qchisq(x, 1, low=FALSE)
	median(da, na.rm=TRUE)/qchisq(0.5, 1)
}

# convert p-value to 1df chisq:
# Get the median p-value, convert to z score, square it
estlambda3 <- function(x)
{
	da <- qnorm(1 - median(x, na.rm=TRUE)/2)
	da^2/qchisq(0.5, 1)
}

# convert p-value to 1df chisq:
# get the median p-value, convert to chisq value
estlambda4 <- function(x)
{
	da <- qchisq(median(x, na.rm=TRUE), 1, low=FALSE)
	da/qchisq(0.5, 1)
}

a <- runif(100000)

estlambda1(a)
estlambda2(a)
estlambda3(a)
estlambda4(a)
```


What happens to the p-values in a 4 df test if not all genotype classes are present (e.g. small sample size and low maf)


```{r}
library(dplyr)

make_geno <- function(nid, af)
{
	rbinom(nid, 2, af)
}

g1 <- make_geno(850, 0.5)
g2 <- make_geno(850, 0.005)
y <- rnorm(850)

table(g1, g2) %>% length
mod1 <- lm(y ~ as.factor(g1) + as.factor(g2))
mod2 <- lm(y ~ as.factor(g1):as.factor(g2))
anova(mod1, mod2)
```




###


library(simulateGP)

g <- make_geno(800, 50001, 0.5)

run <- function(g)
{
	n <- nrow(g)
	m <- ncol(g)
	y <- rnorm(n)
	out <- list()
	for(i in 1:(m-1))
	{
		if(i %% 1000 == 0) message(i)
		mod1 <- anova(lm(y ~ g[,i]))
		mod2 <- lm(y ~ as.factor(g[,i]):as.factor(g[,m]))
		mod3 <- lm(y ~ as.factor(g[,i])+as.factor(g[,m]))
		mod4 <- anova(mod3,mod2)
		out[[i]] <- tibble(df1=mod1$F[1], df2=mod4$F[1])
	}
	bind_rows(out) %>% return()
}



ftest4df <- function(x1, x2, y)
{
	my <- mean(y)
	ty <- table(x1, x2)
	my_fac <- tapply(y, list(x1, x2), mean)
	mean_col <- tapply(y, x1, mean)
	mean_row <- tapply(y, x2, mean)
	mean_col <- rbind(mean_col, mean_col, mean_col)
	mean_row <- cbind(mean_row, mean_row, mean_row)
	SSI <- sum((my_fac - mean_col - mean_row + my)^2 * ty)
	SSB <- sum((my_fac - my)^2 * ty)
	SS <- sum((y - my)^2)
	MSW <- (SS - SSB) / (length(y) - 9)
	f8df <- (SSB/8) / MSW
	f4df <- (SSI/4) / MSW
	return(c(f8df, f4df))

}





	((SSM1 + SSM2) / 4) / (SSB - (SSM1 + SSM2) / 995)


	# Marginal SSB




	anova(lm(y ~ as.factor(x1)))
	anova(lm(y ~ as.factor(x2)))

	anova(lm(y ~ as.factor(x1):as.factor(x2)))


	(SSM1 / 2) / ((SS - SSM1) / 997)
	(SSM2 / 2) / ((SS - SSM2) / 997)


	summary(lm(y ~ as.factor(x1)+as.factor(x2)))

	s <- SSM1 + SSM2
	(SSM / 4) / ((SS-SSM) / 995)


	(((SS - SSM) - (SS - SSB)) / 8)
	((SS-SSM - (SS-SSB)) / (SS-SSB)) / 8
	anova(
		lm(y ~ as.factor(x1)+as.factor(x2)),
		lm(y ~ as.factor(x1)*as.factor(x2))
	)


	(((SS - SSM1) - (SS - SSB)) / 6)
	anova(
		lm(y ~ as.factor(x1)),
		lm(y ~ as.factor(x1)*as.factor(x2))
	)

	(((SS - SSM2) - (SS - SSB)) / 6)
	anova(
		lm(y ~ as.factor(x2)),
		lm(y ~ as.factor(x1)*as.factor(x2))
	)


	(((SS - SSM) - (SS - SSB)) / 4)
	anova(
		lm(y ~ as.factor(x1)+as.factor(x2)),
		lm(y ~ as.factor(x1)*as.factor(x2))
	)


	SSB

	((SSB - s) / 4) / (SSB / 995)


}

s <- 5.11 + 0.28
(s / 4) / ((997 - s) / 995)



ftest4df <- function(x1, x2, y)
{
	my <- mean(y)
	ty <- table(x1, x2)
	my_fac <- tapply(y, list(x1, x2), mean)
	mean_col <- tapply(y, x1, mean)
	mean_row <- tapply(y, x2, mean)
	t1 <- table(x1)
	t2 <- table(x2)
	SSM1 <- sum((mean_col - my)^2 * t1)
	SSM2 <- sum((mean_row - my)^2 * t2)
	SSM <- SSM1 + SSM2
	SSB <- sum((my_fac - my)^2 * ty)
	SS <- sum((y - my)^2)
	SSI <- SSB - SSM
	MSW <- (SS - SSB) / (length(y) - 9)
	f8df <- (SSB/8) / MSW
	f4df <- (SSI/4) / MSW
	fmarg <- (SSM/4) / MSW
	return(c(f8df, f4df, fmarg))
}



X <- x1 + x2 * 3 + 1
table(X)

m1 <- rep(0, 3)
m2 <- rep(0, 3)
for(i in 1:length(X))
{
	m1[(X[i]-1) %% 3 + 1] <- m1[(X[i]-1) %% 3 + 1] + 1
	m2[floor((X[i]-1)/3)+1] <- m2[floor((X[i]-1)/3)+1] + 1
}

m1
m2

table(x1)
table(x2)



summary(lm(y ~ as.factor(x1)+as.factor(x2)))
summary(lm(y ~ as.factor(x1)*as.factor(x2)))
ftest4df(x1,x2,y)
anova(
	lm(y ~ as.factor(x1)+as.factor(x2)),
	lm(y ~ as.factor(x1)*as.factor(x2))
)




out <- list()
for(i in 1:100)
{
	message(i)
	af1 <- runif(1, 0.3, 0.7)
	af2 <- runif(1, 0.3, 0.7)
	x1 <- make_geno(1000, 1, af1)
	x2 <- make_geno(1000, 1, af2)
	# b <- rnorm(1)
	b <- 0
	y <- b * x1 * x2 + rnorm(1000)


	o <- ftest4df(x1, x2, y)

	m1 <- lm(y ~ as.factor(x1):as.factor(x2))
	m2 <- lm(y ~ as.factor(x1)+as.factor(x2))
	
	out[[i]] <- tibble(
		b=b,
		af1=af1,
		af2=af2,
		f8df=o[1],
		f4df=o[2],
		ao8=anova(m1)$F[1],
		ao4=anova(m2, m1)$F[2]
	)
}

out <- bind_rows(out)
plot(out$f4df ~ out$ao4)


ou <- run(g)


a <- read.table("temp.txt", he=F)
a$x1 <- a$V1 %% 3
a$x2 <- floor(a$V1 / 3)
cor(a$x1, a$x2)

table(a$V1)
table(a$x1)
table(a$x2)

ftest4df(a$x1, a$x2, a$V2)

anova(
	lm(V2 ~ as.factor(x1)+as.factor(x2), data=a),
	lm(V2 ~ as.factor(x1)*as.factor(x2), data=a)
)

cor(a$x1, a$x2)


anova(lm(V2 ~ as.factor(x1)+as.factor(x2), data=a))
summary(lm(V2 ~ as.factor(x1)+as.factor(x2), data=a))


####



./episcan -A ~/mr-eve/eQTL-2D/null_replication/data/scratch/MBNL1_rs67903230_rs13069559/rs13069559_merge_disc.egu out -t i -F 0.00000001 -I 0.0000001 -1 1 -2 2 -f ~/mr-eve/eQTL-2D/null_replication/data/scratch/MBNL1_rs67903230_rs13069559/1_disc.fam -T 1


./episcan -De temp.ped temp.map ~/mr-eve/eQTL-2D/null_replication/data/scratch/MBNL1_rs67903230_rs13069559/rs13069559_merge_disc.egu

plink --file temp --recode A --out temp

cut -d " " -f 1-500 temp.raw > temp2.raw




ftest4df <- function(x1, x2, y)
{
	my <- mean(y)
	ty <- table(x1, x2)
	my_fac <- tapply(y, list(x1, x2), mean)
	mean_col <- tapply(y, x1, mean)
	mean_row <- tapply(y, x2, mean)
	t1 <- table(x1)
	t2 <- table(x2)
	SSM1 <- sum((mean_col - my)^2 * t1)
	SSM2 <- sum((mean_row - my)^2 * t2)
	SSM <- SSM1 + SSM2
	SSB <- sum((my_fac - my)^2 * ty)
	SS <- sum((y - my)^2)
	SSI <- SSB - SSM
	MSW <- (SS - SSB) / (length(y) - 9)
	f8df <- (SSB/8) / MSW
	f4df <- (SSI/4) / MSW
	fmarg <- (SSM/4) / MSW
	return(c(f8df, f4df, fmarg))
}


library(dplyr)
library(data.table)
a <- fread("temp2.raw", he=T)
fam <- fread("~/mr-eve/eQTL-2D/null_replication/data/scratch/MBNL1_rs67903230_rs13069559/1_disc.fam",)
stopifnot(all(a$FID == fam$V1))
a <- a[,-c(1:6)]

out <- fread("~/mr-eve/eQTL-2D/null_replication/data/scratch/MBNL1_rs67903230_rs13069559/1_disc_out0", he=FALSE)

out$f1 <- NA
out$f2 <- NA
out$p4 <- NA

fam$V6[is.na(fam$V6)] <- mean(fam$V6, na.rm=T)

res <- list()
a <- as.data.frame(a)
x1 <- a$rs13069559_G
y <- unlist(fam$V6)
for(i in 2:ncol(a))
{
	message(i)
	x2 <- unlist(a[,i, drop=T])
	j <- i - 1
	out$f1[j] <- ftest4df(x1, x2, y)[2]
	mod <- anova(
		lm(y ~ as.factor(x1) + as.factor(x2)),
		lm(y ~ as.factor(x1) : as.factor(x2))
	)
	out$f2[j] <- mod$F[2]
	out$p4[j] <- mod$P[2]
}



out2 <- subset(out, V3 == 8)

cor(out2$f1, out2$f2, use="pair")
cor(out2$V6, out2$f2, use="pair")

summary(lm(f1 ~ f2, out2))



out2$p42 <- pf(out2$f1, 4, out2$V4, low=FALSE)

cor(out2$p42, out2$p4, use="pair")


out2 <- subset(out2, !duplicated(V2))


estlambda1 <- function(x)
{
	da <- qnorm(1 - x/2)
	median(da, na.rm=TRUE)^2/qchisq(0.5, 1)
}

# convert p-value to 1df chisq:
# get the chisq value directly, get the median
estlambda2 <- function(x)
{
	da <- qchisq(x, 1, low=FALSE)
	median(da, na.rm=TRUE)/qchisq(0.5, 1)
}

# convert p-value to 1df chisq:
# Get the median p-value, convert to z score, square it
estlambda3 <- function(x)
{
	da <- qnorm(1 - median(x, na.rm=TRUE)/2)
	da^2/qchisq(0.5, 1)
}

# convert p-value to 1df chisq:
# get the median p-value, convert to chisq value
estlambda4 <- function(x)
{
	da <- qchisq(median(x, na.rm=TRUE), 1, low=FALSE)
	da/qchisq(0.5, 1)
}





estlambda4(out2$p4)
estlambda3(out2$p4)
estlambda2(out2$p4)
estlambda1(out2$p4)

median(out2$V6) / qf(0.5, 4, 840, low=F)

median(out2$V6) / qchisq(0.5, 4) * 4

median(out2$V6 / qchisq(0.5, 4) * 4)

