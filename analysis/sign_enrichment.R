library(noia)
library(plyr)
library(xtable)


# Largest discovery has the same sign in replication

test1 <- function(noiad)
{
	els <- ddply(noiad, .(interaction), function(x){

		x <- mutate(x)
		a <- subset(x, dataset == 1)
		m <- which.min(a$p[5:8])
		s <- sign(x$E[m])

		r <- list()
		for(i in 2:3)
		{
			y <- subset(x, dataset == i)
			r[[i-1]] <- s == sign(y$E[m])
		}

		out <- data.frame(r)
		names(out) <- c("r1", "r2")
		return(out)
	})
	p1 <- binom.test(n = nrow(els), x = sum(els$r1), p = 0.5)$p.value
	p2 <- binom.test(n = nrow(els), x = sum(els$r2), p = 0.5)$p.value
	p12 <- binom.test(n = nrow(els), x = sum(els$r2 & els$r1), p = 0.25)$p.value
	l <- list()
	l$els <- els
	l$p1 <- p1
	l$p2 <- p2
	l$p12 <- p12
	l$e1 <- 0.5 * nrow(els)
	l$e2 <- 0.5 * nrow(els)
	l$e12 <- 0.25 * nrow(els)
	return(l)
}

test2 <- function(noiad)
{
	els <- ddply(noiad, .(interaction), function(x){

		x <- mutate(x)
		a <- subset(x, dataset == 1)
		m <- which.min(a$p[5:8])
		s <- sign(x$E[m])

		r <- list()
		for(i in 2:3)
		{
			y <- subset(x, dataset == i)
			m1 <- which.min(y$p[5:8])
			r[[i-1]] <- s == sign(y$E[m]) & m == m1
		}

		out <- data.frame(r)
		names(out) <- c("r1", "r2")
		return(out)
	})
	p1 <- binom.test(n = nrow(els), x = sum(els$r1), p = 1/8)$p.value
	p2 <- binom.test(n = nrow(els), x = sum(els$r2), p = 1/8)$p.value
	p12 <- binom.test(n = nrow(els), x = sum(els$r2 & els$r1), p = 1/64)$p.value
	l <- list()
	l$els <- els
	l$p1 <- p1
	l$p2 <- p2
	l$p12 <- p12
	l$e1 <- 1/8 * nrow(els)
	l$e2 <- 1/8 * nrow(els)
	l$e12 <- 1/64 * nrow(els)
	return(l)
}

test3 <- function(noiad)
{
	els <- ddply(noiad, .(interaction), function(x)
	{
		x <- mutate(x)
		a <- subset(x, dataset == 1)
		sig <- which(x$p[5:8] < 0.1)
		s <- sign(x$E[sig])

		r <- list()
		for(i in 2:3)
		{
			y <- subset(x, dataset == i)
			r[[i-1]] <- s == sign(y$E[sig])
		}
		r[[3]] <- length(sig)

		out <- data.frame(r)
		names(out) <- c("r1", "r2", "l")
		return(out)
	})
	p1 <- binom.test(n = nrow(els), x = sum(els$r1), p = 0.5)$p.value
	p2 <- binom.test(n = nrow(els), x = sum(els$r2), p = 0.5)$p.value
	p12 <- binom.test(n = nrow(els), x = sum(els$r2 & els$r1), p = 0.25)$p.value
	l <- list()
	l$els <- els
	l$p1 <- p1
	l$p2 <- p2
	l$p12 <- p12
	l$e1 <- 0.5 * nrow(els)
	l$e2 <- 0.5 * nrow(els)
	l$e12 <- 0.25 * nrow(els)
	return(l)
}

testMult <- function(noiad)
{
	els <- ddply(noiad, .(interaction), function(x)
	{
		x <- mutate(x)
		a <- subset(x, dataset == 1)
		s <- sign(x$E[5:8])

		r <- list()
		for(i in 2:3)
		{
			y <- subset(x, dataset == i)
			r[[i-1]] <- sum(s == sign(y$E[5:8]))
		}

		out <- data.frame(r)
		names(out) <- c("r1", "r2")
		return(out)
	})


	l <- list()
	l$els <- els
	l$p1 <- NA
	l$p2 <- NA
	l$p12 <- NA
	l$e1 <- 1/16 * nrow(els)
	l$e2 <- 1/16 * nrow(els)
	l$e12 <- 1/256 * nrow(els)
	return(l)

}

makeTable <- function(els, bels, nom)
{
	# 1, 2, and 1+2
	# number of trials, number of successes, p value

	n1 <- rep(nrow(els$els), 3)
	x1 <- c(sum(els$els$r1), sum(els$els$r2), sum(els$els$r1 & els$els$r2))
	e1 <- c(els$e1, els$e2, els$e12)
	p1 <- c(els$p1, els$p2, els$p12)

	n2 <- rep(nrow(bels$els), 3)
	x2 <- c(sum(bels$els$r1), sum(bels$els$r2), sum(bels$els$r1 & bels$els$r2))
	e2 <- c(bels$e1, bels$e2, bels$e12)
	p2 <- c(bels$p1, bels$p2, bels$p12)

	a1 <- data.frame(Test = rep(nom, 3), Interactions = rep("All", 3), Dataset = c("EGCUT", "Fehrmann", "Both"), n = n1, Expected = e1, Observed = x1, p = p1)
	a2 <- data.frame(Test = rep(nom, 3), Interactions = rep("Significant", 3), Dataset = c("EGCUT", "Fehrmann", "Both"), n = n2, Expected = e2, Observed = x2, p = p2)

	return(rbind(a1,a2))
}

multTest <- function(obs, pval)
{
	k <- length(obs)
	n <- sum(obs)
	chisq <- -2 * sum(obs * log(pval / (obs/n)))
	q <- 1 + sum(1/pval - 1) / (6 * n * (k - 1))
	pchisq(chisq / q, k-1, lower.tail=FALSE)
}

multTable <- function(noiad, bnoiad)
{
	els <- testMult(noiad)
	bels <- testMult(bnoiad)

	ecounts <- rbind(
		c(table(els$els$r1)),
		c(table(els$els$r2)),
		c(table(els$els$r1) + table(els$els$r2))
	)

	bcounts <- rbind(
		c(table(bels$els$r1)),
		c(table(bels$els$r2)),
		c(table(bels$els$r1) + table(bels$els$r2))
	)

	ep <- array(0, 3)

	expected <- c(1,4,6,4,1) / 16

	ne <- c(nrow(els$els), nrow(els$els), 2 * nrow(els$els))
	nb <- c(nrow(bels$els), nrow(bels$els), 2 * nrow(bels$els))

	ep[1] <- multTest(ecounts[1,], expected)
	ep[2] <- multTest(ecounts[2,], expected)
	ep[3] <- multTest(ecounts[3,], expected)

	bp <- array(0, 3)

	bp[1] <- multTest(bcounts[1,], expected)
	bp[2] <- multTest(bcounts[2,], expected)
	bp[3] <- multTest(bcounts[3,], expected)

	dat <- data.frame(c("Expected", rep("All", 3), rep("Significant", 3)), c(NA, rep(c("EGCUT", "Fehrmann", "Combined"), 2)), c(NA, ne, nb), rbind(expected, ecounts / ne, bcounts / nb), c(NA, ep, bp))

	names(dat) <- c("Interactions", "n", "Dataset", "0", "1", "2", "3", "4", "p")
	return(dat)
}



load("~/repo/eQTL-2D/analysis/noiad.RData")

els1 <- test1(noiad)
els2 <- test2(noiad)
els3 <- test3(noiad)

bels1 <- test1(bnoiad)
bels2 <- test2(bnoiad)
bels3 <- test3(bnoiad)

dat <- rbind(
	makeTable(els1, bels1, "1"),
	makeTable(els2, bels2, "2"),
	makeTable(els3, bels3, "3")
)

dat


print(xtable(dat, digits = c(1, 1, 1, 1, 1, 2, 1, -2)), include.rownames=FALSE)

dat2 <- multTable(noiad, bnoiad)

print(xtable(dat2, digits = c(0, 0, 0, 0, 2, 2, 2, 2, 2, -2)), include.rownames=FALSE)

