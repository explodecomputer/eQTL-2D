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
		s <- sign(x$E[5:8])

		r <- list()
		for(i in 2:3)
		{
			y <- subset(x, dataset == i)
			r[[i-1]] <- sum(s == sign(y$E[5:8])) == 4
		}

		out <- data.frame(r)
		names(out) <- c("r1", "r2")
		return(out)
	})
	p1 <- binom.test(n = nrow(els), x = sum(els$r1), p = 1/16)$p.value
	p2 <- binom.test(n = nrow(els), x = sum(els$r2), p = 1/16)$p.value
	p12 <- binom.test(n = nrow(els), x = sum(els$r1 & els$r2), p = 1/256)$p.value

	l <- list()
	l$els <- els
	l$p1 <- p1
	l$p2 <- p2
	l$p12 <- p12
	l$e1 <- 1/16 * nrow(els)
	l$e2 <- 1/16 * nrow(els)
	l$e12 <- 1/256 * nrow(els)
	return(l)

}

test4 <- function(noiad)
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


load("~/repo/eQTL-2D/analysis/noiad.RData")

els1 <- test1(noiad)
els2 <- test2(noiad)
els3 <- test3(noiad)
els4 <- test4(noiad)

bels1 <- test1(bnoiad)
bels2 <- test2(bnoiad)
bels3 <- test3(bnoiad)
bels4 <- test4(bnoiad)

dat <- rbind(
	makeTable(els1, bels1, "1"),
	makeTable(els2, bels2, "2"),
	makeTable(els3, bels3, "3"),
	makeTable(els4, bels4, "4")
)

dat


print(xtable(dat, digits = c(1, 1, 1, 1, 1, 2, 1, -2)), include.rownames=FALSE)
