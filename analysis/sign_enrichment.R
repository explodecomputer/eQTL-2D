library(noia)
library(plyr)


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
	l$p1 <- p2
	l$p12 <- p12
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
	l$p1 <- p2
	l$p12 <- p12
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
			r[[i-1]] <- sum(s == sign(y$E[5:8]))
		}

		out <- data.frame(r)
		names(out) <- c("r1", "r2")
		return(out)
	})
	p1 <- binom.test(n = nrow(els), x = sum(els$r1 == 4), p = 1/16)$p.value
	p2 <- binom.test(n = nrow(els), x = sum(els$r2 == 4), p = 1/16)$p.value
	p12 <- binom.test(n = nrow(els), x = sum(els$r1 == 4 & els$r2 == 4), p = 1/256)$p.value
	l <- list()
	l$els <- els
	l$p1 <- p1
	l$p1 <- p2
	l$p12 <- p12
	return(l)

}

test4 <- function(noiad)
{
	els <- ddply(noiad, .(interaction), function(x)
	{
		x <- mutate(x)
		a <- subset(x, dataset == 1)
		sig <- which(x$p[5:8] < 0.05)
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
	l$p1 <- p2
	l$p12 <- p12
	return(l)
}


load("~/repo/eQTL-2D/analysis/noiad.RData")


els1 <- test1(noiad)
els2 <- test2(noiad)
els3 <- test3(noiad)
els4 <- test4(noiad)

els1$p12
els2$p12
els3$p12
els4$p12

bels1 <- test1(bnoiad)
bels2 <- test2(bnoiad)
bels3 <- test3(bnoiad)
bels4 <- test4(bnoiad)

bels1$p12
bels2$p12
bels3$p12
bels4$p12
