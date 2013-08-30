#===========================================================#
#															#
# model_simulation.R										#
#															#
# simulations to test for an inflation in FPR when 			#
# selecting from a high initial p-value for the nested		#
# test. (i.e. 4df. vs 8df)									#
#															#
# joseph.powell@uq.edu.au  29.08.2013						#		
#															#
#===========================================================#



# Initial scan of correlation of 
# Jian's initial method


sim1.fun <- function(
	n=1000,
	e 	# Effect
	) {


	# create objects
	r <- c()
	pval <- c()
	a_d <- c()
	b_d <- c()
	pval_a <- c()
	pval_b <- c()

	b <- 0

	for(i in 1:1000) {

		a <- rnorm(1,0,e) 
		x <- rnorm(n)
		y <- rnorm(n)
		z <- x*a + y*b + rnorm(n)

		r[i] <- cor(x,y)

		reg 		<- summary(lm(z~x+y))
		a_d[i] 		<- coefficients(reg)[2,1]-a
		b_d[i] 		<- coefficients(reg)[3,1]
		pval_a[i] 	<- coefficients(reg)[2,4]
		pval_b[i] 	<- coefficients(reg)[3,4]
		pval[i] 	<- pf(reg$fstatistic[1], reg$fstatistic[2], reg$fstatistic[3], lower.tail=F)
	}

	out <- cbind(r,pval,a_d,b_d,pval_a,pval_b)	
	names(out) <- c("r","pval", "a_d", "b_d", "pval_a", "pval_b")
	return(out)
}


# Need to apply and have a testable outcome for lack of correlation at the tails. 



tmp=order(-log10(pval), decreasing=T)
out <- c()
for(i in 1:(length(pval)-29)) {
	out[i] <- cor(a_d[tmp[i:(i+29)]], b_d[tmp[i:(i+29)]])
}


# NCP of chi sq
# Peter's suggestion

sim2.fun <- function(
	thres.8df, 
	sim,
	NCP
	){


	#NCP <- n*(h2/(1-h2))
	chi8 <- array(0, sim)
	chi4 <- array(0, sim)
	p8 <- array(0, sim)
	p4 <- array(0, sim)

	for(i in 1:sim) {

		# simulate z^2 normal deviates
		zs <- rnorm(7, 0, 1)^2
		z8 <- rnorm(1, sqrt(NCP), 1)^2

		chi8[i] <- sum(c(zs,z8))
		chi4[i] <- sum(zs[1:4])

		p8[i] <- -log10(pchisq(chi8[i], df=8, lower.tail=F))
		p4[i] <- -log10(pchisq(chi4[i], df=4, lower.tail=F))
	}

	out <- as.data.frame(cbind(chi8, chi4, p8, p4))
	names(out) <- c("chi8", "chi4", "p8", "p4")

	thres.index <- which(out$p8 > thres.8df)
	l <- length(thres.index)	
	type1 <- length(which(out$p4[thres.index] > 1.3)) / l
	type2 <- (l/sim)
	list(type1=type1, type2=type2, ncp=NCP, fullresults=out)
}



tmp <- sim2.fun(16, 100000, 85)

#tmp$type2
#tmp$type1

ty1 <- c()
ty2 <- c()
for(i in seq(50, 150, 5)) {
	tmp <- sim2.fun(16, 100000, i)
	ty1[i] <- tmp$type1
	ty2[i] <- tmp$type2
}


plot(ty2,ty1,xlab="power (1st stage)", ylab="type 1 error rate 2ns stage")























