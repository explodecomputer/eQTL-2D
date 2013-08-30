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
	h2, 					# heritability of cis effect
	thres.8df, 
	n,
	sim
	
	){

	NCP <- n*(h2/(1-h2))
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
#	return(out)

	thres.index <- which(out$p8 > thres.8df)
	type1 <- length(which(out$p4[thres.index] > -log10(0.05))) / length(thres.index)
	list(above.thres=length(thres.index), type1=type1, ncp=NCP, fullresults=out)
}



tmp <- sim2.fun(0.09, 16, 846, 1000000)
tmp$ncp
tmp$type1

qchisq(10^-16, 8, lower.tail=F)



#################
# Peter's email.
#################

Dear both,

Re(1). I have thought about this a bit more and think that this can be approximated by simulation but without having to simulate individual trait values (i.e. can be done easily and fast). The 8df threshold is about 91. If we assume that the null model is that there is one additive cis-effect only (reasonable for a null model I think), then both the null and alternative (4 df test) can be done by simulating chisquare statistics (from normal deviates):

Null mode: chi(8) = sum(z^2), with z ~N(0,1) variates for 7 value, and for 1 value z ~ N(mu,1), with mu = sqrt(NCP). 

NCP is the chisquare non-centrality parameter pertaining to the cis-effect. If the cis-effects explains h^2 of phenotypic variance for the trait, then NCP = N * h^2 / (1 - h^2), with N the sample size (~ 500 in our case). I'm not sure what a reasonable value of NCP is for the study (maybe 50?) but we can get this from the previous BSGS publications. 

Reduced model: chi(4) = chi(8) - sum[z^2], where chi(8) is the already simulated value and the sum is over the simulated z-statistics from 5 to 8, including the z-statistics pertaining to the cis-effect.

So a simple simulation algorithm would be:
1. Set 8df threshold and set NCP
2. Simulate 8 z-statistics and perform test
3. If chi(8) > threshold then perform 4 df test
4. Count if 4df test statistic is above a 0.05 threshold for a 4df chisquare distribution (= 9.49)

Please let me know if you think that this would work and if it is worthwhile doing. I'm curious about the actual type-I error rate among  the 500 selected pairs.

Cheers,

Peter



















