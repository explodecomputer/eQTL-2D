# Simulate some intermediate frequency SNPs with known R
# Calculate empirical R
# Calculate sampling variance
# Calculate empirical R^8
# Calculate sampling variance of R^8


# r = (x11 - p1q1) / sqrt(p1q1p2q2)
# r = (x12 - p1q2) / sqrt(p1q1p2q2)

# given p1, q1, r we can calculate expected x11, x12, x21, x22


library(plyr)
library(ggplot2)


#=================================================================#


sampleSnp <- function(r, n, p1, q1)
{
	p2 <- 1 - p1
	q2 <- 1 - q1

	x11 <- r * sqrt(p1*q1*p2*q2) + p1*q1
	x12 <- p1*q2 - r * sqrt(p1*q1*p2*q2)
	x21 <- p2*q1 - r * sqrt(p1*q1*p2*q2)
	x22 <- r * sqrt(p1*q1*p2*q2) + p2*q2

	xs <- c(x11, x12, x21, x22)

	if(any(xs < 0))
	{
		return(NA)
	}

	x <- sample(1:4, n*2, replace=T, prob=xs)
	
	p1 <- (sum(x == 1) + sum(x == 2)) / (n*2)
	p2 <- (sum(x == 3) + sum(x == 4)) / (n*2)
	q1 <- (sum(x == 1) + sum(x == 3)) / (n*2)
	q2 <- (sum(x == 2) + sum(x == 4)) / (n*2)
	
	r_obs <- (sum(x == 1)/(n*2) - p1*q1) / sqrt(p1*q1*p2*q2)
	
	return(r_obs)
}


doSim <- function(param)
{
	for(i in 1:nrow(param))
	{
		param$r_obs[i] <- with(param, sampleSnp(r[i], n[i], p1[i], q1[i]))
	}
	return(param)
}


#=================================================================#

# Create parameters
param <- expand.grid(
	r = c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99), 
	n = c(50, 900, 1200, 20000), 
	p1 = c(0.1, 0.2, 0.3, 0.4, 0.5), 
	q1 = c(0.1, 0.2, 0.3, 0.4, 0.5), 
	sim = 1:50, r_obs = NA
)
param <- subset(param, p1 >= q1)

# Run simulations
param <- doSim(param)


#=================================================================#


# Summarise results
p <- subset(param, !is.na(r_obs), select=c(r, r_obs, n, p1, q1))
ps <- rbind(p, p, p, p, p)
ps$pow <- rep(c(1,2,4,6,8), each=nrow(p))
ps$r_obs <- ps$r_obs^ps$pow

ps_summary <- ddply(subset(ps), .(pow, r, n), summarise, 
	N    = length(r_obs),
	sd   = sd(r_obs),
	se   = sd(r_obs) / sqrt(length(r_obs)),
	qu1  = quantile(r_obs, prob=0.25, na.rm=T),
	mean = mean(r_obs),
	qu3  = quantile(r_obs, prob=0.75, na.rm=T)	
)

subset(ps_summary, r==0.5 & pow==2)

#=================================================================#

# Plot
ggplot(subset(ps_summary, n > 100), aes(x = pow, y = sd)) +
geom_point() +
geom_line(aes(colour=factor(n))) +
facet_grid(. ~ r) + 
labs(x = "Power term", y = "Standard deviation", colour = "Sample size") +
scale_colour_brewer(type="qual")
ggsave(file="~/repo/eQTL-2D/analysis/images/ld_sampling_sd.pdf", width=10, height=6)

ggplot(subset(ps_summary, n > 100), aes(x = pow, y = mean)) +
geom_point() +
geom_line(aes(colour=factor(n))) +
facet_grid(. ~ r) + 
labs(x = "Power term", y = "Mean", colour = "Sample size") +
scale_colour_brewer(type="qual")
ggsave(file="~/repo/eQTL-2D/analysis/images/ld_sampling_mean.pdf", width=10, height=6)

ggplot(subset(ps_summary, r != 0), aes(x = pow, y = sd^2 / mean)) +
geom_point() +
geom_line(aes(colour=factor(n))) +
facet_grid(. ~ r) + 
labs(x = "Power term", y = "Coefficient of variance", colour = "Sample size") +
scale_colour_brewer() +
theme_bw()

ggplot(subset(ps, r != 0 & n > 100), aes(y = r_obs, x = factor(pow))) +
geom_boxplot(aes(fill = factor(n))) +
facet_grid( . ~ r)

ggplot(subset(ps, r != 0 & n > 100), aes(y = r_obs, x = factor(pow))) +
geom_point(aes(colour = factor(n)), position="dodge") +
facet_grid( . ~ r)
