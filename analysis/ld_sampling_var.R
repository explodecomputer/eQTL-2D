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


doSimTwo <- function(param)
{
	for(i in 1:nrow(param))
	{
		param$r_obs1[i] <- with(param, sampleSnp(r1[i], n[i], p1[i], q1[i]))
		param$r_obs2[i] <- with(param, sampleSnp(r2[i], n[i], y1[i], z1[i]))
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
	sim = 1:20, 
	r_obs = NA
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



#=================================================================#
#=================================================================#



# Two loci
paramtwo <- expand.grid(
	r1 = c(0.75, 0.8, 0.85, 0.9, 0.95), 
	r2 = c(0.75, 0.8, 0.85, 0.9, 0.95), 
	n = c(891, 1240), 
	p1 = c(0.1, 0.3, 0.5), 
	q1 = c(0.1, 0.3, 0.5), 
	y1 = c(0.1, 0.3, 0.5), 
	z1 = c(0.1, 0.3, 0.5), 
	sim = 1:20, 
	r_obs1 = NA,
	r_obs2 = NA
)
paramtwo <- subset(paramtwo, p1 >= q1 & y1 >= z1 & r1 >= r2)

paramtwo <- doSimTwo(paramtwo)


paramtwo$r1_1 <- paramtwo$r_obs1
paramtwo$r1_2 <- paramtwo$r_obs2

paramtwo$r2_1 <- paramtwo$r_obs1^2
paramtwo$r2_2 <- paramtwo$r_obs2^2

p_r <- paramtwo
p_r$rx <- paramtwo$r1_2
p_r$type <- "Marginal"
p_r$x <- 1

p_a <- paramtwo
p_a$rx <- paramtwo$r2_2
p_a$type <- "Marginal"
p_a$x <- 2

p_d <- paramtwo
p_d$rx <- paramtwo$r2_2^2
p_d$type <- "Marginal"
p_d$x <- 4

p_aa <- paramtwo
p_aa$rx <- paramtwo$r2_1 * paramtwo$r2_2
p_aa$type <- "Epistatic"
p_aa$x <- 4

p_ad <- paramtwo
p_ad$rx <- paramtwo$r2_1^2 * paramtwo$r2_2
p_ad$type <- "Epistatic"
p_ad$x <- 6

p_dd <- paramtwo
p_dd$rx <- paramtwo$r2_1^2 * paramtwo$r2_2^2
p_dd$type <- "Epistatic"
p_dd$x <- 8


p <- rbind(p_r, p_a, p_d, p_aa, p_ad, p_dd)
p <- subset(p, !is.na(rx))
with(p, table(r1, r2, n, type, x))


p_summary <- ddply(p, .(x, r1, r2, n, type), summarise,
	nsim = length(rx),
	mean = mean(rx),
	sd   = sd(rx)
)

p_summary

ggplot(subset(p_summary, n > 100), aes(x = x, y = sd)) +
	geom_point() +
	geom_line(aes(colour = factor(n), linetype = factor(type))) +
	facet_grid(r1 ~ r2)

ggplot(p_summary, aes(x = x, y = mean)) +
	geom_point() +
	geom_line(aes(colour = factor(n), linetype = factor(type))) +
	facet_grid(r1 ~ r2)


