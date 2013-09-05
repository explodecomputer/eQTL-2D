# Simulate some intermediate frequency SNPs with known R
# Calculate empirical R
# Calculate sampling variance
# Calculate empirical R^8
# Calculate sampling variance of R^8


r = (x11 - p1q1) / sqrt(p1q1p2q2)
r = (x12 - p1q2) / sqrt(p1q1p2q2)


# given p1, q1, r we can calculate expected x1, x2, x3, x4



sampleSnp <- function(r, n, p1, p2)
{
	q1 <- 1 - p1
	q2 <- 1 - p2

	x11 <- r * sqrt(p1*q1*p2*q2) + p1*q1
	x12 <- p1*q2 - r * sqrt(p1*q1*p2*q2)
	x21 <- p2*q1 - r * sqrt(p1*q1*p2*q2)
	x22 <- r * sqrt(p1*q1*p2*q2) + p2*q2

	x <- sample(1:4, n, replace=T, prob=c(x11, x12, x21, x22))

	r_empirical <- (sum(x == 1)/n - p1*q1) / sqrt(p1*q1*p2*q2)
	return(r_empirical)
}


sampleSnp(0.8, 1200, 0.5, 0.5)
