

i <- which(info$Probe=="ILMN_2358626")
#i <- 21
info[i,]
index <- which(meta$probename==as.character(info$Probe[i]) & meta$snp1==as.character(info$SNP1[i]) & meta$snp2==as.character(info$SNP2[i]))
p8 <- meta$pfull[index]
var_8_explain.fun(p8)


index <- which(meta$probename==as.character(info$Probe[i]) & meta$snp1==as.character(info$SNP1[i]) & meta$snp2==as.character(info$SNP2[i]))
p4 <- meta$pnest_fehr[index]
p8 <- meta$pfull_fehr[index]

var_4_explain.fun(p4)
var_8_explain.fun(p8)


####
index <- which(meta$probename==as.character(info$Probe[i]) & meta$snp1==as.character(info$SNP1[i]) & meta$snp2==as.character(info$SNP2[i]))
p4 <- meta$pnest_egcut[index]
p8 <- meta$pfull_egcut[index]

var_4_explain_fehr.fun(p4)
var_8_explain_fehr.fun(p8)






var_8_explain.fun <- function(p) {
	p <- 10^-p
	ch <- qchisq(p, df=8, lower.tail=FALSE)
	v <- (ch/846)/(1+(ch/846))*100
	print(v)

}

var_4_explain.fun <- function(p) {
	p <- 10^-p
	ch <- qchisq(p, df=4, lower.tail=FALSE)
	v <- (ch/846)/(1+(ch/846))*100
	print(v)

}


var_8_explain_fehr.fun <- function(p) {
	p <- 10^-p
	ch <- qchisq(p, df=8, lower.tail=FALSE)
	v <- (ch/1250)/(1+(ch/1250))*100
	print(v)

}

var_4_explain_fehr.fun <- function(p) {
	p <- 10^-p
	ch <- qchisq(p, df=4, lower.tail=FALSE)
	v <- (ch/1250)/(1+(ch/1250))*100
	print(v)

}



var_1_explain.fun <- function(p) {
	#p <- 10^-p
	ch <- qchisq(p, df=1, lower.tail=FALSE)
	v <- (ch/846)/(1+(ch/846))*100
	return(v)

}

inc$R2 <- var_1_explain.fun(inc$PVALUE)

