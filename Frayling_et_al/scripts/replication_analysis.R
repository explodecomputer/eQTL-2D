# analysis of the egcut and ferhman replication datasets

# Read in data (one at a time)

# egcut
load("/Users/jpowell/repo/eQTL-2D/Frayling_et_al/replication_analysis/results/egcut/egcut-logtransformedquantilenormalized/replication_plus.RData")
egcut1 <- fit
load("/Users/jpowell/repo/eQTL-2D/Frayling_et_al/replication_analysis/results/egcut/egcut-40pcs/replication_plus.RData")
egcut2 <- fit

# ferhmann
load("/Users/jpowell/repo/eQTL-2D/Frayling_et_al/replication_analysis/results/Groningen/groningen-logtransformedquantilenormalized/replication_plus.RData")
feh1 <- fit

load("/Users/jpowell/repo/eQTL-2D/Frayling_et_al/replication_analysis/results/Groningen/groningen-40pcs/replication_plus.RData")
reh2 <- fit



combine <- function(x,y) return(pchisq(-2*log(x/2)-2*log(y/2),4,low=F))
i <- 26

egcut1[[i]]$snp1
egcut1[[i]]$snp2
egcut1[[i]]$probe
egcut1[[i]]$incsnp

egcut1[[i]]$replication_pnest
egcut2[[i]]$replication_pnest
feh1[[i]]$replication_pnest
feh2[[i]]$replication_pnest



out <- combine(a,b)
combine(out,c)





