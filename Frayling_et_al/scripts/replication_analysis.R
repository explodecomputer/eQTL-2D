# analysis of the egcut and ferhman replication datasets

# Read in data (one at a time)

# egcut
load("/Users/jpowell/repo/eQTL-2D/Frayling_et_al/replication_analysis/results/egcut/egcut-logtransformedquantilenormalized/replication_plus.RData")

load("/Users/jpowell/repo/eQTL-2D/Frayling_et_al/replication_analysis/results/egcut/egcut-40pcs/replication_plus.RData")


# ferhmann
load("/Users/jpowell/repo/eQTL-2D/Frayling_et_al/replication_analysis/results/Groningen/groningen-logtransformedquantilenormalized/replication_plus.RData")

load("/Users/jpowell/repo/eQTL-2D/Frayling_et_al/replication_analysis/results/Groningen/groningen-40pcs/replication_plus.RData")





for(i in 1:26) {

print(fit[[i]]$replication_pnest)

}