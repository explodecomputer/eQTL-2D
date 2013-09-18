load("../analysis/data_lists.RData")
repids_rel <- read.csv("DATA/relation2probe_one_mb_window.csv")

meta_replicated$pair1 <- sapply(meta_replicated$snp1, function(x) repids_rel[repids_rel$rs_id == x, ]$relation2probe)
meta_replicated$pair2 <- sapply(meta_replicated$snp2, function(x) repids_rel[repids_rel$rs_id == x, ]$relation2probe)

table(meta_replicated[ ,c("pair1", "pair2")])
#        pair2
# pair1   cis trans
#   cis    23   125
#   trans 187    10