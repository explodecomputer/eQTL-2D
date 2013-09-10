load("../analysis/data_lists.RData")

annotateMy <- function(snp_df) {
    snp_df$pair1 <- ifelse(snp_df$chr1 != snp_df$probechr, "trans", "cis")
    snp_df$pair2 <- ifelse(snp_df$chr2 != snp_df$probechr, "trans", "cis")
    snp_df
}

meta_replicated <- annotateMy(meta_replicated)

a <- list(meta_replicated[ ,c("snp1", "pair1")], meta_replicated[ ,c("snp2", "pair2")])
a <- lapply(a, function(df) {
  colnames(df) <- c("rs_id", "in")
  df
})

a <- do.call(rbind, a)
a <- unique(a)
rownames(a) <- NULL
colnames(a)[2] <- "relation2probe"
write.csv(a, file = "DATA/replicated_ids_relation2probe.txt", quote = FALSE, row.names = FALSE)
