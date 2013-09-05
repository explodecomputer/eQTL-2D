load("../analysis/data_lists.RData")
# R>lsos()
#                            Type    Size PrettySize Rows Columns
# sig                  data.frame 4446928     4.2 Mb  501      35
# meta                 data.frame 4179800       4 Mb  434      44
# meta_replicated      data.frame 3327288     3.2 Mb  345      44
# chr_interaction_list data.frame   11768    11.5 Kb   44      10

chr_int_ids <- unique(chr_interaction_list$snp1, chr_interaction_list$snp2)
write.table(chr_int_ids, file = "DATA/chr_int_ids.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

annotateMy <- function(snp_df) {
    snp_df$pair1 <- ifelse(snp_df$chr1 != snp_df$probechr, "trans", "cis")
    snp_df$pair2 <- ifelse(snp_df$chr2 != snp_df$probechr, "trans", "cis")
    snp_df
}

x <- meta_replicated
x <- subset(x, chr1 <= 22 & chr2 <= 22)
x <- subset(x, snp1 != "rs7405659" & snp2 != "rs7405659")
x <- annotateMy(x)
head(x)
z <- list(x[ ,c("snp1", "pair1")], x[ ,c("snp2", "pair2")])
z <- lapply(z, function(x) {
  colnames(x) <- c("snp", "regulates_in")
  x
})

z <- do.call(rbind, z)
z <- unique(z)
rownames(z) <- NULL
trans_ids <- subset(z, regulates_in == "trans")
cis_ids <- subset(z, regulates_in == "cis")
write.table(trans_ids[ ,1], file = "DATA/trans_ids.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(cis_ids[ ,1], file = "DATA/cis_ids.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)