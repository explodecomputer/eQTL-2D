# re-annotate snp list for Trynka

library(biomaRt)
#-------------------------------------------------------------------------------
snpmart <- useMart("snp", dataset = "hsapiens_snp")
load("DATA/sig_with_one_mb_window_annotation.RData")

all_snps <- function(x) {
    unique(c(x$snp1, x$snp2))
}

trans_snps <- function(x) {
    unique(c(subset(x, pair1 == "trans")$snp1, subset(x, pair2 == "trans")$snp2))
}

cis_snps <- function(x) {
    unique(c(subset(x, pair1 == "cis")$snp1, subset(x, pair2 == "cis")$snp2))
}

get_snp_pos <- function(rs_ids) {
  z <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                filters = c("snp_filter"),
                  value = rs_ids,
                   mart = snpmart)
  z$chr_name <- str_c("chr", z$chr_name)
  # remove mapping to non placed chromosomes
  autosomes <- str_c("chr", 1:22)
  z <- subset(z, chr_name %in% autosomes)
  z
}

save_mappings <- function(snp_pos, myname) {
    write.table(snp_pos, str_c("DATA/one_mb_rule_final_h3k4me3_lists/", myname, ".snpmappings.txt"), 
                         col.names = FALSE, 
                         row.names = FALSE, 
                             quote = FALSE, 
                               sep = "\t")
}

per_df <- function(set_name) {
    df <- all_together[[set_name]]
    df_all <- all_snps(df)
    df_cis <- cis_snps(df)
    df_trans <- trans_snps(df)
    all_mappings <- get_snp_pos(df_all)
    cis_mappings <- all_mappings[all_mappings$refsnp_id %in% df_cis, ]
    trans_mappings <- all_mappings[all_mappings$refsnp_id %in% df_trans, ]
    save_mappings(all_mappings, str_c(set_name, "_all"))
    save_mappings(cis_mappings, str_c(set_name, "_cis"))
    save_mappings(trans_mappings, str_c(set_name, "_trans"))
}

all_together <- list(levis = sig)
# remove the snp that are not in the 1KG (see below)
all_together <- lapply(all_together, function(x) {
    x <- subset(x, snp1 != "rs7405659" & snp2 != "rs7405659")
})

all_together <- lapply(all_together, function(x) {
    x <- subset(x, chr1 <= 22 & chr2 <= 22)
})

lapply(names(all_together), per_df)
