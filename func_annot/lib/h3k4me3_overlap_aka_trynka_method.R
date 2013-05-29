#-------------------------------------------------------------------------------
# Trynka et al analysis
# blame k.shakhbazov@uq.edu.au
#-------------------------------------------------------------------------------
# email from Gib:
# eQTL-2D/analysis/interaction_list_replication_summary.RData
# Then this has 4 objects.
# sig_all = all SNPs including control SNPs that should be ignored (filter == 3)
# sig = SNPs that were significant in our dataset
# sig_rep1 = SNPs that are replicated at least once in an independent dataset
# sig_rep2 = SNPs that are replicated in both independent datasets
# So ignore sig_all
#-------------------------------------------------------------------------------
library(biomaRt)
#-------------------------------------------------------------------------------
snpmart <- useMart("snp", dataset = "hsapiens_snp")
load("../analysis/interaction_list_replication_summary.RData")
#-------------------------------------------------------------------------------
annotateMy <- function(snp_df) {
    snp_df$pair1 <- ifelse(snp_df$chr1 != snp_df$probechr, "trans", "cis")
    snp_df$pair2 <- ifelse(snp_df$chr2 != snp_df$probechr, "trans", "cis")
    snp_df
}

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
  z
}

save_mappings <- function(snp_pos, myname) {
    write.table(snp_pos, str_c("DATA/", myname, ".snpmappings.txt"), 
                         col.names = FALSE, 
                         row.names = FALSE, 
                             quote = FALSE, 
                               sep = "\t")
}

per_df <- function(set_name) {
    df <- all_together[[set_name]]
    df <- annotateMy(df)
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
#-------------------------------------------------------------------------------
all_together = list(original = sig, once = sig_rep1, twice = sig_rep2)
# remove the snp that not in the 1KG (see below)
all_together <- lapply(all_together, function(x) {
    x <- subset(x, snp1 != "rs7405659" & snp2 != "rs7405659")
})

all_together <- lapply(all_together, function(x) {
    x <- subset(x, chr1 <= 22 & chr2 <= 22)
})


lapply(names(all_together), per_df)
#-------------------------------------------------------------------------------
# in the DATA folder 
grep NA *.snpmappings.txt

# once_all.snpmappings.txt:rs7405659  chrNA   NA
# once_cis.snpmappings.txt:rs7405659  chrNA   NA
# original_all.snpmappings.txt:rs7405659  chr NA
# original_cis.snpmappings.txt:rs7405659  chr NA

# rs7405659 is not in the 1KG
#-------------------------------------------------------------------------------





























