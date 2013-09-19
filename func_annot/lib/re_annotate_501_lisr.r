# 501 list is in sig dataframe
load("../analysis/data_lists.RData")
one_mb = 1000000 # bp
annotateMy <- function(snp_df) {
    snp_df$pair1 <- ifelse(snp_df$chr1 != snp_df$probechr, "trans", "cis")
    snp_df$pair2 <- ifelse(snp_df$chr2 != snp_df$probechr, "trans", "cis")
    snp_df
}

sig <- annotateMy(sig)

a <- list(sig[ ,c("snp1", "probename", "pair1", "probegene")], sig[ ,c("snp2", "probename", "pair2", "probegene")])
a <- lapply(a, function(df) {
  colnames(df) <- c("rs_id", "probename", "relation_old", "probegene")
  df
})

a <- do.call(rbind, a)
a <- unique(a)
rownames(a) <- NULL


library(biomaRt)
snpmart <- useMart("snp", dataset = "hsapiens_snp")
epi_ids <- a$rs_id

snp_pos <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                    filters = c("snp_filter"),
                      value = epi_ids,
                       mart = snpmart)

ensembl_gene <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

exp_probe_ann <- getBM(attributes = c("chromosome_name", "start_position", "illumina_humanht_12",
                                               "hgnc_symbol", "transcript_start"),
         filters = "illumina_humanht_12",
           value = unique(a$probename),
            mart = ensembl_gene)

a_cis <- subset(a, relation_old == "cis")

colnames(snp_pos)[2:3] <- str_c(colnames(snp_pos)[2:3], "_snp")
colnames(exp_probe_ann)[c(1,2,5)] <- str_c(colnames(exp_probe_ann)[c(1,2,5)], "_expr")

a_cis <- merge(a_cis, snp_pos, by.x = "rs_id", by.y = "refsnp_id", all.x = TRUE)

# # drop transcript information keep per gene single TSS
# a_naive <- merge(a_cis, exp_probe_ann[ ,1:4], by.x = "probename", by.y = "illumina_humanht_12", all.x = TRUE)
# a_naive$chromosome_name_expr <- as.numeric(a_naive$chromosome_name_expr)


no_biomart_annot <- setdiff(a_cis$probename ,exp_probe_ann$illumina_humanht_12)
genes_no_annot <- sig[sig$probename %in% no_biomart_annot ,c("probegene")]
genes_no_annot <- unique(genes_no_annot)

ann_via_gene_symbol <- getBM(attributes = c("chromosome_name", "start_position", "illumina_humanht_12",
                                               "hgnc_symbol", "transcript_start"),
         filters = "hgnc_symbol",
           value = c(genes_no_annot, "RAB44", "LINC00339"),
            mart = ensembl_gene)

# "FLJ43093" -> "RAB44"

setdiff(genes_no_annot, ann_via_gene_symbol$hgnc_symbol)


colnames(ann_via_gene_symbol)[c(1,2,5)] <- str_c(colnames(ann_via_gene_symbol)[c(1,2,5)], "_expr")

exp_probe_ann <- rbind(exp_probe_ann, ann_via_gene_symbol)
exp_probe_ann[exp_probe_ann$hgnc_symbol == "C1orf86", ]$hgnc_symbol <- "C1ORF86"


a_naive <- merge(a_cis, unique(exp_probe_ann[ ,1:4]), by.x = "probegene", by.y = "hgnc_symbol", all.x = TRUE)

a_naive$distance <- abs(a_naive$start_position_expr - a_naive$chrom_start_snp)

a_naive$relation_new <- ifelse(a_naive$distance > one_mb, "trans", "cis")






new_trans <- unique(subset(a_naive, relation_new == "trans", select = c("rs_id"), drop = TRUE))
a[a$rs_id %in% new_trans, ]$relation_old <- "trans"

a <- unique(a[ ,c(1,3)])

colnames(a)[2] <- "relation2probe"
write.csv(a, "DATA/all_aka_from_sig_relation2probe_one_mb_window.csv", quote = FALSE, row.names = FALSE)



