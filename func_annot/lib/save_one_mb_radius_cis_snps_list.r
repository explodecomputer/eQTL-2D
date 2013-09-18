# Decide which probes in cis and which probes in trans
# based on distance between TSS related to a expression probe and a SNP

# CAREFULL magic number

one_mb = 1000000 # bp

load("../analysis/data_lists.RData")

# Old annotation
annotateMy <- function(snp_df) {
    snp_df$pair1 <- ifelse(snp_df$chr1 != snp_df$probechr, "trans", "cis")
    snp_df$pair2 <- ifelse(snp_df$chr2 != snp_df$probechr, "trans", "cis")
    snp_df
}

meta_replicated <- annotateMy(meta_replicated)

# New annotation is a rule applied to the old-cis snps
# ie trans will stay trans and cis might become a trans only if 
# further than one_mb from corresponding TSS

a <- list(meta_replicated[ ,c("snp1", "probename", "pair1", "probegene")], meta_replicated[ ,c("snp2", "probename", "pair2", "probegene")])
a <- lapply(a, function(df) {
  colnames(df) <- c("rs_id", "probename", "relation_old", "probegene")
  df
})

a <- do.call(rbind, a)
a <- unique(a)
rownames(a) <- NULL



### Just to make sure we working on hg19, query online db 
# instead of using local files 
library(biomaRt)
snpmart <- useMart("snp", dataset = "hsapiens_snp")
epi_ids <- a$rs_id

snp_pos <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                    filters = c("snp_filter"),
                      value = epi_ids,
                       mart = snpmart)

#------------------------------------------------------------------------------
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

# drop transcript information keep per gene single TSS
a_naive <- merge(a_cis, exp_probe_ann[ ,1:4], by.x = "probename", by.y = "illumina_humanht_12", all.x = TRUE)
a_naive$chromosome_name_expr <- as.numeric(a_naive$chromosome_name_expr)


#------------------------------------------------------------------------------
dt_all <- read.table("../annot/all_in.ld", header = TRUE)
k562_gr <- import.bed("../chr_seg/k562.combined.bed", asRangedData = FALSE, genome = "hg19")
k562_gr <- keepSeqlevels(k562_gr, str_c("chr", 1:22))
hg19_seq_info <- seqinfo(Hsapiens)[str_c("chr", 1:22), ]


df2gr <- function(dfr) {
  GRanges(seqnames = str_c("chr", dfr$CHR_B),
          IRanges(start = dfr$BP_B, width = 1, name = dfr$SNP_B),
          strand = "*",
          r2 = dfr$R2,
          tag_snp = dfr$SNP_A,
          rs_id = dfr$SNP_B,
          seqinfo = hg19_seq_info)
}

all_gr <- df2gr(dt_all)

all_ol <- findOverlaps(all_gr, k562_gr)
# 
all_gr_f <- all_gr[queryHits(all_ol)]
# 
mcols(all_gr_f)$LOCATION <- mcols(k562_gr)$name[subjectHits(all_ol)]
#------------------------------------------------------------------------------


no_biomart_annot <- setdiff(a_cis$probename ,exp_probe_ann$illumina_humanht_12)
genes_no_annot <- meta_replicated[meta_replicated$probename %in% no_biomart_annot ,c("probegene")]
ann_via_gene_symbol <- getBM(attributes = c("chromosome_name", "start_position", "illumina_humanht_12",
                                               "hgnc_symbol", "transcript_start"),
         filters = "hgnc_symbol",
           value = c(genes_no_annot, "RAB44"),
            mart = ensembl_gene)

# "FLJ43093" -> "RAB44"


colnames(ann_via_gene_symbol)[c(1,2,5)] <- str_c(colnames(ann_via_gene_symbol)[c(1,2,5)], "_expr")

exp_probe_ann <- rbind(exp_probe_ann, ann_via_gene_symbol)


a_naive <- merge(a_cis, unique(exp_probe_ann[ ,1:4]), by.x = "probegene", by.y = "hgnc_symbol", all.x = TRUE)

a_naive$distance <- abs(a_naive$start_position_expr_expr - a_naive$chrom_start_snp_snp)

a_naive$relation_new <- ifelse(a_naive$distance > one_mb, "trans", "cis")






new_trans <- unique(subset(a_naive, relation_new == "trans", select = c("rs_id"), drop = TRUE))
a[a$rs_id %in% new_trans, ]$relation_old <- "trans"

a <- unique(a[ ,c(1,3)])

colnames(a)[2] <- "relation2probe"
write.csv(a, "DATA/relation2probe_one_mb_window.csv", quote = FALSE, row.names = FALSE)
#------------------------------------------------------------------------------

repids_rel <- read.csv("relation2probe_one_mb_window.csv")
dt <- read.table("replicated_ld_1mb.txt.ld", header = TRUE)

snp_gr <- df2gr(dt)

snp_gr_high_ld <- snp_gr[mcols(snp_gr)$r2 >= 0.8]
mcols(snp_gr_high_ld)$relation2probe <- repids_rel[match(mcols(snp_gr_high_ld)$tag_snp, repids_rel$rs_id), ]$relation2probe

cis_gr <- snp_gr_high_ld[mcols(snp_gr_high_ld)$relation2probe == "cis"]
trans_gr <- snp_gr_high_ld[mcols(snp_gr_high_ld)$relation2probe == "trans"]


all_ids <- unique(mcols(all_gr)$tag_snp)
cis_ids <- unique(mcols(cis_gr)$tag_snp)
trans_ids <- unique(mcols(trans_gr)$tag_snp)
#--------

fin <- unique(as.data.frame(mcols(all_gr_f)[ ,c("tag_snp", "LOCATION")]))

id_list <- list(all = all_ids, cis = cis_ids, trans = trans_ids)

res_count <- lapply(id_list, function(x) {
  dt <- ddply(fin[fin$tag_snp %in% x, ], .(LOCATION), nrow)
})

res_count <- lapply(names(res_count), function(x) {
  z <- res_count[[x]]
  colnames(z)[2] <- x
  z
})

count_df <- Reduce(function(x, y) merge(x, y, all = TRUE), res_count)
rownames(count_df) <- count_df$LOCATION

tag_snp_cnt <- sapply(id_list, length)

final_counts <- sweep(count_df[, 2:4], 2, tag_snp_cnt, "/")
final_counts$loc <- rownames(final_counts)

bin_test_res <- ddply(count_df, .(LOCATION), function(x) {
  cis_test   <- binom.test(x = x$cis, n = tag_snp_cnt[["cis"]], p = x$all / tag_snp_cnt[["all"]])$p.value
    trans_test <- binom.test(x = x$trans, n = tag_snp_cnt[["trans"]], p = x$all / tag_snp_cnt[["all"]])$p.value
    data.frame(c_pval = cis_test, t_pval = trans_test)
})

write.csv(bin_test_res, file = "binTestCombinedK562_one_mb_window.csv", quote = FALSE, row.names = FALSE)
write.csv(final_counts, file = "propCombinedK562_one_mb_window.csv", quote = FALSE, row.names = FALSE)













































