#########################################################################
# TSS Predicted promoter region including TSS
# PF  Predicted promoter flanking region
# E Predicted enhancer
# WE  Predicted weak enhancer or open chromatin cis regulatory element
# CTCF  CTCF enriched element
# T Predicted transcribed region
# R Predicted Repressed or Low Activity region
#########################################################################

#--------------------------------------------------------------------------------
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
#--------------------------------------------------------------------------------
# GF annotation
k562_gr <- import.bed("../DATA/k562.combined.bed", asRangedData = FALSE, genome = "hg19")
k562_gr <- keepSeqlevels(k562_gr, str_c("chr", 1:22))

# SNP annotation
repids_rel <- read.csv("../DATA/replicated_ids_relation2probe.txt")
hg19_seq_info <- seqinfo(Hsapiens)[str_c("chr", 1:22), ]
dt <- read.table("../DATA/replicated_ld_1mb.txt.ld", header = TRUE)
# data.frame to GenomicRanges
df2gr <- function(dfr) {
  GRanges(seqnames = str_c("chr", dfr$CHR_B),
          IRanges(start = dfr$BP_B, width = 1, name = dfr$SNP_B),
          strand = "*",
          r2 = dfr$R2,
          tag_snp = dfr$SNP_A,
          rs_id = dfr$SNP_B,
          seqinfo = hg19_seq_info)
}

snp_gr <- df2gr(dt)

snp_gr_high_ld <- snp_gr[mcols(snp_gr)$r2 >= 0.8]
mcols(snp_gr_high_ld)$relation2probe <- repids_rel[match(mcols(snp_gr_high_ld)$tag_snp, repids_rel$rs_id), ]$relation2probe

cis_gr <- snp_gr_high_ld[mcols(snp_gr_high_ld)$relation2probe == "cis"]
trans_gr <- snp_gr_high_ld[mcols(snp_gr_high_ld)$relation2probe == "trans"]
# all genotyped snps passed QC (null-distribution)
dt_all <- read.table("../annot/all_in.ld", header = TRUE)
all_gr <- df2gr(dt_all)

all_ids <- unique(mcols(all_gr)$tag_snp)
cis_ids <- unique(mcols(cis_gr)$tag_snp)
trans_ids <- unique(mcols(trans_gr)$tag_snp)

#--------------------------------------------------------------------------------
all_ol <- findOverlaps(all_gr, k562_gr)
all_gr_f <- all_gr[queryHits(all_ol)]
mcols(all_gr_f)$LOCATION <- mcols(k562_gr)$name[subjectHits(all_ol)]

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

write.csv(bin_test_res, file = "binTestCombinedK562.csv", quote = FALSE, row.names = FALSE)
write.csv(final_counts, file = "propCombinedK562.csv", quote = FALSE, row.names = FALSE)
#--------------------------------------------------------------------------------