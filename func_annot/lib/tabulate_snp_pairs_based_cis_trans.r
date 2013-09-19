load("../analysis/data_lists.RData")
repids_rel <- read.csv("DATA/relation2probe_one_mb_window.csv")

meta_replicated$pair1 <- sapply(meta_replicated$snp1, function(x) repids_rel[repids_rel$rs_id == x, ]$relation2probe)
meta_replicated$pair2 <- sapply(meta_replicated$snp2, function(x) repids_rel[repids_rel$rs_id == x, ]$relation2probe)

table(meta_replicated[ ,c("pair1", "pair2")])
#        pair2
# pair1   cis trans
#   cis    23   125
#   trans 187    10


#--------------------------------------
# NOW FOR ALL SNP PAIRS AKA 501 PAIRS LIST
#--------------------------------------

repids_rel_all <- read.csv("DATA/all_aka_from_sig_relation2probe_one_mb_window.csv")
# drop SNP which both cis and trans
repids_rel_all <- repids_rel_all[-which(repids_rel_all$rs_id == "rs10120023"), ]

sig$pair1 <- repids_rel_all[match(sig$snp1, repids_rel_all$rs_id), "relation2probe"]
sig$pair2 <- repids_rel_all[match(sig$snp2, repids_rel_all$rs_id), "relation2probe"]

# R>subset(repids_rel_all, rs_id == "rs10120023")
#          rs_id relation2probe
# 3   rs10120023          trans
# 730 rs10120023            cis

# Done based on difference between chromosomes names between the SNP and a probe, cis one
# was manually queried at Biomart see below
sig[sig$snp1 == "rs10120023", ]$pair1 <- "trans"
sig[sig$snp2 == "rs10120023", ]$pair2 <- c("cis", "trans")
save(sig, file = "DATA/sig_with_one_mb_window_annotation.RData")
# R>table(sig[ ,c("pair1", "pair2")])
#        pair2
# pair1   cis trans
#   cis    26   185
#   trans 277    13
# R>sum(table(sig[ ,c("pair1", "pair2")]))
# [1] 501

#---------------------------------------------------------------------
library(biomaRt)
ensembl_gene <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
getBM(attributes = c("chromosome_name", "start_position", "illumina_humanht_12",
                                               "hgnc_symbol", "transcript_start"),
         filters = "hgnc_symbol",
           value = c("FCN1"),
            mart = ensembl_gene)

#   chromosome_name start_position illumina_humanht_12 hgnc_symbol
# 1               9      137801431        ILMN_1668063        FCN1
#   transcript_start
# 1        137801431

snpmart <- useMart("snp", dataset = "hsapiens_snp")
getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                    filters = c("snp_filter"),
                      value = "rs10120023",
                       mart = snpmart)
#    refsnp_id chr_name chrom_start
# 1 rs10120023        9   137810259

