#------------------------------------------------------------------------------
library(GenomicRanges)
library(VariantAnnotation)
library(GenomicFeatures)
#------------------------------------------------------------------------------
dt <- read.table("all_in.ld", header = TRUE)
repids_rel <- read.csv("replicated_ids_relation2probe.txt")
snp_gr <- GRanges(seqnames = str_c("chr", dt$CHR_B),
                  IRanges(start = dt$BP_B, width = 1, name = dt$SNP_B),
                  strand = "*",
                  r2 = dt$R2,
                  tag_snp = dt$SNP_A,
                  seqinfo = my_seq_info)


gc16 <- loadDb("Gencode16.sqlite")
all_var <- locateVariants(snp_gr, gc16, AllVariants())

mcols(all_var)$query_rs_id <- names(snp_gr)[mcols(all_var)$QUERYID]
mcols(all_var)$tag_snp <- mcols(snp_gr)$tag_snp[mcols(all_var)$QUERYID]

#------------------------------------------------------------------------------
# save(all_var, file = "AllGenoAndInLD_OverlapAnnotation.RData")
# load("AllGenoAndInLD_OverlapAnnotation.RData")
#------------------------------------------------------------------------------

avar <- as.data.frame(mcols(all_var)[ ,c("LOCATION", "tag_snp")])
avar <- unique(avar)
cis <- subset(avar, tag_snp %in% subset(repids_rel, relation2probe == "cis", select = "rs_id", drop = TRUE))
trans <- subset(avar, tag_snp %in% subset(repids_rel, relation2probe == "trans", select = "rs_id", drop = TRUE))

per_tag_all <- table(avar$LOCATION)
per_tag_cis <- table(cis$LOCATION)
per_tag_trans <- table(trans$LOCATION)

final_counts <- as.data.frame(t(rbind(per_tag_all, per_tag_cis, per_tag_trans)))

nm_cnt <- c(table(repids_rel$relation2probe), all = length(unique(avar$tag_snp)))[c(3,1,2)]
fin2 <- sweep(final_counts[ ,1:3], 2, nm_cnt, `/`)
fin2$loc <- rownames(fin2)

write.csv(fin2, file = "DATA/gene_centric_overlap_with_index_snps_proportion.csv", quote = FALSE, row.names = FALSE)