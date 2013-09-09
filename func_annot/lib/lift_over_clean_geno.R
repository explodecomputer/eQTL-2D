library(rtracklayer)
#----------------------------------------------------------------------------------------
all_ids <- read.table(pipe("cut -f1,2,4 clean_geno_final.bim"), header = FALSE, sep = "\t")
colnames(all_ids) <- c("CHR", "SNP", "BP")
all_ids <- subset(all_ids, CHR %in% 1:22)

df2gr <- function(df) {
  df.gr <- GRanges(seqnames = Rle(str_c("chr", df$CHR)),
                     ranges = IRanges(df$BP, width = 1, names = df$SNP),
                     strand = "*",
                    seqinfo = SeqinfoForUCSCGenome("hg18"),
                     rs_ids = df$SNP)
  chain <- import.chain("/hox/u/uqkshakb/seq/ucsc_chain/hg18ToHg19.over.chain")
  df.gr.hg19 <- liftOver(df.gr, chain)
  df.gr.hg19 <- unlist(df.gr.hg19)
  full <- SeqinfoForUCSCGenome("hg19")
  seqlevels(df.gr.hg19) <- seqlevels(full)
  seqinfo(df.gr.hg19) <- full
  df.gr.hg19 <- keepSeqlevels(df.gr.hg19, str_c("chr", c(as.character(1:22))))
  names(mcols(df.gr.hg19)) <- "rs_ids"
  df.gr.hg19
}

all_ids_gr.hg19 <- df2gr(all_ids)
save(all_ids_gr.hg19, file = "Clean_geno_final_GRangesLifted2hg19.RData")
write.table(unique(mcols(all_ids_gr.hg19)$rs_ids), file = "clean_geno_ids.txt")