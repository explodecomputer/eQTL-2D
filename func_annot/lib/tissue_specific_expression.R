#---------------------------------------------------
# blame: k.shakhbazov@uq.edu.au
# Do genes regulated via epi snps have particular tissue specificity properties???
#---------------------------------------------------

# expression data for human tissues comes from biogps
# wget http://plugins.biogps.org/download/gnf1h-gcrma.zip
# wget http://plugins.biogps.org/download/human_sample_annot.csv
# wget http://plugins.biogps.org/download/gnf1h-anntable.zip
# lives in the DATA folder but not commited via .gitignore

library(biomaRt)

dt <- read.csv("../../replication/results/matched_replication_results.csv", header = TRUE)
exprs <- read.csv(pipe("gunzip -c ../DATA/gnf1h-gcrma.zip"), header = TRUE)
annot <- read.delim(pipe("gunzip -c ../DATA/gnf1h-anntable.zip"), header = TRUE)

#---------------------------------------------------
epi_genes2probes <- subset(annot, Symbol %in% unique(dt$probegene))
# epi_exprs <- subset(exprs, X %in% )


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
probes <- getBM(attributes = c("ensembl_gene_id", "affy_hg_u133a", "hgnc_symbol", "illumina_humanht_12"),
                   filters = "illumina_humanht_12",
                    values = unique(dt$probename),
                      mart = ensembl)

#---------------------------------------------------

merge(annot, exprs, )
