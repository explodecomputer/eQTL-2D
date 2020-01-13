library(data.table)
library(dplyr)
library(tidyr)

h2014 <- fread("../../investigation/data/total_analysis_data.txt")

# need a list of genes
# their cis snp
# their finemap snp

eqtlgen <- fread("https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")

bim <- fread("../data/combined.bim")
disc <- fread("../data/disc.bim")

tab <- rbind(
	subset(h2014, chr1 == probechr) %>%
		dplyr::select(gene=gene, cissnp=snp1),
	subset(h2014, chr2 == probechr) %>%
		dplyr::select(gene=gene, cissnp=snp2)
) %>% filter(!duplicated(paste(gene, cissnp)))
tab <- subset(tab, cissnp %in% disc$V2)

genelist <- unique(tab$gene)

eqtl <- subset(eqtlgen, GeneSymbol %in% genelist)
eqtl <- subset(eqtl, SNP %in% bim$V2)
eqtl$rsq <- eqtl$Zscore^2 / (eqtl$NrSamples + eqtl$Zscore^2)
eqtl <- arrange(eqtl, desc(rsq))
finemap <- eqtl %>% filter(!duplicated(GeneSymbol)) %>% dplyr::select(gene = GeneSymbol, finemap=SNP, rsq=rsq) %>% merge(., tab)

write.table(finemap, file="../data/finemap.txt", row=F, col=F, qu=F)
