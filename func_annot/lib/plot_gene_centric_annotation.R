fin2 <- read.csv("gene_centric_overlap_with_index_snps_proportion.csv")
m_fin2 <- melt(fin2)
gfeatures_order <- c("intergenic", "intron", "promoter", "fiveUTR", "threeUTR", "coding", "spliceSite")
m_fin2 <- within(m_fin2, loc <- factor(loc, levels = gfeatures_order))

qplot(x = loc, y = value, fill = variable, data = m_fin2, geom = "bar", position = "dodge", stat = "identity") + theme_bw() + ylab("Proportion of index SNPs\n overlaping with a genomic feature") + xlab("Genomic Features") + scale_fill_brewer(palette = "Dark2", labels = c(per_tag_all = "Null distribution", per_tag_cis = "Cis SNPs", per_tag_trans = "Trans SNPs")) +
scale_x_discrete(labels = c(coding = "Coding", fiveUTR = "5\' UTR", intergenic = "Intergenic", intron = "Intron", promoter = "Promoter", spliceSite = "Splice site", threeUTR = "3\' UTR" )) + 
theme(
      # Legend
      legend.position="right",
      legend.title=element_blank(),
      legend.background=element_blank(),
      legend.key=element_blank(),
      # Text in general
      text=element_text(family="Helvetica", size=14), 
      # Strip aka facet lables
      strip.background=element_blank(),
      # inside of plot
      panel.grid=element_line(size=0),
      panel.border=element_rect(size=.6),
      panel.grid.minor.x=element_blank(),
      # Give me back my x axis
      axis.line.x=element_line(colour="black", size=2)) 
