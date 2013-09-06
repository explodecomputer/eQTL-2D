dt <- read.csv("propCombinedK562.csv")

segmenet_order <- c("E", "CTCF", "WE", "PF", "TSS", "T", "R")

m_fin2 <- melt(dt)
m_fin2 <- within(m_fin2, loc <- factor(loc, levels = segmenet_order))

qplot(x = loc, y = value, fill = variable, data = m_fin2, geom = "bar", position = "dodge", stat = "identity") + theme_bw()
+ ylab("Proportion of index SNPs\n overlaping with a chromatin segment") + xlab("Chromatin states")
+ scale_fill_brewer(palette = "Dark2", labels = c(all = "Null distribution", cis = "Cis SNPs", trans = "Trans SNPs")) +
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

ggsave("~/Desktop/ChromatinStatesOverlapFixedOrder.pdf")
