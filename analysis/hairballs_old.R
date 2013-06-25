library(ggbio)
library(grid)
library(gridExtra)
library(plyr)
library(GenomicRanges)
library(qgraph)


# This is code to make hairball plots
# Redundant - Kostya did it better


key <- subset(sig, !duplicated(probename), select=c(probename, probegene))
nom <- subset(key, duplicated(probegene))$probename

sig2 <- subset(sig, ! probename %in% nom)

# More circles
a <- subset(sig2, select=c(snp1, snp2))
names(a) <- c("from", "to")
a$thickness <- 1
a$col <- as.character(sig2$rep)
a$col[a$col == 0] <- "white"
a$col[a$col == 1] <- "black"
a$col[a$col == 2] <- "black"
a$col2 <- a$col
a$col2[sig2$rep == 1] <- "white"
a$col3 <- as.character(sig2$rep)
a$col3[a$col3 == 0] <- ""
a$col3[a$col3 == 1] <- "black"
a$col3[a$col3 == 2] <- "red"
a$thickness2 <- a$thickness
a$thickness2[a$col3 == "red"] <- 1.5
a$cistrans <- "blue"
a$cistrans[sig2$c]


pdf(file="~/repo/eQTL-2D/analysis/images/hairballs_all.pdf")
qgraph(a, curve=0.1, borders=FALSE, esize=1, color="grey", edge.color=a$thickness, labels=FALSE, vsize=0.2, arrows=FALSE)
dev.off()
pdf(file="~/repo/eQTL-2D/analysis/images/hairballs_rep.pdf")
qgraph(a, curve=0.1, borders=FALSE, esize=1, color="grey", edge.color=a$col, labels=FALSE, vsize=0.2, arrows=FALSE)
dev.off()
pdf(file="~/repo/eQTL-2D/analysis/images/hairballs_rep2.pdf")
qgraph(a, curve=0.1, borders=FALSE, esize=1, color="grey", edge.color=a$col2, labels=FALSE, vsize=0.2, arrows=FALSE)
dev.off()
pdf(file="~/repo/eQTL-2D/analysis/images/hairballs_allrep.pdf")
qgraph(a, curve=0.1, borders=FALSE, esize=a$thickness2, color="grey", edge.color=a$col3, labels=FALSE, vsize=0.2, arrows=FALSE)
dev.off()

pdf(file="~/repo/eQTL-2D/analysis/images/hairballs_all_reps.pdf", width=18, height=6)
par(mfrow=c(1,3))
qgraph(a, curve=0.1, borders=FALSE, esize=1, color="red", edge.color=a$thickness, labels=FALSE, vsize=0.3, arrows=FALSE)
qgraph(a, curve=0.1, borders=FALSE, esize=1, color="red", edge.color=a$col, labels=FALSE, vsize=0.2, arrows=FALSE)
qgraph(a, curve=0.1, borders=FALSE, esize=1, color="red", edge.color=a$col2, labels=FALSE, vsize=0.2, arrows=FALSE)
dev.off()


# Are there SNP pairs that affect more than one probe?
m <- multipleSnps(sig)

# in all cases this is because the different probes tag the same gene
