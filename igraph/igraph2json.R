#-------------------------------------------------------------------------------------
library(RJSONIO) 
library(igraph)
#-------------------------------------------------------------------------------------
load("../analysis/data_lists.RData")

uniq_pairs_by_gene <- unique(meta[ ,c("snp1", "snp2", "probegene")])
rownames(uniq_pairs_by_gene) <- NULL

dir_all <- list(uniq_pairs_by_gene[ ,c("snp1", "snp2")], uniq_pairs_by_gene[ ,c("snp1", "probegene")], uniq_pairs_by_gene[ ,c("snp2", "probegene")])
dir_all <- lapply(dir_all, function(x) {
  colnames(x) <- c("From", "To")
  x
})

dir_all <- do.call(rbind, dir_all)
dir_all <- unique(dir_all)
rownames(dir_all) <- NULL
# data frame to igraph object
all_g <- graph.data.frame(dir_all, directed = TRUE)
# SNPs are black, genes are red
V(all_g)$color <- ifelse(grepl(x = V(all_g)$name, pattern = "^rs"), "black", "red")
# precompute ans store layout
my_lay <- layout.fruchterman.reingold(all_g)
rownames(my_lay) <- V(all_g)$name

all_together <- list(meta = meta, replicated = meta_replicated)

node_names_per_list <- lapply(all_together, function(x) {
    unique(c(x$snp1, x$snp2, x$probegene))
})

not_replicated <- setdiff(V(all_g)$name, node_names_per_list$replicated)
# Grey out not replicated nodes and edges from/to not replicated nodes
replicated <- all_g - not_replicated
E(all_g)$color <- "grey50" # Grey50 for all edges and Grey90 for not replicated
E(all_g)[from(not_replicated)]$color <- "grey90"
E(all_g)[to(not_replicated)]$color <- "grey90"

V(all_g)$color <- ifelse(V(all_g)$name %in% not_replicated, "grey90", V(all_g)$color)

#-------------------------------------------------------------------------------------
# Convert igraph object to JSON file to be loaded into D3.js
#-------------------------------------------------------------------------------------
# map node color(not planning on chaging edge color with D3) to number
# (will be treated as factor in D3) D3 will map back to color
lkup <- list("black" = 1, "red" = 2, "grey90" = 3) 
V(all_g)$group <- sapply(V(all_g)$color, function(x) lkup[[x]])
temp <- cbind(V(all_g)$name, V(all_g)$group)
colnames(temp) <- c("name", "group")
js1 <- toJSON(temp)

write.graph(all_g, "/tmp/edgelist.csv", format = "edgelist")
edges <- read.csv("/tmp/edgelist.csv", sep = " ", header = FALSE)
colnames(edges) <- c("source", "target")
edges <- as.matrix(edges)
js2 <- toJSON(edges)
asn <- paste('{"nodes":', js1,',"links":', js2, '}', sep = "")
write(asn, file = "../epi_web/epistasis.json")
#-------------------------------------------------------------------------------------