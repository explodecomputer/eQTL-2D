# Hair ball plot
library(igraph)
load("../analysis/data_lists.RData")
uniq_pairs_by_gene <- unique(meta[ ,c("snp1", "snp2", "probegene")])
rownames(uniq_pairs_by_gene) <- NULL
head(uniq_pairs_by_gene)
dir_all <- list(uniq_pairs_by_gene[ ,c("snp1", "snp2")], uniq_pairs_by_gene[ ,c("snp1", "probegene")], uniq_pairs_by_gene[ ,c("snp2", "probegene")])
dir_all <- lapply(dir_all, function(x) {
	colnames(x) <- c("From", "To")
	x
})

dir_all <- do.call(rbind, dir_all)
dir_all <- unique(dir_all)
rownames(dir_all) <- NULL

all_g <- graph.data.frame(dir_all, directed = TRUE)

V(all_g)$color <- ifelse(grepl(x = V(all_g)$name, pattern = "^rs"), "black", "red")

my_lay <- layout.fruchterman.reingold(all_g)
rownames(my_lay) <- V(all_g)$name

all_together <- list(meta = meta, replicated = meta_replicated)

node_names_per_list <- lapply(all_together, function(x) {
    unique(c(x$snp1, x$snp2, x$probegene))
})

not_replicated <- setdiff(V(all_g)$name, node_names_per_list$replicated)

replicated <- all_g - not_replicated
E(all_g)$color <- "grey50"
E(all_g)[from(not_replicated)]$color <- "grey90"
E(all_g)[to(not_replicated)]$color <- "grey90"

V(all_g)$color <- ifelse(V(all_g)$name %in% not_replicated, "grey90", V(all_g)$color)

pdf("pale_gray_graph_of_interactions_2_lists_GrayEdges.pdf", width = 10, height = 10)
plot(all_g, vertex.size = 1.5, edge.width = 1, edge.arrow.size = 0,
     vertex.color = V(all_g)$color, vertex.label = NA,
     layout = my_lay[V(all_g)$name, ],
     vertex.frame.color = NA, edge.color = E(all_g)$color)
dev.off()
