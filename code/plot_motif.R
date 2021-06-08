# plot mini motifs from graph6

#setwd("FNM/code/")
library(igraph)


gmat <- as.matrix(read.table("tmp.inp"))
g <- graph_from_adjacency_matrix(gmat, 'undirected', diag=F)


svg(file = "motif.svg", height=4, width=4)
plot(g, vertex.color="black", vertex.label=NA, vertex.size=45, edge.width=7, edge.color="black")
dev.off()

q(save='no')
