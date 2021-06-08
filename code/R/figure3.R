# Figure 3 FNM
# SP (sebastian@pechmannlab.net)
# 2021


library(ggplot2)
library(reshape2)
library(cowplot)
library(igraph)
library(ggraph)

setwd("~/FNM/")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 3A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gi <- as.data.frame(read.table("data/figures/figure2/gi_summary.txt", header=T))


svg(file = "figures/Figure3/A_gihisto.svg", height = 5, width = 4)

plot_pos <- ggplot(gi, aes(pct_pos) ) + 
  geom_histogram(binwidth=0.05, position='dodge', fill="gold", color="white") + 
  labs(x="% positive GI", y="Count") +
  theme_classic()+
  theme(
    text = element_text(size=18) 
  )

plot_neg <- ggplot(gi, aes(pct_neg) ) + 
  geom_histogram(binwidth=0.05, position='dodge', fill="darkblue", color="white") + 
  labs(x="% negative GI", y="Count") +
  theme_classic()+
  theme(
    text = element_text(size=18) 
  )

plot_grid(plot_pos, plot_neg, labels ="", ncol = 1, align = 'v')

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 3B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

complex_gi_pos  <- as.matrix(read.table("data/figures/figure3/complex_network_gi_pos.txt"))
complex_gi_neg  <- as.matrix(read.table("data/figures/figure3/complex_network_gi_neg.txt"))

res <- c()
for (i in 1:nrow(complex_gi_pos)){
  for (j in i:ncol(complex_gi_pos)){
    if (i < j){
      current_pos <- complex_gi_pos[i,j]
      current_neg <- complex_gi_neg[i,j]
      if (current_pos > 0 || current_neg > 0){
        current_score <- current_pos/(current_pos+current_neg)
        res <- c(res, current_score)
      }
    }
  }
}
complex_pct_gi <- data.frame(pct=res)


svg(file = "figures/Figure3/B_pctgi.svg", height = 2.5, width = 4)

ggplot(complex_pct_gi, aes(x=pct)) + 
  geom_histogram(binwidth = 0.1, color="white") + 
  labs(x="", y="Count") +
  scale_x_continuous(breaks=c(0, 0.5, 1), labels=c("100% neg", "50%", "100% pos")) +
  theme_classic() +
  theme(
    text = element_text(size=18)
  )

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 3C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# the code follows the edge-bundling code example on R-Gallery

hier <- as.data.frame(read.table("data/figures/figure3/network_hierarchy.txt", header=T))
conn_ppi <- as.data.frame(read.table("data/figures/figure3/network_connection.txt", header=T))
conn_gi <- as.data.frame(read.table("data/figures/figure3/network_connection_gi.txt", header=T))
vertices <- as.data.frame(read.table("data/figures/figure3/network_vertex.txt", header=T))

vertices$id <- NA
myleaves <- which(is.na( match(vertices$name, edges$from) ))
myleaves2 <-  match(vertices$name, hier$to) 
nleaves <- length(myleaves)
vertices$id <- myleaves2
vertices$angle <- 90 - 360 * vertices$id / nleaves
vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)
vertices$angle <- ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

mygraph <- graph_from_data_frame(hier, vertices=vertices)

from_ppi <- match( conn$from, vertices$name )
to_ppi   <- match( conn$to,   vertices$name)
weight_ppi <- rep( log(abs(conn$weight))/4, each=100)

from_gi <- match( conn_gi$from, vertices$name )
to_gi   <- match( conn_gi$to,   vertices$name)
weight_gi <- rep( abs(conn_gi$weight)*5, each=100)
color_gi <- rep( ifelse( conn_gi$weight < 0, "darkblue", "gold"), each=100)




svg(file = "figures/Figure3/C_full6_gi_names.svg", height = 5, width = 5)
ggraph(mygraph, layout='dendrogram', circular=T) + 
  geom_conn_bundle(data=get_con(from=from_ppi, to=to_ppi), edge_width=weight_ppi, alpha=0.6, color="#222222", tension=0) + 
  geom_conn_bundle(data=get_con(from=from_gi, to=to_gi), edge_width=weight_gi, alpha=0.2, color=color_gi, tension=0.9) + 
  geom_node_text(aes(x=x*1.15, y=y*1.15, filter=leaf, angle=angle, hjust=hjust, label=short), size=3.3, alpha=1) + 
  geom_node_point(aes(filter=leaf, x=x*1.05, y=y*1.05, size=vert$degree)) + 
  theme_void() +
  theme(
    legend.position='none',
    plot.margin=unit(c(0,0,0,0), 'cm')
  )
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 3D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

complex_ppi <- as.matrix(read.table("data/figures/figure3/complex_network_ppi.txt"))
complex_gi  <- as.matrix(read.table("data/figures/figure3/complex_network_gi.txt"))

res <- data.frame(ppi=c(), gi=c())
for (i in 1:nrow(complex_ppi)){
  for (j in i:ncol(complex_ppi)){
    if (i < j){
      current_ppi <- complex_ppi[i,j]
      current_gi <- complex_gi[i,j]
      if (current_ppi > 0 || current_gi > 0){
        res <- rbind(res, c(current_ppi, current_gi))
      }
    }
  }
}
colnames(res) <- c("PPI", "GI")


complex_pct <- data.frame(type=c("PPI+GI", "PPI", "GI"),
                          value=c( sum( res[res$PPI > 0,2] != 0 )/nrow(res), 
                                   sum( res[res$PPI > 0,2] == 0 )/nrow(res), 
                                   sum( res[res$GI > 0,1] == 0 )/nrow(res) ) )

svg(file = "figures/Figure3/D_pctint.svg", height = 2.5, width = 4)

ggplot(complex_pct, aes(x=type, y=value)) + 
  geom_col(aes(fill=type)) +
  scale_fill_manual(values=c("#BBBBBB", "#777777", "#222222")) + 
  labs(y="Fraction of interactions", x="") + 
  coord_flip() + 
  theme_classic() + 
  theme(
    text = element_text(size=18),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(), 
    legend.position = 'none'
  )

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 3E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

names <- read.table("data/figures/figure3/bestmotif_orfs.txt", header=T)
names$complex <- paste(names$name, ' (', names$short, ')', sep='')

df_net <- as.data.frame(read.table("data/figures/figure3/bestmotif_df.txt", header=T))
graph <- graph_from_data_frame(df_net, vertices=names$ORF)

df_net$alpha <- ifelse(df_net$type=="PPI", 1, 0.95)

svg("figures/Figure3/E_motif.svg", height=5, width=4)

ggraph(graph, 'stress') + 
  geom_edge_arc(aes(colour=df_net$type, alpha=df_net$alpha), width=1.5, strength=df_net$strength) + 
  geom_node_point(size=8 ) + 
  geom_node_text(aes(x=x*1.15, y=y*1.15, label=names$complex), vjust=2, size=5, alpha=1) + 
  scale_edge_color_manual(values=c("darkblue", "gold", "#222222")) +  # "#2dae37"
  theme_void() + 
  theme(
    legend.position='none'
  )

dev.off()









