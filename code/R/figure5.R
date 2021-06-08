# Figure 5 FNM
# SP (sebastian@pechmannlab.net)
# 2021

library(ggplot2)
library(reshape2)
library(scales)

setwd("~/FNM/")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 5A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cossimi <- as.data.frame(read.table("data/figures/figure5/cossimi.txt", header=T, sep='\t'))
cossimi <- cossimi[rowSums(cossimi>0.5)>0,]

cossimi.df <- data.frame(cbind(id=c(1:nrow(cossimi)), cossimi))
cm <- melt(cossimi.df, id='id')
cm$breaks <- cut(as.numeric(cm$value), c(-1, 0.55, 0.8, 1), right=F)

cm2 <- cm[cm$value > 0.5,] # reduce the file size for plotting


svg("figures/Figure5/A_heatmap.svg", height=5, width=4.5) #file size too big!

ggplot(cm2, aes(x=id, y=variable)) + 
  geom_tile(aes(fill=breaks)) + 
  coord_flip() + 
  labs(x="Motifs", y="Stress conditions") + 
  scale_fill_manual(values=c("white", "grey50", "red")) +
  theme_classic() +
  theme(
    text = element_text(size=20),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks=element_blank(), 
    legend.title=element_blank()
  )
  
dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 5B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cossimi2 <- as.numeric(as.vector(as.matrix(read.table("data/figures/figure5/cossimi.txt", header=T, sep='\t')) ))
cossimi2 <- data.frame( value=cossimi2)


svg("figures/Figure5/B_histogram.svg", height=2.5, width=4.5)

ggplot(cossimi2, aes(value)) + 
  geom_vline(xintercept = 0.55, size=0.5, linetype="twodash" ) +
  geom_vline(xintercept = 0.8, size=0.5, linetype="twodash" ) +
  geom_histogram(binwidth=0.025, fill="red", color="white") + 
  scale_x_continuous(limits=c(-0.5, 1)) + 
  labs(x="Co-regulation similarity", y="Count") +

  theme_classic() +
  theme(
    text = element_text(size=18)
  )

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 5C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cossimi3 <- as.matrix(read.table("data/figures/figure5/cossimi.txt", header=T, sep='\t'))
cossimiL <- scan("data/figures/figure5/motifL.txt")

cossimiDF <- rbind( data.frame(value=as.numeric(cossimi3[cossimiL==3,]), k=3),
                    data.frame(value=as.numeric(cossimi3[cossimiL==4,]), k=4),
                    data.frame(value=as.numeric(cossimi3[cossimiL==5,]), k=5),
                    data.frame(value=as.numeric(cossimi3[cossimiL==6,]), k=6) )
cossimiDF$k <- as.factor(cossimiDF$k)


svg("figures/Figure5/C_cdf.svg", height=2.5, width=4.5)

ggplot(cossimiDF, aes(x=value, colour=k)) + 
  stat_ecdf(size=1.25) + 
  labs(x="Co-regulation similarity", y="Fraction") + 
  scale_x_continuous(limits=c(-0.5, 1)) + 
  scale_color_manual(values=c("grey10", "grey30", "grey50", "grey70")) +
  theme_classic() + 
  theme(
    text = element_text(size=18),
    legend.position = c(0.9, 0.45)
  )

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 5D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


df_net <- as.data.frame(read.table("data/figures/figure5/motifs/Heat_shock_1_motif1_network.txt", header=T))
df_net$alpha <- ifelse(df_net$weight == 1, 1, 0.4)
df_net$width <- ifelse(df_net$weight == 1, 1, 0.9)
df_net$type <- as.factor(df_net$type)
vertices <- as.data.frame(read.table("data/figures/figure5/motifs/Heat_shock_1_motif1_vertices.txt", header=T))
vertices$angle <- ifelse(vertices$type =="other", 0, 90)
vertices$hjust <- ifelse(vertices$type =="other", 0, 1)
graph <- graph_from_data_frame(df_net, vertices=vertices$name)

edge_alpha <- rep(df_net$alpha, each=100)
edge_width <- rep(df_net$width, each=100)


svg("figures/Figure5/D_spliceosome.svg", height=3, width=4.5)

ggraph(graph, layout='stress') + 
  geom_edge_arc(aes(colour=df_net$type), width=edge_width, edge_alpha=edge_alpha, strength=df_net$strength) + 
  scale_edge_color_manual(values=c("gold", "skyblue", "#222222")) + 
  geom_node_point(aes(color=factor(vertices$type) ), size=vertices$size) +
  scale_color_manual(values=c("black")) +  
  geom_node_text(aes(x=x*1.15, y=y*1.15, label=vertices$gene), angle=vertices$angle, hjust=vertices$hjust, size=5, alpha=1) + 
  theme_void() + 
  theme(
    legend.position='none'
  )

dev.off()







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 5E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


df_net <- as.data.frame(read.table("data/figures/figure5/motifs/Heat_shock_1_motif11_network.txt", header=T))
df_net$alpha <- ifelse(df_net$weight == 1, 1, 0.6)
df_net$width <- ifelse(df_net$weight == 1, 1, 0.9)
vertices <- as.data.frame(read.table("data/figures/figure5/motifs/Heat_shock_1_motif11_vertices.txt", header=T))
vertices$angle <- ifelse(vertices$type =="other", 0, 90)
vertices$hjust <- ifelse(vertices$type =="other", 0, 1)
graph <- graph_from_data_frame(df_net, vertices=vertices$name)

edge_alpha <- rep(df_net$alpha, each=100)
edge_width <- rep(df_net$width, each=100)

svg("figures/Figure5/E_surveillance.svg", height=2, width=4)

ggraph(graph, layout='stress') + 
  geom_edge_arc(aes(colour=df_net$type), width=edge_width, edge_alpha=edge_alpha, strength=df_net$strength) + 
  scale_edge_color_manual(values=c("skyblue", "#222222")) + 
  geom_node_point(aes(color=factor(vertices$type) ), size=vertices$size) +
  scale_color_manual(values=c("black")) +  
  geom_node_text(aes(x=x*1.15, y=y*1.15, label=vertices$gene), angle=vertices$angle, hjust=vertices$hjust, size=5, alpha=1) + 
  theme_void() + 
  theme(
    legend.position='none'
  )

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ plot all stress coregulation motifs

files <- dir('data/figures/figure5/motifs/')
sel.net <- grep("network", files)
sel.node <- grep("vertices", files)

files.net <- files[sel.net]
files.node <- files[sel.node]

for (i in 1:length(files.net)){
  f_net <- files.net[i]
  f_node <- files.node[i]
  f_out <- paste("data/figures/figure5/motifs/", strsplit(f_net, "_network.txt"), ".svg", sep="")

  
  df_net <- as.data.frame(read.table(paste("data/figures/figure5/motifs/", f_net, sep=""), header=T))
  df_net$alpha <- ifelse(df_net$weight == 1, 1, 0.6)
  df_net$width <- ifelse(df_net$weight == 1, 1, 0.8)
  vertices <- as.data.frame(read.table(paste("data/figures/figure5/motifs/", f_node, sep=""), header=T))
  graph <- graph_from_data_frame(df_net, vertices=vertices$name)
  
  edge_alpha <- rep(df_net$alpha, each=100)
  edge_width <- rep(df_net$width, each=100)
  
  types_edges <- unique(df_net$type)
  if (length(grep("GIpos", types_edges))>0 && length(grep("GIneg", types_edges))>0 ){edge_col <- c("gold", "skyblue", "#222222") }
  if ( length(grep("GIpos", types_edges))>0 && length(grep("GIneg", types_edges))==0 ){edge_col <- c("gold", "#222222") }
  if ( length(grep("GIpos", types_edges))==0 && length(grep("GIneg", types_edges))>0 ){edge_col <- c("skyblue", "#222222") }
  
  
  g <- ggraph(graph, layout="stress") +  
    geom_edge_arc(aes(colour=df_net$type), width=edge_width, edge_alpha=edge_alpha, strength=df_net$strength) + 
    scale_edge_color_manual(values=edge_col) + 
    geom_node_point(aes(color=factor(vertices$type) ), size=vertices$size) +
    scale_color_manual(values=c("black")) +  
    geom_node_text(aes(x=x*1.15, y=y*1.15, label=vertices$gene), size=5, alpha=1) + 
    theme_void() + 
    theme(
      legend.position='none'
    )
  
  svg(file = f_out, height=5, width=5)
  plot(g)
  dev.off()
  
}








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 5F ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_net <- as.data.frame(read.table("data/figures/figure5/ubr_motif_network.txt", header=T))
df_net$alpha <- ifelse(df_net$weight == 1, 1, 0.6)
df_net$width <- ifelse(df_net$weight == 1, 1, 0.9)
vertices <- as.data.frame(read.table("data/figures/figure5/ubr_motif_vertices.txt", header=T))
vertices$angle <- ifelse(vertices$type =="other", 0, 90)
vertices$hjust <- ifelse(vertices$type =="other", 0, 1)
graph <- graph_from_data_frame(df_net, vertices=vertices$name)

edge_alpha <- rep(df_net$alpha, each=100)
edge_width <- rep(df_net$width, each=100)

lay <- create_layout(graph, layout='auto')
lay[match("YGR184C", lay$name),c(1,2)] <- c(3, 7)
lay[match("YGL058W", lay$name),c(1,2)] <- c(3, 5)
lay[match("YLR024C", lay$name),c(1,2)] <- c(3, 3)
lay[match("YBL041W", lay$name),c(1,2)] <- c(1, 2)
lay[match("YOR236W", lay$name),c(1,2)] <- c(3, 10)
lay[match("YDL020C", lay$name),c(1,2)] <- c(3, 0)

svg("figures/Figure5/F_ubr.svg", height=5, width=1.5)

ggraph(graph, layout=lay) + 
  geom_edge_arc(aes(colour=df_net$type), width=edge_width, edge_alpha=edge_alpha, strength=df_net$strength) + 
  scale_edge_color_manual(values=c("gold", "skyblue", "#222222")) + 
  geom_node_point(aes(color=factor(vertices$type) ), size=vertices$size) +
  scale_color_manual(values=c("black")) +  
  geom_node_text(aes(x=x*1.15, y=y*1.15, label=vertices$gene), angle=vertices$angle, hjust=vertices$hjust, size=5, alpha=1) + 
  theme_void() + 
  theme(
    legend.position='none'
  )

dev.off()



ubr <- as.data.frame(read.table("data/figures/figure5/ubr_df.txt", header=T))
ubr_profiles <- as.matrix( t(ubr[,1:8]) )
colnames(ubr_profiles) <- ubr$gene
ubr_profiles <- as.data.frame(ubr_profiles)
ubr_profiles$position <- 1:nrow(ubr_profiles)
ubr_profiles <- ubr_profiles[,c(5,1,2,3,4,6,7)]

profiles_df <- melt(ubr_profiles, id="position")
profiles_df$direction <- ifelse(profiles_df$value >= 0, "up", "down")


svg("figures/Figure5/F_profiles.svg", height =5, width = 3)

p0 <- ggplot(data=profiles_df) +
  geom_col(aes(x=position, y=value, fill=direction), alpha=0.5) +
  geom_hline(yintercept = 0, size=0.5, linetype="twodash" ) +
  geom_line(aes(x=position, y=value), size=1.2) +
  scale_fill_manual(values=c("green", "red")) +
  scale_x_continuous(breaks=c(1:8), labels=c("5min", "10min", "15min", "20min", "30min", "40min", "60min", "80min")) + 
  theme_classic() +
  labs(x = "", y = "Expression change (log2)") +
  theme(
    text = element_text(size=18), 
    axis.text.x = element_text(angle=90, vjust=1),
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )

p0 + facet_grid(rows = vars(variable))

dev.off()



