# Figure 4 FNM
# SP (sebastian@pechmannlab.net)
# 2021

setwd("~/FNM/")

library(ggplot2)
library(reshape2)
library(MASS)
library(igraph)
library(ggraph)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 4A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ganno <- as.data.frame(read.table("data/figures/figure4/network_annotation_5.txt", header=T))
gmat <- read.table("data/figures/figure4/network_5.txt")
gmat <- as.matrix(gmat)
g <- graph_from_adjacency_matrix(gmat, 'undirected', diag=F)
#coords <- layout_with_graphopt(g, niter=500, charge=0.05, mass=20, spring.length=0, spring.constant=1)


ganno$color <- rep(NA, nrow(ganno))
ganno$color[ganno$type=="C"] <- "#222222"
ganno$color[ganno$type=="S"] <- "#AA00AA"

ganno$border <- rep(NA, nrow(ganno))   #rep("#222222", nrow(ganno))
#ganno$border[ganno$type=="C"] <- "#555555"

ganno$label <- ganno$name
ganno$label[which(ganno$type=="S")] <- NA

ganno$size <- ganno$betw * 150
ganno$size[ganno$type=="C"] <- 5

ganno$shape <- rep("circle", nrow(ganno))
ganno$shape[which(ganno$type=="C")] <- "square"


svg(file = "figures/Figure4/A_complex_network.svg", height=5, width=5)
plot(g, vertex.color=ganno$color, vertex.shape=ganno$shape, vertex.frame.color=ganno$border, vertex.label=NA, vertex.size=ganno$size, edge.width=0.3, edge.color="#222222") #, layout=coords)
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 4B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ganno_fnm_2 <- as.data.frame(read.table("data/figures/figure4/network_annotation_2.txt", header=T))
ganno_0_2 <- as.data.frame(read.table("data/figures/figure4/refnet_annotation_0_2.txt", header=T))
ganno_0_5 <- as.data.frame(read.table("data/figures/figure4/refnet_annotation_0_5.txt", header=T))
ganno_50_2 <- as.data.frame(read.table("data/figures/figure4/refnet_annotation_50_2.txt", header=T))
ganno_50_5 <- as.data.frame(read.table("data/figures/figure4/refnet_annotation_50_5.txt", header=T))

betw.fnm.2 <- ganno_fnm_2[ganno_fnm_2$type!='C',3]
betw.fnm <- ganno[ganno$type!='C',3]
betw.0.2 <- ganno_0_2[ganno_0_2$type!='C',3]
betw.0.5 <- ganno_0_5[ganno_0_5$type!='C',3]
betw.50.2 <- ganno_50_2[ganno_50_2$type!='C',3]
betw.50.5 <- ganno_50_5[ganno_50_5$type!='C',3]


betw.DF <- rbind(
    data.frame(betw=betw.fnm.2, type=rep("FNM", length(betw.fnm.2))),
    data.frame(betw=betw.50.2, type=rep("PPI", length(betw.50.2))),
    data.frame(betw=betw.0.2, type=rep("PPI_all", length(betw.0.2)))
)
betw.DF <- betw.DF[betw.DF$betw>0,]



svg(file = "figures/Figure4/B_betw.svg", height=3, width=5)

ggplot(betw.DF) + 
  geom_boxplot(aes(x=type, y=betw, fill=factor(type))) + 
  #geom_jitter(aes(x=type, y=betw)) +
  labs(x="", y="Betweenness") + 
  scale_fill_manual(values=c("#AA00AA", "#555555", "#BBBBBB")) + 
  scale_y_log10() +
  theme_classic() + 
  theme(
    text = element_text(size=20), 
    legend.title = element_blank(), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.text.x = element_blank()
  )

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 4C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data.betw <- data.frame(class=factor(c("FNM", "FNM", "PPI", "PPI", "PPI_all", "PPI_all")),
                        type=factor(c("all", "path", "all", "path", "all", "path")),
                        count=c(nrow(ganno_fnm_2), nrow(ganno_fnm_2[ganno_fnm_2$betw>0,]),
                                nrow(ganno_50_2[ganno_50_2$type!='C',]), nrow(ganno_50_2[ganno_50_2$betw>0,]),
                                nrow(ganno_0_2[ganno_0_2$type!='C',]), nrow(ganno_0_2[ganno_0_2$betw>0,]) )
)



svg(file = "figures/Figure4/C_counts.svg", height=3, width=5)

ggplot(data.betw) +
  geom_col(aes(x=class, y=count, fill=class), color="black", position="dodge2", size=0.8, linetype=c(1,5,1,5,1,5)) +
  scale_fill_manual(values=c("#AA00AA", "#555555", "#BBBBBB")) +
  labs(x="", y="Count") +
  theme_classic() +
  theme(
    text = element_text(size=20),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.title = element_blank()
  )

dev.off()








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 4D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


bn <- ganno[ganno$type=="S",]
idx.sorted <- sort(bn$betw,index.return=T, decreasing=F)
bn.top <- bn[idx.sorted$ix,]
bn.top$name <- factor(bn.top$name, levels=bn.top$name)

bn.25 <- bn.top[(nrow(bn.top)-24):nrow(bn.top),]



svg(file = "figures/Figure4/D_bottlencks.svg", height=5, width=3)

ggplot(bn.25) +
  geom_col(aes(x=name, y=betw)) + 
  labs(x="", y="Betweenness") + 
  coord_flip() + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(size=14, angle=90),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=20),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 4E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_net <- as.data.frame(read.table("data/figures/figure4/tfmetab_network.txt", header=T))
df_net$alpha <- ifelse(df_net$weight == 1, 1, 0.6)
df_net$width <- ifelse(df_net$weight == 1, 1, 0.8)
vertices <- as.data.frame(read.table("data/figures/figure4/tfmetab_vertices.txt", header=T))
vertices$angle <- ifelse(vertices$type =="other", 0, 90)
vertices$hjust <- ifelse(vertices$type =="other", 0, 1)
graph <- graph_from_data_frame(df_net, vertices=vertices$name)

edge_alpha <- rep(df_net$alpha, each=100)
edge_width <- rep(df_net$width, each=100)

lay <- create_layout(graph, layout='auto')  # create random layout and then change coords manually
###
lay[match("YCR073C", lay$name),c(1,2)] <- c(66, 6)
lay[match("YHL020C", lay$name),c(1,2)] <- c(66, 4)
lay[match("YDR007W", lay$name),c(1,2)] <- c(66, 10)
lay[match("YOL108C", lay$name),c(1,2)] <- c(66, 0)
#
lay[match("YNL325C", lay$name),c(1,2)] <- c(60, 7)
lay[match("YFR019W", lay$name),c(1,2)] <- c(60, 10)
lay[match("YOL001W", lay$name),c(1,2)] <- c(63, 6)
lay[match("YMR165C", lay$name),c(1,2)] <- c(63, 10)
lay[match("YNL027W", lay$name),c(1,2)] <- c(63, 0)
#
lay[match("YDR168W", lay$name),c(1,2)] <- c(60, 2)
lay[match("YDR463W", lay$name),c(1,2)] <- c(60, 0)
#
lay[match("YNL098C", lay$name),c(1,2)] <- c(57, 8)
lay[match("YOL081W", lay$name),c(1,2)] <- c(57, 6)
lay[match("YCR052W", lay$name),c(1,2)] <- c(54, 6)
lay[match("YBR035C", lay$name),c(1,2)] <- c(54, 10)
lay[match("YFL033C", lay$name),c(1,2)] <- c(57, 4)
lay[match("YJL005W", lay$name),c(1,2)] <- c(57, 10)
lay[match("YOR162C", lay$name),c(1,2)] <- c(57, 0)
#
lay[match("YNL157W", lay$name),c(1,2)] <- c(54, 2)
lay[match("YFR053C", lay$name),c(1,2)] <- c(51, 3)
lay[match("YOL097C", lay$name),c(1,2)] <- c(51, 10)
lay[match("YLR131C", lay$name),c(1,2)] <- c(51, 0)
lay[match("YGL253W", lay$name),c(1,2)] <- c(51, 5)
lay[match("YKL038W", lay$name),c(1,2)] <- c(48, 0)
#
lay[match("YKR064W", lay$name),c(1,2)] <- c(45, 7)
lay[match("YOR249C", lay$name),c(1,2)] <- c(45, 5)
lay[match("YDR110W", lay$name),c(1,2)] <- c(45, 3)
lay[match("YGR286C", lay$name),c(1,2)] <- c(45, 10)
lay[match("YKR055W", lay$name),c(1,2)] <- c(42, 7)
lay[match("YOL128C", lay$name),c(1,2)] <- c(42, 5)
lay[match("YER009W", lay$name),c(1,2)] <- c(42, 3)
lay[match("YKL216W", lay$name),c(1,2)] <- c(42, 10)
lay[match("YBR049C", lay$name),c(1,2)] <- c(42, 0)
#
lay[match("YGR184C", lay$name),c(1,2)] <- c(39, 7)
lay[match("YGL058W", lay$name),c(1,2)] <- c(39, 5)
lay[match("YLR024C", lay$name),c(1,2)] <- c(39, 3)
lay[match("YBL041W", lay$name),c(1,2)] <- c(36, 3)
lay[match("YOR236W", lay$name),c(1,2)] <- c(39, 10)
lay[match("YDL020C", lay$name),c(1,2)] <- c(39, 0)
#
lay[match("YKL048C", lay$name),c(1,2)] <- c(33, 7)
lay[match("YJR097W", lay$name),c(1,2)] <- c(33, 5)
lay[match("YHR205W", lay$name),c(1,2)] <- c(33, 3)
lay[match("YJR130C", lay$name),c(1,2)] <- c(33, 10)
lay[match("YCR065W", lay$name),c(1,2)] <- c(33, 0)
#
lay[match("YCR063W", lay$name),c(1,2)] <- c(30, 7)
lay[match("YDR482C", lay$name),c(1,2)] <- c(30, 5)
lay[match("YMR240C", lay$name),c(1,2)] <- c(30, 3)
lay[match("YIL162W", lay$name),c(1,2)] <- c(27, 10)
lay[match("YGL035C", lay$name),c(1,2)] <- c(27, 0)
#
lay[match("YOL005C", lay$name),c(1,2)] <- c(21, 7)
lay[match("YIL027C", lay$name),c(1,2)] <- c(24, 5)
lay[match("YCR040W", lay$name),c(1,2)] <- c(24, 3)
lay[match("YMR272C", lay$name),c(1,2)] <- c(24, 10)
lay[match("YMR043W", lay$name),c(1,2)] <- c(24, 0)
#
lay[match("YNL294C", lay$name),c(1,2)] <- c(18, 7)
lay[match("YML029W", lay$name),c(1,2)] <- c(18, 5)
lay[match("YGL045W", lay$name),c(1,2)] <- c(18, 3)
lay[match("YNL219C", lay$name),c(1,2)] <- c(18, 10)
lay[match("YBL021C", lay$name),c(1,2)] <- c(18, 0)
#
lay[match("YAR031W", lay$name),c(1,2)] <- c(12, 7)
lay[match("YCL025C", lay$name),c(1,2)] <- c(12, 5)
lay[match("YOR328W", lay$name),c(1,2)] <- c(15, 3)
lay[match("YLR100W", lay$name),c(1,2)] <- c(15, 10)
lay[match("YLR256W", lay$name),c(1,2)] <- c(15, 0)
#
lay[match("YLR305C", lay$name),c(1,2)] <- c(9, 6)
lay[match("YOR047C", lay$name),c(1,2)] <- c(9, 4)
lay[match("YHL032C", lay$name),c(1,2)] <- c(9, 10)
lay[match("YOR358W", lay$name),c(1,2)] <- c(9, 0)
#
lay[match("YGR063C", lay$name),c(1,2)] <- c(6, 6)
lay[match("YBL008W", lay$name),c(1,2)] <- c(6, 4)
lay[match("YBR166C", lay$name),c(1,2)] <- c(6, 10)
lay[match("YMR021C", lay$name),c(1,2)] <- c(6, 0)
#
lay[match("YDR245W", lay$name),c(1,2)] <- c(3, 7)
lay[match("YGL167C", lay$name),c(1,2)] <- c(0, 10)
lay[match("YER040W", lay$name),c(1,2)] <- c(0, 0)
###



svg("figures/Figure4/E_TFmetab.svg", height=4, width=16)

ggraph(graph, layout=lay) + 
  geom_edge_arc(aes(colour=df_net$type), width=edge_width, edge_alpha=edge_alpha, strength=df_net$strength) + 
  scale_edge_color_manual(values=c("gold", "skyblue", "#222222")) + 
  geom_node_point(aes(color=factor(vertices$type) ), size=vertices$size) +
  scale_color_manual(values=c("#760095", "black", "#4d74aa")) +  
  geom_node_text(aes(x=x*1.15, y=y*1.15, label=vertices$gene), angle=vertices$angle, hjust=vertices$hjust, size=5, alpha=1, repel=T) + 
  theme_void() + 
  theme(
    legend.position='none'
  )

dev.off()




