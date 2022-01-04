# Figure 2 FNM
# SP (sebastian@pechmannlab.net)
# 2021


library(ggplot2)
library(reshape2)
library(cowplot)
library(igraph)
library(ggraph)

setwd("~/FNM/")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 2A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

supp <- read.table("data/figures/figure2/yeast_suppressors.txt", header=T)
supp <- supp[,c(1,2)]
supp <- data.frame(cat=factor(c("FNM", "PPI", "CTRL"), levels=c("FNM", "PPI", "CTRL")), pct=(supp[,2]/supp[,1]) )


svg(file = "figures/Figure2/A_suppressor.svg", height=4, width=2.5)

ggplot(supp) + 
  geom_col(aes(x=cat, y=pct*100, fill=cat), position="dodge2") + 
  theme_classic() +
  scale_fill_manual(values=c("#2988E2", "#555555", "#222222")) + 
  labs(x="", y="% with suppressor interactions") + 
  theme(
    text = element_text(size=18),
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )

dev.off()






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 2B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gencat <- read.table("data/figures/figure2/gencatDF.txt", header=T)

scores.ess.fnm <- scan("data/figures/figure2/FNM.scores.ess.txt")
scores.ess.all <- scan("data/figures/figure2/ALL.scores.ess.txt")
scores.sub.fnm <- scan("data/figures/figure2/FNM.scores.sub.txt")
scores.sub.all <- scan("data/figures/figure2/ALL.scores.sub.txt")

# only take random 10% of the data for plotting to reduce the plot file size! distributions are sufficiently identical
scores.ess.fnm2 <- scores.ess.fnm[ sample(c(1:length(scores.ess.fnm)), round(0.1*length(scores.ess.fnm) ,0) ) ]
scores.ess.all2 <- scores.ess.all[ sample(c(1:length(scores.ess.all)), round(0.1*length(scores.ess.all) ,0) ) ]
scores.sub.fnm2 <- scores.sub.fnm[ sample(c(1:length(scores.sub.fnm)), round(0.1*length(scores.sub.fnm) ,0) ) ]
scores.sub.all2 <- scores.sub.all[ sample(c(1:length(scores.sub.all)), round(0.1*length(scores.sub.all) ,0) ) ]


gencat.ess <- rbind(
      data.frame(value=scores.ess.fnm2, cat=rep("FNM", length(scores.ess.fnm2)) ),
      data.frame(value=scores.ess.all2, cat=rep("PPI", length(scores.ess.all2)) )
  )
  
gencat.sub <- rbind(
  data.frame(value=scores.sub.fnm2, cat=rep("FNM", length(scores.sub.fnm2)) ),
  data.frame(value=scores.sub.all2, cat=rep("PPI", length(scores.sub.all2)) )
)


plot.ess <- ggplot(gencat.ess) + 
  geom_boxplot(aes(x=cat, y=value*100, fill=cat)) + 
  labs(x="", y="% essential") +
  scale_fill_manual(values=c("#2988E2", "#555555")) + 
  theme_classic() +
  theme(
    text = element_text(size=20),
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )

plot.sub <- ggplot(gencat.sub) + 
  geom_boxplot(aes(x=cat, y=value*100, fill=cat)) + 
  labs(x="", y="% complex subunit") +
  scale_fill_manual(values=c("#2988E2", "#555555")) + 
  theme_classic() +
  theme(
    text = element_text(size=20),
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )

svg(file = "figures/Figure2/B_gencat.svg", height=4, width=4)

plot_grid(plot.ess, plot.sub, labels ="", ncol = 2, align = 'h')

dev.off()






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 2C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

summary <- read.table("data/figures/figure2/motif_summary.txt", header=T, sep='\t')
summary$cat <- factor(summary$cat, levels=summary$cat)


plot.orf <- ggplot(summary) + 
  geom_col(aes(x=cat, y=orf, fill=cat)) + 
  labs(x="", y="Number of ORFs") +
  scale_fill_manual(values=c("grey50", "red")) + 
  theme_classic() +
  theme(
    text = element_text(size=20),
    axis.text.x = element_text(angle=45, size=20, vjust=0.5, hjust=0.75),
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )

plot.motif <- ggplot(summary) + 
  geom_col(aes(x=cat, y=motif, fill=cat)) + 
  labs(x="", y="Number of motifs") + 
  scale_fill_manual(values=c("grey50", "red")) + 
  theme_classic() +
  theme(
    text = element_text(size=20),
    axis.text.x = element_text(angle=45, size=20, vjust=0.5, hjust=0.75),
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )



svg(file = "figures/Figure2/C_summary.svg", height=4, width=3.8)

plot_grid(plot.motif, plot.orf, labels ="", ncol = 2, align = 'h')

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL 2D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


degree <- as.data.frame(read.table("data/figures/figure2/degree.txt", header=T, sep='\t'))
orfs_noclust <- read.csv("data/figures/figure2/list_orfs_noclust.txt", header=T)

sel_noclust <- rep(NA, nrow(orfs_noclust))
for (i in 1:nrow(orfs_noclust)){ sel_noclust[i] <- which(degree$ORF == as.character(orfs_noclust$ORF[i]) ) }

clustdf <- rbind(
  data.frame(class=rep("clust", nrow(degree)-length(sel_noclust)), N=degree$N[-sel_noclust], deg=degree$degG[-sel_noclust]), 
  data.frame(class=rep("noclust", length(sel_noclust)), N=degree$N[sel_noclust], deg=degree$degG[sel_noclust])
  )
           
    

plot.N <- ggplot(clustdf, aes(y=N)) + 
  geom_boxplot(aes(fill=class)) + 
  labs(x="", y="Count")+
  scale_fill_manual(values=c("grey50", "red") ) + 
  scale_y_log10() +
  theme_classic() + 
  theme(
    text = element_text(size=18), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(), 
    legend.title = element_blank(),
    legend.position = 'bottom'
  )

plot.deg <- ggplot(clustdf, aes(y=deg)) + 
  geom_boxplot(aes(fill=class)) +
  labs(x="", y="Degree") +
  scale_fill_manual(values=c("grey50", "red") ) + 
  scale_y_continuous(breaks=c(0, 5, 10, 20, 30, 40, 50), labels=c(0, 10000, 10, 20, 30, 40, 50)) +  #weird hack to get spacing, delete the 10000 in illustrator/inkscape
  theme_classic() + 
  theme(
    text = element_text(size=18), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.title = element_blank(), 
    legend.position = 'none'
  )



svg(file = "figures/Figure2/D_clust.svg", height=4, width=4)

plot_grid(plot.N, plot.deg, labels ="", ncol = 2, align = 'h')

dev.off()




#wilcox.test(clustdf$deg[clustdf$class=="clust"], clustdf$deg[clustdf$class=="noclust"])$p.value
#[1] 2.315483e-30
#wilcox.test(clustdf$N[clustdf$class=="clust"], clustdf$N[clustdf$class=="noclust"])$p.value
#[1] 5.300762e-134






##################################################################################
# SUPPLEMENT
##################################################################################


library(igraph)

for (i in 1:12){

  f_mat <- paste("data/figures/figure2/clusters/cluster.", i, ".txt", sep="")
  f_names <- paste("data/figures/figure2/clusters/cluster.", i, "_orfs.txt", sep="")
  f_out <- paste("figures/Figure2/clusters/cluster_", i, ".svg", sep="")

  gmat <- as.matrix(read.table(f_mat))
  names <- read.table(f_names, header=T)
  colnames(gmat) <- names$ORF
  rownames(gmat) <- names$ORF
  g <- graph_from_adjacency_matrix(gmat, 'undirected', diag=F, add.rownames="code")

  svg(file = f_out, height=5, width=5)
  plot(g, vertex.color="black", vertex.size=10, edge.width=3, edge.color="black", vertex.label.dist=1.8, vertex.label.color="red", vertex.label.cex=1.5)
  dev.off()

}




