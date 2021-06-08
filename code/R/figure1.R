# Figure 1 FNM
# SP (sebastian@pechmannlab.net)
# 2021

library(ggplot2)
library(cowplot)
library(reshape2)

setwd("~/FNM/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Schematic.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c_d50_k3 <- scan('data/motif/result_d50_k3_counts.txt')
c_d50_k4 <- scan('data/motif/result_d50_k4_counts.txt')
c_d50_k5 <- scan('data/motif/result_d50_k5_counts.txt')
c_d50_k6 <- scan('data/motif/result_d50_k6_counts.txt')

total <- rbind(c_d50_k3, c_d50_k4, c_d50_k5, c_d50_k6)
deg <- c(50, 50, 50, 50)
size <- c(3,4,5,6)

counts <- rbind(
  data.frame(counts=total[,1], deg=deg, size=size, class=rep("PPI", ncol(total)) ), 
  data.frame(counts=total[,2], deg=deg, size=size, class=rep("FNM", ncol(total)) ) )


svg(file = "figures/Figure1/B_counts.svg", height = 4, width = 3)

ggplot(counts, aes(x=factor(size), y=counts, fill=factor(class) ) ) + 
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("black", "#2988E2")) + 
  scale_y_log10() + 
  labs(y="Count", x="Size k") +
  theme_classic() + 
  theme(
    axis.text = element_text(size=16),
    axis.title = element_text(size=20), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(), 
    legend.text = element_text(size=14), 
    legend.position = c(0.256, 0.9)
  )

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


motif <- read.table("data/figures/figure1/motifDF2.txt", header=T, sep='\t')

m3 <- motif[,c(8, 1,5,2,3)]   # subselect columns for plotting
mm <- melt(m3)


svg(file = "figures/Figure1/C_motifs.svg", height = 8, width = 8)

ggplot(mm, aes(x=topo, y=value )) + 
  geom_col(aes(fill=variable), position="dodge2" ) +
  labs(x="", y="Count") + 
  scale_y_log10() +
  coord_polar() + 
  scale_fill_manual(values=c("black", "#777777", "#2988E2", "#8cd1ff")) + 
  theme_classic() +
  theme(
    text = element_text(size=16),
    axis.line.x = element_blank(), 
    legend.title = element_blank(), 
    legend.text = element_text(size=12)
  )
  
dev.off()









##################################################################################
# SUPPLEMENT
##################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL S1A ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


gi_df <- as.data.frame(read.table("data/figures/figure1/gihistogram.txt", header=T, sep='\t'))


svg(file = "figures/Figure1/supplement_A_pctgi.svg", height = 3, width = 4.5)

ggplot(gi_df, aes(x=bin*100, y=value ) ) + 
  geom_bar(aes( fill=factor(deg)), stat="identity", position="dodge") +
  geom_vline( aes(xintercept = 50), size=0.5, linetype="twodash" ) + 
  scale_y_continuous(breaks=c(0, 0.3, 0.6), labels=c(0, 0.3, 0.6)) + 
  theme_classic() + 
  facet_wrap(~k, ncol=1, nrow=4, strip.position="right") +
  scale_fill_manual(values=c("#333333", "#777777", "#BBBBBB"))+
  labs(y="Frequency", x="% Genetic interactions") +
  theme(
    axis.text = element_text(size=16),
    axis.title = element_text(size=20), 
    legend.title = element_blank()
  )

dev.off()


# "#00204DFF" "#00326FFF" "#31446BFF"




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL S1B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


rand_df <- as.data.frame(read.table("data/figures/figure1/rand.txt", header=T, sep='\t'))
rand_df$k <- as.factor(rand_df$k)


svg(file = "figures/Figure1/supplement_B_rand.svg", height = 4, width = 6)

ggplot(rand_df,  aes(x=N, y=rand_mean) ) + 
  geom_ribbon(aes(ymin=rand_mean-rand_std, ymax=rand_mean+rand_std), alpha=0.3) + 
  geom_line(size=1.2, color="black") + 
  facet_wrap(~k, ncol=2, nrow=2, strip.position="right", scales="free_y") +
  labs(y="Average count", x="Number of randomizations") + 
  theme_classic() + 
  theme(
    text = element_text(size=16)
  )

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL S1C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


c_d25_k3 <- scan('data/motif/result_d25_k3_counts.txt')
c_d25_k4 <- scan('data/motif/result_d25_k4_counts.txt')
c_d25_k5 <- scan('data/motif/result_d25_k5_counts.txt')
c_d25_k6 <- scan('data/motif/result_d25_k6_counts.txt')

c_d50_k3 <- scan('data/motif/result_d50_k3_counts.txt')
c_d50_k4 <- scan('data/motif/result_d50_k4_counts.txt')
c_d50_k5 <- scan('data/motif/result_d50_k5_counts.txt')
c_d50_k6 <- scan('data/motif/result_d50_k6_counts.txt')

c_d100_k3 <- scan('data/motif/result_d100_k3_counts.txt')
c_d100_k4 <- scan('data/motif/result_d100_k4_counts.txt')
c_d100_k5 <- scan('data/motif/result_d100_k5_counts.txt')
c_d100_k6 <- scan('data/motif/result_d100_k6_counts.txt')


total <- rbind(c_d25_k3, c_d50_k3, c_d100_k3, c_d25_k4,  c_d50_k4, c_d100_k4, c_d25_k5, c_d50_k5, c_d100_k5, c_d25_k6, c_d50_k6, c_d100_k6)
deg <- c(25, 50, 100, 25, 50, 100, 25, 50, 100, 25, 50, 100)
size <- c( rep("k=3", 3), rep("k=4", 3), rep("k=5", 3), rep("k=6", 3))

counts <- rbind(
  data.frame(counts=total[,1], deg=deg, size=size, class=rep("PPI", ncol(total)) ), 
  data.frame(counts=total[,2], deg=deg, size=size, class=rep("FNM", ncol(total)) ) )



svg(file = "figures/Figure1/supplement_C_counts.svg", height = 4, width = 6.5)

ggplot(counts, aes(x=factor(deg), y=counts, fill=factor(class) ) ) + 
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~size, ncol=4, nrow=1) +
  scale_fill_manual(values=c("black", "#2988E2")) + 
  scale_y_log10() + 
  labs(y="Count", x="Max degree") +
  theme_classic() + 
  theme(
    axis.text = element_text(size=16),
    axis.title = element_text(size=20), 
    axis.text.x = element_text(angle=90, hjust=1),
    legend.title = element_blank(), 
    legend.text = element_text(size=14)
  )

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PANEL S1D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


svg(file = "figures/Figure1/supplement_D_pct.svg", height = 4, width = 4)

counts.rel <- data.frame( perc=total[,2]/total[,1]*100, deg=deg, size=size)

ggplot(counts.rel, aes(x=factor(deg), y=perc ) ) + 
  geom_bar(stat="identity", position="dodge", fill="#2988E2") +
  facet_wrap(~size, ncol=4, nrow=1) +
  labs(y="% FNM", x="Max degree") +
  theme_classic() + 
  theme(
    axis.text = element_text(size=16),
    axis.title = element_text(size=20), 
    axis.text.x = element_text(angle=90, hjust=1),
    legend.title = element_blank()
  )

dev.off()





