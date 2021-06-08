# BiRewire

#setwd("FNM/code/")
library(BiRewire)


data <- as.matrix(read.table("tmp.inp"))
diag(data) <- 0

data.rewired <- birewire.rewire.undirected(data, max.iter="n",accuracy=0.00005, verbose=TRUE,MAXITER_MUL=10,exact=FALSE)

write.table(data.rewired, "tmp.rewired", col.names=F, row.names=F)

q(save='no')
