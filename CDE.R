# Astrid M Manuel
# Cross Dataset Evaluation of MS NAWM Gene Expression Profiles
# 05/28/2019
# Updated: 09/20/2020

#read in results of EW_dmGWAS permutation, for discovery set
#***********Please edit you working directory path
setwd("C:/Users/amanuel1/Dev/ZhaoLab//MS/GeneExpressionProfiles/EW/RNAseq/RNAseq_MS_NAWM_L1")

#***********Please edit your file destination of the EW_dmGWAS discovery set results
res.mat <- read.delim("RNAseq_MSNAWM_compiledPPI_L1.txt",as.is=T)
sapply(res.mat[,1], function(u)strsplit(u, split=" ")[[1]] ) -> module.genes

#Evaluation set edges
# ********Please edit your file destination of Evaluation Set Edge Weight and Node Weights
edge.weight <- read.delim("C:/Users/amanuel1/Dev/ZhaoLab/MS/GeneExpressionProfiles/EW/Microarray/Microarray_PLNAWM_L1/EdgeWeights_PLNAWMvsControl.txt",header=F,as.is=T)
node.weight <- read.delim("C:/Users/amanuel1/Dev/ZhaoLab/MS/GeneExpressionProfiles/EW/Microarray/Microarray_PLNAWM_L1/NodeWeights_PLNAWMvsControl.txt",header=F,as.is=T)

#recalculating module edge scores
eval_edgeSc <- c() 

for(i in seq(1:length(module.genes))){
  gene.list <- unlist(module.genes[i], use.names = F)
  scores <- edge.weight[which(edge.weight[,1] %in% gene.list & edge.weight[,2] %in% gene.list),3]
  eval_edgeSc <- c(eval_edgeSc, sum(scores))
}

eval_edgeSc_abs <- abs(eval_edgeSc)


#recalculating module scores
eval_moduleSC <- c()
dis_nodeSc <- res.mat$node_score

for(i in seq(1, length(dis_nodeSc))){
  eval_moduleSC <- c(eval_moduleSC, dis_nodeSc[i] + eval_edgeSc_abs[i])
}


dualEval <- res.mat[,c(1,7,4,5,3)]
dualEval$eval_edgeSc <- eval_edgeSc_abs
dualEval$eval_moduleSc <- eval_moduleSC

# based on module score
idx.edge.dis <- which(dualEval$edge_score > quantile(dualEval$edge_score, probs = 0.975))
idx.edge.eval <- which(dualEval$eval_edgeSc > quantile(dualEval$eval_edgeSc, probs = 0.975))
idx.edge.975percent <- idx.edge.dis[which(idx.edge.dis %in% idx.edge.eval)]

idx.dis <- which(dualEval$z_perm > quantile(dualEval$z_perm, probs = 0.975))
idx.eval <- which(dualEval$eval_moduleSc > quantile(dualEval$eval_moduleSc, probs = 0.975))
idx.975percent <- idx.dis[which(idx.dis %in% idx.eval)]

diff <- abs(dualEval$eval_edgeSc - dualEval$edge_score)
dualEval$edgeSc_diff <- diff
dualEval_sig <- dualEval[idx.975percent,]
write.table(dualEval[idx.975percent,], "DualEval_top97.5.txt", quote = F, row.names = F, sep = "\t")

par(mfrow=c(1,2))
plot(dualEval$z_perm, dualEval$eval_moduleSc, 
     col=ifelse(dualEval$z_perm > quantile(dualEval$z_perm, probs = 0.975) & dualEval$eval_moduleSc > quantile(dualEval$eval_moduleSc,probs = 0.975), "red", "black"), 
     pch=ifelse(dualEval$z_perm > quantile(dualEval$z_perm, probs = 0.975) & dualEval$eval_moduleSc > quantile(dualEval$eval_moduleSc, probs = 0.975), 19, 1), 
     xlab = "Discovery set z-score (RNAseq Permutation)", 
     ylab = "Evaluation Module Score")
plot(dualEval$edge_score, dualEval$eval_edgeSc,
     ylab = "Evaluation set edge score",
     xlab="Discovery set edge score")
points(dualEval$edge_score[idx.975percent], dualEval$eval_edgeSc[idx.975percent], col="red", pch=19)
points(dualEval$edge_score[c(1312)], dualEval$eval_edgeSc[c(1312)], col="red", pch=1, cex = 3.5)

#Graphing with Eval (Microarray) Edges
#Graphing most significant by evaluation score with Eval Edges
i <- 1312
  
  #genes of top module #i
  top.genes <- unlist(module.genes[i], use.names = F)
  #edges and edges weights for module i
  idx.ew <- which(edge.weight[,1] %in% top.genes & edge.weight[,2] %in% top.genes)
  ew <- edge.weight[idx.ew,]
  ew[,4] <- 2*pnorm(-abs(ew[,3]))
  ew[,5] <- -log10(ew[,4])
  ew.file <- ew[,c(1,2,5)]
  colnames(ew.file) <- c("source", "target", "edge weight")
  write.table(ew.file, paste("EW_DualEval_lowDiff_EvalEdge", i, ".txt", sep = ""), quote=F, row.names = F, sep = "\t")
  #nodes and node weights for module i
  idx.nw <- which(node.weight[,1] %in% top.genes)
  nw <- node.weight[idx.nw,]
  nw[,3] <- 2*pnorm(-abs(nw[,2]))
  nw[,4] <- -log10(nw[,3])
  nw.file <- nw[,c(1,4)]
  colnames(nw.file) <- c("key", "node weight")
  write.table(nw.file, paste("NW_DualEval_lowDiff_EvalEdge", i, ".txt", sep = ""), quote=F, row.names = F, sep = "\t")
  
  
  #Graphing with Discovery (RNAseq) edges
  #***********Please edit file directories for edge weight and node weight files of the discovery set
  edge.weight <- read.delim("C:/Users/amanuel1/Dev/ZhaoLab//MS/GeneExpressionProfiles/EW/RNAseq/RNAseq_MS_NAWM_L1/EdgeWeights_WM.txt",header=F,as.is=T)
  node.weight <- read.delim("C:/Users/amanuel1/Dev/ZhaoLab//MS/GeneExpressionProfiles/EW/RNAseq/RNAseq_MS_NAWM_L1/NodeWeights_WM.txt",header=F,as.is=T)
  #genes of top module #i
  top.genes <- unlist(module.genes[i], use.names = F)
  #edges and edges weights for module i
  idx.ew <- which(edge.weight[,1] %in% top.genes & edge.weight[,2] %in% top.genes)
  ew <- edge.weight[idx.ew,]
  ew[,4] <- 2*pnorm(-abs(ew[,3]))
  ew[,5] <- -log10(ew[,4])
  ew.file <- ew[,c(1,2,5)]
  colnames(ew.file) <- c("source", "target", "edge weight")
  write.table(ew.file, paste("EW_DualEval_lowDiff_DisEdge", i, ".txt", sep = ""), quote=F, row.names = F, sep = "\t")
  #nodes and node weights for module i
  idx.nw <- which(node.weight[,1] %in% top.genes)
  nw <- node.weight[idx.nw,]
  nw[,3] <- 2*pnorm(-abs(nw[,2]))
  nw[,4] <- -log10(nw[,3])
  nw.file <- nw[,c(1,4)]
  colnames(nw.file) <- c("key", "node weight")
  write.table(nw.file, paste("NW_DualEval_lowDiff_DisEdge", i, ".txt", sep = ""), quote=F, row.names = F, sep = "\t")

