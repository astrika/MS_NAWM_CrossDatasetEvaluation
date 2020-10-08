#Preparing RNA seq data
RNAseq <- read.table("C:/users/amanuel1/Dev/ZhaoLab/MS/GSE111972_norm_data.txt", header=T, as.is = T)
row.names(RNAseq) <- RNAseq[,1]
RNAseq <- RNAseq[,-c(1)]
RNAseq <- data.matrix(RNAseq)

#quantile normalization by Dave Tang, source: https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

#Log 2 transform, this is common for RNA seq analysis:
RNA_log2T <- (log2(RNAseq+1))

#performing quantile normalization, this is not common for RNA seq analysis, but may help for our analysis:
RNA_normal <- quantile_normalisation(RNA_log2T)

#plotting distributions to show the log2 transform, and normalized distributions:

par(mfrow=c(2,1))
boxplot(RNA_log2T, main="After log2 Tranformation")
boxplot(RNA_normal, main= "After log2 Transform & Quantile Normalization")


#separating cases vs control
GE.control <- RNA_normal[,c(seq(1,16))]
GE.case <- RNA_normal[,c(seq(17,31))]

#seperating sub-phenotypes (white matter vs grey matter samples)
GE.control.WM <- GE.control[,c(2,4,6,8,10,11,12,13,14,15,16)]
GE.control.GM <- GE.control[,-c(2,4,6,8,10,11,12,13,14,15,16)]

GE.case.WM <- GE.case[,c(2,4,6,8,10,11,12,13,14,15)]
GE.case.GM <- GE.case[,-c(2,4,6,8,10,11,12,13,14,15)]


