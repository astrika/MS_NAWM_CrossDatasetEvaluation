#Run forMicroarray.R

#must install these packages:
install.packages("ggplot2")
library(ggplot2)
install.packages("gridExtra")
library(gridExtra)


#consider X case
case_GE <- PL_ChronicActive_GEmatrix
case_coexpression <- cor(t(case_GE))
case <- "Microarray_PLNAWM_"

# control gene expression and co-expression:
Control_GEmatrix <- Control_GEmatrix
Control_coexpression <- cor(t(Control_GEmatrix))

# consider top module of X module
XtopModule <- read.delim("C:/Users/amanuel1/Dev/ZhaoLab//MS/GeneExpressionProfiles/EW/Microarray/Microarray_PLNAWM_L1/EW_MS_NAWM_TopModule62.txt", as.is = T)

# where do you want images to save?
setwd("C:/Users/amanuel1/Dev/ZhaoLab//MS/GeneExpressionProfiles/EW/Microarray/Microarray_PLNAWM_L1")

# for each edge
for(i in seq(1:nrow(XtopModule))){
  gene1 <- XtopModule[i,1]
  gene2 <- XtopModule[i,2]
  
  gene1_control <- Control_GEmatrix[XtopModule[i,1], ]
  gene1_control <- unlist(c(gene1_control[,1:ncol(gene1_control)]), use.names = F)
  gene2_control <- Control_GEmatrix[gene2, ]
  gene2_control <- unlist(c(gene2_control[,1:ncol(gene2_control)]), use.names = F)
  
  gene1_case <- case_GE[gene1, ]
  gene1_case <- unlist(c(gene1_case[,1:ncol(gene1_case)]), use.names = F)
  gene2_case <- case_GE[gene2, ]
  gene2_case <- unlist(c(gene2_case[,1:ncol(gene2_case)]), use.names = F)
  
  dat_control <- data.frame(cbind(gene1_control, gene2_control))
  
  dat_case <- data.frame(cbind(gene1_case, gene2_case))
  
  
  
  p1 <- ggplot(dat_case, aes(x=gene1_case, gene2_case)) + geom_point(shape=19, col="red") + geom_smooth(method=lm, col="red") + labs(title=paste("PCC =", substr(as.character(case_coexpression[gene1,gene2]), 1, 5),"\n"), x=gene1, y=gene2)

  p2 <- ggplot(dat_control, aes(x=gene1_control, y=gene2_control)) + geom_point(shape=19, col="blue") + geom_smooth(method=lm, col="blue") + labs(title=paste("PCC =", substr(as.character(Control_coexpression[gene1,gene2]), 1, 5),"\n"), x=gene1, y=gene2)
  
  full <- grid.arrange(p1,p2, ncol=1)
  
  ggsave(filename = paste(case, "M62_Edge", i, ".jpg"), plot=full, width= 2, height = 8, unit="in", dpi=150)

}
