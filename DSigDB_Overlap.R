# Astrid M Manuel
# February 27, 2020
# Overlapping 55 MS genes from prioritized gene network to Drug Signatures DataBase (DSigSB: http://dsigdb.tanlab.org/DSigDBv1.0/collection.html)

# Libraries used for reading gmt files:
install.packages("remotes")
remotes::install_github("mw201608/msigdb")
library(msigdb)

# Reading gmt files:

# FSA signature genes
FDA_DSigDB <- read.gmt("C:/Users/amanuel1/Dev/ZhaoLab/Drug_Repurposing/DSigDB_Tan/FDA.gmt.txt")
FDA <- unique(unlist(FDA_DSigDB$genesets))
# Kinase Inhibitors
KinaseI_DSigDB <- read.gmt("C:/Users/amanuel1/Dev/ZhaoLab/Drug_Repurposing/DSigDB_Tan/KinaseInhibitors.gmt.txt")
KinaseI <- unique(unlist(KinaseI_DSigDB$genesets))

CDE <- read.delim("C:/Users/amanuel1/Dev/ZhaoLab/MS/MSNAWM_DualEval/ResultsData/NW_TopModuleDualEval.txt", as.is = T)
CDE <- CDE$key

KinaseI_CDE <- genes_CDE[which(genes_CDE %in% KinaseI)]
FDA_CDE <- genes_CDE[which(genes_CDE %in% FDA)]


geneset_FDA <- FDA_DSigDB$genesets
FDA_MSassociatedPathways <- list()
geneset_list <- list()

i <- 1
for(i in seq(1:length(FDA_CDE))){
  FDA_gene <- FDA_CDE[i]
  j <- 1
  for(j in seq(1:length(geneset_FDA))){
    geneset <- unlist(geneset_FDA[j])
    if(FDA_gene %in% geneset){
      geneset_list[FDA_gene] <- list(geneset, geneset_list) 
    }
  }
  FDA_MSassociatedPathways <- list(geneset_list, FDA_MSassociatedPathways)
}

geneset_KI <- KinaseI_DSigDB$genesets

i <- 1
for(i in seq(1:length(KinaseI_CDE))){
  KI_gene <- KinaseI_CDE[i]
  j <- 1
  for(j in seq(1:length(geneset_KI))){
    geneset <- unlist(geneset_KI[j])
    if(KI_gene %in% geneset){
      print(geneset_KI[j])
    }
  }
}





#----------------------------------


drug_names <- c()

i <- 1
for(i in seq(1, length(geneset_FDA))){
  geneset <- unlist(geneset_FDA[i])
  j <- 1
  for(j in seq(1, length(genes_CDE))){
    if(genes_CDE[j] %in% geneset){
      drug <- names(geneset_FDA[i])
      drug_names <- c(drug_names, drug)
      j <- j + 1
    }
  }

 i <- i + 1
}

#------------------------------------

drug_names <- unique(drug_names)

M <- 55
N <- 13964
genes_per_drug <- list()
p <- c()

i <- 1
for(i in seq(1,length(drug_names))){
  drug <- drug_names[i]
  geneset <- unlist(geneset_FDA[drug], use.names = F)
  genes <- list(genes_CDE[which(genes_CDE %in% geneset)])
  genes_per_drug <- c(genes_per_drug, genes)
  k <- length(unlist(genes))
  s <- length(geneset)
  
  CT <- matrix(c(k,  s-k,M-k,N-(M-k)), nrow = 2, ncol = 2)
  FT <- fisher.test(CT)
  p <- c(p, FT$p.value)
}


bonferroni_p <- p.adjust(p, method = "bonferroni")
BH_p <- p.adjust(p, method = "BH")

FDA_enrichment <- cbind(drug_names, genes_per_drug, p, BH_p, bonferroni_p)


#-------------------------------------------


kinaseInhibitors <- c()
i <- 1
for(i in seq(1, length(geneset_KI))){
  geneset <- unlist(geneset_KI[i])
  j <- 1
  for(j in seq(1, length(genes_CDE))){
    if(genes_CDE[j] %in% geneset){
      kinase <- names(geneset_KI[i])
      kinaseInhibitors <- c(kinaseInhibitors, kinase)
      j <- j + 1
    }
  }
  
  i <- i + 1
}


kinaseInhibitors <- unique(kinaseInhibitors)

M <- 55
N <- 13964
genes_per_drug <- list()
p <- c()

i <- 1
for(i in seq(1,length(kinaseInhibitors))){
  kinase <- kinaseInhibitors[i]
  geneset <- unlist(geneset_KI[kinase], use.names = F)
  genes <- genes_CDE[which(genes_CDE %in% geneset)]
  genes_per_drug <- c(genes_per_drug, list(genes))
  k <- length(genes)
  s <- length(geneset)
  
  CT <- matrix(c(k,  s-k,M-k,N-(M-k)), nrow = 2, ncol = 2)
  FT <- fisher.test(CT)
  p <- c(p, FT$p.value)
}


KI_bonferroni_p <- p.adjust(p, method = "bonferroni")
KI_BH_p <- p.adjust(p, method = "BH")

KI_enrichment <- cbind(kinaseInhibitors, genes_per_drug, p, KI_BH_p, KI_bonferroni_p)

#---------------------------------------------------------------------------

CDE <- read.delim("C:/Users/amanuel1/Dev/ZhaoLab/MS/MSNAWM_DualEval/ResultsData/NW_TopModuleDualEval.txt", as.is = T)
CDE <- CDE$key

G_pathways <- read.delim("C:/Users/amanuel1/Dev/ZhaoLab/Drug_Repurposing/pathways.tsv", header =F, as.is = T)
drug_names_G <- unique(G_pathways$V1)

M <- 55
N <- 13964
genes_per_drug_G <- list()
genes_CDE_per_drug <- list()
p <- c()

i <- 1
for(i in seq(1,length(drug_names_G))){
 
  idx <- which(G_pathways$V1 == drug_names_G[i])
  geneset <- unique(c(G_pathways$V2[idx], G_pathways$V3[idx]))
  genes_per_drug_G <- c(genes_per_drug_G, list(geneset))
  
  CDE_geneset <- CDE[which(CDE %in% geneset)]
  genes_CDE_per_drug <- c(genes_CDE_per_drug, list(CDE_geneset))
  
  k <- length(CDE_geneset)
  s <- length(geneset)
  
  CT <- matrix(c(k,  s-k,M-k,N-(M-k)), nrow = 2, ncol = 2)
  FT <- fisher.test(CT)
  p <- c(p, FT$p.value)
  
}


G_bonferroni_p <- p.adjust(p, method = "bonferroni")
G_BH_p <- p.adjust(p, method = "BH")

G_enrichment <- cbind(drug_names_G, genes_CDE_per_drug, p, G_BH_p, G_bonferroni_p)


Drug_genesets <- cbind(drug_names_G, genes_per_drug_G)
write.table(Drug_genesets, "Genesets_per_drug.txt", quote = F, sep = "\t", row.names = F)