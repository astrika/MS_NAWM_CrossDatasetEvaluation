# Astrid M Manuel 1/16/19
# Handling duplicate probes from gene expression matrix

# Loading the normalized gene expression values with matched gene symbols
ExprsNM <- read.table("c:\\users\\amanuel1\\Dev\\ZhaoLab\\MS\\GeneExpressionProfiles\\GSE108000_GeneExprsionMatrix_Normalized_Matched.txt")

# determining duplicates
# help from: https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r
geneProbeCounts <- ExprsNM[duplicated(ExprsNM$tT...7.) | duplicated(ExprsNM$tT...7., fromLast = TRUE),]
numdups <- table(geneProbeCounts$tT...7.)
actualdups <- numdups[numdups > 0]
dupgenes <- rownames(actualdups)

# removing duplicate gene IDs (keeping a the first gene record for each gene that is duplicated, later to be replaced by aquired average)
removeduplicates <- ExprsNM[!duplicated(ExprsNM$tT...7.),]
nodupGeneRecords <- removeduplicates[,2:41]
row.names(nodupGeneRecords) <- removeduplicates$tT...7.

# obtaining average values per duplicate genes
for(gene in dupgenes){
  
  #creating matrix of duplicate records for gene
  duplicateRecords <- ExprsNM[ExprsNM$tT...7. == gene, ]
  
  # removing non-numeric value (gene IDs in 1st column)
  duplicateRecords <- duplicateRecords[,2:41]
  
  # average for duplicate records:
  geneProbeAverage <- colMeans(duplicateRecords)
  
  # replacing the previously held gene record
  nodupGeneRecords[gene,] <- geneProbeAverage
  
}

# write.table(nodupGeneRecords, "GSE108000_NoDupGeneRecords.txt")


Control_GEmatrix <- nodupGeneRecords[,1:10]
RIM_ChronicActive_GEmatrix <- nodupGeneRecords[,11:17]
PL_ChronicActive_GEmatrix <- nodupGeneRecords[,18:24]
RIM_Inactive_GEmatrix <- nodupGeneRecords[,25:32]
PL_Inactive_GEmatrix <- nodupGeneRecords[,33:40]