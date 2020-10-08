# Astrid M Manuel
# Calculating edge weights for EW_dmGWAS input, matching PPI edges to node weights.
# 05/10/2019

#gene-level p-values
node_weights <- read.table("C:/Users/amanuel1/Dev/ZhaoLab/MS/GeneExpressionProfiles/EW/MS_EUR_2018/MS_EUR_2018_pg.txt", as.is = T)

NW <-  read.table("C:/Users/amanuel1/Dev/ZhaoLab/MS/GeneExpressionProfiles/EW/MS_EUR_2018/NW_MS_EUR_2018.txt", as.is = T)
NW <- NW[which(NW[,2] < 1 & NW[,2] > 0),]

node_weights <- NW

#numeric matrix of gene expression values for control samples (gene IDs as row names & sample IDs as column names):
#GE.control <- read.table("PATH_TO_GENE_EXPRESSION_DATA_FOR_CONTROL", as.is = T)

#numeric matrix of gene expression values for control samples (gene IDs as row names & sample IDs as column names):
#GE.case <- read.table("PATH_TO_GENE_EXPRESSION_DATA_FOR_CASE", as.is = T)

GE.case <- PL_ChronicActive_GEmatrix
GE.control <- Control_GEmatrix

#Protein-prontein interaction network (2 columns for each interaction between 2 genes)
PPI <- read.table("C:/Users/amanuel1/Dev/ZhaoLab/MS/GeneExpressionProfiles/EW/EW_dmGWAS Documentation/EW_dmGWAS-Java/new.PPI_compiled.txt", as.is = T)


#Matching PPI connections to available gene expression values and node weights
genes.with.expression <- intersect(row.names(GE.case), row.names(GE.control))
genes.PPI <- unique(c(PPI[,1],PPI[,2]))
genes.with.NW <- NW[,1]
genes.with.expression <- intersect(row.names(GE.case), row.names(GE.control))

PPI <- PPI[which(PPI[,1] %in% genes.with.expression & PPI[,2] %in% genes.with.expression),]
PPI <- PPI[which(PPI[,1] %in% node_weights[,1] & PPI[,2] %in% node_weights[,1]),]


#Obtaining PCC
obtainPCC <- function(GE.matrix){
  PCCvalues <- cor(t(GE.matrix))
  return(PCCvalues)
}

#matching to PPI
matchPCCtoPPI <- function(PCCvalues, PPI){

  PPI_PCC <- c()
  
  for(i in seq(1:nrow(PPI))){
    PPI_PCC[i] <- PCCvalues[PPI[i, 1],PPI[i, 2]]
  }
  
  PPI_PCC_rmNA <- na.omit(PPI_PCC)
  
  return(PPI_PCC)
}


# Obtaining F(x): Fisher's transformation (Jia et al., 2014)
F.equation = function(x){
  (1/2) * log((1+x)/(1-x))
}



message("Obtaining CASE co-expression values...")
casePCC <- obtainPCC(GE.case)
message("Obtaining CONTROL co-expression values...")
controlPCC <- obtainPCC(GE.control)

message("Iterating through PPI network and matching co-expression values...")
control_PPI_PCC <- matchPCCtoPPI(controlPCC, PPI)
case_PPI_PCC <- matchPCCtoPPI(casePCC, PPI)  
case_PPI_PCC <- cbind(PPI, case_PPI_PCC)
control_PPI_PCC <- cbind(PPI, control_PPI_PCC)

holder_case <- case_PPI_PCC
holder_control <- control_PPI_PCC

# matching case vs control shared edges for edge-weight calculations

case_PPI_PCC <- case_PPI_PCC[-c(which(is.na(case_PPI_PCC[,3]))),]
control_PPI_PCC <- control_PPI_PCC[-c(which(is.na(control_PPI_PCC[,3]))),]

message("Matching case and control shared edges...")
apply(control_PPI_PCC, 1, function(u){
    paste( sort(c(u[1], u[2])), collapse="|")
  }) -> ctrl.str

  apply(case_PPI_PCC, 1, function(u){
    paste(sort(c(u[1], u[2])), collapse="|")
  }) -> case.str

  shared.edge.str = intersect(case.str, ctrl.str)
  
  match(shared.edge.str, ctrl.str) -> idx.1
  control_PPI_PCC.shared = control_PPI_PCC[idx.1, ] 
  
  match(shared.edge.str, case.str) -> idx.2
  case_PPI_PCC.shared = case_PPI_PCC[idx.2, ]
  
  
  
# obtaining edge weights
  
message("Calculating edge weights...")
  Ncase <- ncol(GE.case)
  Ncontrol <- ncol(GE.control)
  
  denominator = sqrt((1/(Ncase-3)) + (1/(Ncontrol-3))) 
  
  X <- c()
  
  for(k in 1:nrow(case_PPI_PCC.shared)){

    r.case <- case_PPI_PCC.shared[k,3]
    r.control <- control_PPI_PCC.shared[k,3]
    f.case = F.equation(r.case)
    f.control = F.equation(r.control)
    X[k] = (f.case - f.control)/denominator

  }

shared.edges <- case_PPI_PCC.shared[,c(1,2)]
EW_file <- cbind(shared.edges, X)
if(length(which(is.infinite(EW_file[,3]))) > 0){
EW_file <- EW_file[-c(which(is.infinite(EW_file[,3]))), ]
}
EW_file$X <- scale(EW_file$X)
EW_file$X <- abs(EW_file$X)
EW_file$X <- qnorm(1-(pnorm(EW_file$X,lower.tail=F)*2))



#getting the node weights and matching to ew
("Matching gene-level p-values...")
if(length(which(is.infinite(EW_file[,3]))) > 0){
  EW_file <- EW_file[-c(which(is.infinite(EW_file[,3]))), ]
}
if(length(which(is.na(EW_file[,3]))) > 0){
  EW_file <- EW_file[-c(which(is.na(EW_file[,3]))), ]
}

shared.genes <- intersect(c(EW_file[,1],EW_file[,2]),node_weights[,1])
shared.nodes.weight <- node_weights[match(shared.genes, node_weights[,1]), ]
qnorm(shared.nodes.weight[,2]/2, lower.tail = F) -> node.z
shared.nodes.weight$z <- node.z

NW_file <- shared.nodes.weight[,c(1,3)]
write.table(NW_file, "NodeWeights_PLNAWMvsControl.txt", row.names = F, col.names = F, quote = F, sep="\t")

shared.edge.weights <- EW_file[which(EW_file[,1] %in% shared.genes & EW_file[,2] %in% shared.genes), ]
EW_file <- shared.edge.weights
write.table(EW_file, "EdgeWeights_PLNAWMvsControl.txt", row.names = F, col.names = F, quote = F, sep="\t")


message("Success! Edge weight and node weight files have been created... You are ready to run EW_dmGWAS!")