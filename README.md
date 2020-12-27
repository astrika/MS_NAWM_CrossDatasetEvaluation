# MS_NAWM_CrossDatasetEvaluation
Custom R scripts used to perform Cross Dataset Evaluation (CDE) are contained in this repository. The CDE has been applied to gene expression profiles from samples of multiple sclerosis (MS) normal appearing white matter (NAWM). 

CDE considers two gene expression datasets of a phenotype of interest, along with GWAS-based gene-level association scores for the phenotype of interest [1]. One gene expression dataset must be considered the discovery dataset, while the other is the evaluation dataset. In this repository, the phenotype of interest is MS. The discovery dataset includes RNA-seq expression profiles from MS NAWM microglia (and matched controls). The evaluation dataset considered is microarray profiles from perilesional MS NAWM (and matched controls).

# Necessary input for CDE: 
(1) Module list of discovery set
* The module list of discovery set is aquired from the output of the Edge-Weighted Dense Module Search of GWAS (EW_dmGWAS) tool [2]. Please see below for details about EW_dmGWAS.

(2) Edge scores of evaluation set
* Custom R scripts for calculation of edge scores of evaluation set, as well as edge scores for the MS NAWM evaluation set, are contained within the "EvaluationSet_Microarray" folder above.

# Running EW_dmGWAS
In order to run CDE, you must first aquire a module list of the discovery set. The module list of the discovery set is aquired from EW_dmGWAS output. EW_dmGWAS is a Java program, which is run through the command line. The input for EW_dmGWAS includes (1) GWAS-based node scores, (2) Expression-based edge weights. The "DiscoverySet_RNAseq" folder contains NodeWight_WM.txt file and EdgeWeight_WM.txt files utilized for running EW_dmGWAS in the context of MS NAWM.

You may visit https://github.com/fangfang0906/EW_dmGWAS-Java for instructions on how to run EW_dmGWAS on your command prompt.

For more information about EW_dmGWAS, you may visit https://bioinfo.uth.edu/dmGWAS/. If you use EW_dmGWAS, please cite [2].

# References
[1] Manuel AM, Dai Y, Freeman LA, Jia P, Zhao Z. A weighted network approach to combining genome-wide genetic variants and brain tissue expression identified intracellular viral pathways and potential drug targets in multiple sclerosis. Journal of Medical Genetics (Manuscript In Submission)

[2] Wang Q, Yu H, Zhao Z, et al. EW-dmGWAS: Edge-weighted dense module search for genome-wide association studies and gene expression profiles. Bioinformatics 2015;31:2591–4. doi:10.1093/bioinformatics/btv150
