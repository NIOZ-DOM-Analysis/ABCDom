# First script to run
library(readr)
library(tidyverse)


#change depending on your working directory
wd.project <- "D:/Milou/Github/Projects/MGIA/2019/ABCDOM/R"


#define folder references
dirCyto <- "../Cytoscape"
dirDoc <- "../doc"
dirFigs <- "../figures"
dirOutput <- "../output"
dirR <- "../R"
dirRAW <- "../RAWdata"
dirWrite <- "../output/write_read"

setwd(wd.project)

#load metadata for all analysis

metadata <- read_csv(paste0(dirRAW, "/Metadata.csv"))

if( "Sample Name" %in% colnames(metadata) ){
  print("metadatasheet accepted")
}else{
  stop("No column with \"Sample Name\" as columnname. check metadata.csv ")}


# Scripts to run:
# Temperature
# Symbiont data
# FCM Microbial growth
# DOC data
# 16S_ASV_fasta_workup
# 16S_ASV_Analysis
# 16S_multivariate_Analysis
# 16S_multivariate_technical_replicate_Analysis ??
# Metabolomics_cleanup_analysis
# Procrustes

