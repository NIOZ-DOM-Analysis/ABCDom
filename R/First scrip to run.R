# First script to run
library(readr)
library(tidyverse)

getwd()
#change depending on your working directory
wd.project <- "D:/Milou/Github/ABCDom/R"


#define folder references
dirCyto <- "../Cytoscape"
dir_analogs_on <- file.path(dirCyto, "Analogs_on")
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


# Scripts to run in order:
# 1. Temperature
# 2. Metabolomics_cleanup_analysis   <- this will create the full_metadata object needed later
# 3. Symbiont data
# 4. FCM Microbial growth
# 5. DOC data
# 6. DO WE NEED THIS 16S_ASV_fasta_workup
# 6. 16S_multivariate_technical_replicate_Analysis
# 7. 16S_multivariate_Analysis
# 8. 16S_ASV_Analysis
# 9. 16S_multivariate_technical_replicate_Analysis ??
# 10. Procrustes


