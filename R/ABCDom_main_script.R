'ABCDom main script'


#check if your wd is the R folder
getwd()
#change the wd.project to the current R folder/current folder
wd.project<-getwd()

# activate packages
source(paste0(dirR,'/packages.R'))

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
BiocManager::install("phyloseq")}
library(phyloseq)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

#first workup DOC data and symbiondiniacae
source(paste0(dirR, '/DOC_workup.R'))
source(paste0(dirR, '/Symbiodiniaceae_density_workup.R'))
source(paste0(dirR, '/FCM_workup.R'))


#then analyze data
source(paste0(dirR, '/DOC_analysis.R'))
source(paste0(dirR, '/Symbiodiniaceae_analysis.R'))
source(paste0(dirR, '/FCM_analysis.R'))

# 16S data

# prep the blank removed data.
file.edit(paste0(dirR, '/Cleaned_metabolite_Data_prep.R'))

# open the ordination of metabolite data here
file.edit(paste0(dirR, '/metabolite_ordination.R'))

file.edit(paste0(dirR, '/metabolite_diversity.R'))

#ABCDOM 2 analysis
# analyze the difference over time in the metabolite data
