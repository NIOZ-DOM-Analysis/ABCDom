"packages.R


Goal:  install and load all packages

Written by:
Milou Arts, NIOZ, NL, 2019

List of alterations:


"
# Specify your packages
package_list <- c(
  "ggplot2",
  "textclean",
  "rmarkdown",
  "knitr",
  "kableExtra",
  "tictoc",
  "expss",
  "vegan",
  "stringi",
  "psych",
  "nortest",
  "binom",
  "epitools",
  "car",
  "ape",
  "wesanderson",
  "RColorBrewer",
  "data.table",
  "DescTools",
  "broom",
  "readxl",
  "multcomp",
  "summarytools",
  "scales",
  "reshape",
  "reshape2",
  "cluster",
  "ggfortify",
  "rfPermute",
  "plyr",
  "tidyverse",
  "tibble",
  "dplyr",
  "svglite",
  "dunn.test",
  "UpSetR",
  "gridExtra",
  "grid",
  "ggpubr",
  "rstatix",
  "ComplexUpset",
  "cowplot",
  "scatterplot3d",
  "pdftools",
  "png",
  "magick",
  "devtools",
  "FSA",
  "pairwiseAdonis",
  "stats",
  "dendextend",
  "scales",
  "phyloseq",
  "DESeq2",
  "pheatmap",
  "openxlsx",
  "BiodiversityR")






not_installed <- package_list[!(package_list %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed)) install.packages(not_installed)                               # Install not installed packages
if(length(not_installed)>1){
  print(paste(length(not_installed), "packages had to be installed."))                    # Print information to console if packages had to be installed.
}

#unpack all packages to use
lapply(package_list, library, character.only = TRUE)

