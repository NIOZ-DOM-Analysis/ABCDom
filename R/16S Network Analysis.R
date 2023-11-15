#Have to restart R before running this script due to package incompatabilities with phyloseq and SpiecEasi.

#load libraries
library(devtools)
library(SpiecEasi)
library(igraph)

#unload phyloseq
#detach("package:phyloseq", unload=TRUE)
remove.packages("phyloseq)")

#-----------------------------------------------------------------------------------------------------------------------------
#load custom functions
source("make.relabund.R")

#-----------------------------------------------------------------------------------------------------------------------------
#Work up data and metadata
abund.nosub.tend <- abund.nosub2[,colnames(abund.nosub2) %in% metadata.tend$Sample_Name_Unique] #subset for just tend samples, including water samples

abund.nosub.tend.cull <- abund.nosub.tend[rownames(abund.nosub.tend) %in% taxonomy.cull$OTUNumber,] #subset for just 159 ASVs run thru DESEQ2

#-----------------------------------------------------------------------------------------------------------------------------
#Construct a network in spieceasi and igraph
net.c1 <- spiec.easi(as.matrix(t(abund.nosub.tend.cull)), method='mb')

net.c1.ig <- adj2igraph(as.matrix(getRefit(net.c1)), vertex.attr=list(name=taxonomy.cull$OTUNumber))

#export for use in supplemental figure 8 in cytoscape
#write.graph(net.c1.ig, file="FigS8_network.txt",format="ncol")
