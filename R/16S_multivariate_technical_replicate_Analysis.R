#load libraries
library(vegan)
library(ggplot2)

#-----------------------------------------------------------------------------------------------------------------------------
#load custom functions
source(file="generate.square.dist.R")

#-----------------------------------------------------------------------------------------------------------------------------
#load metadata and data

unifrac <- read.csv(file.path(dirRAW, "16S/ASRAAMP_MCR2019_16S_Mar_2021_2k_subsample_04212021", "otu_repr_100.tre1.wsummary.csv"))

metadata_tech_replicates <- read.csv(file.path(dirRAW, "16S/ASRAAMP_MCR2019_16S_Mar_2021_2k_subsample_04212021", "metadata_tech_replicates_16S.csv"))

#-----------------------------------------------------------------------------------------------------------------------------
#Work up the abcDOM tech replicate metadata.
metadata_tech_replicates$Experiment <- as.factor(metadata_tech_replicates$Experiment) #convert to factor

metadata_abc_replicates <- subset(metadata_tech_replicates, Experiment=="ABC") #subset for just abcDOM samples

metadata_16S <- subset(metadata, Sample_Type=="Sterivex") #subset metadata for just sterivex samples.

metadata_abc_replicates1 <- merge(metadata_abc_replicates, metadata_16S, by.x="Sample_Name", by.y="Sample Name", all.x=T, all.y=F) #merge metadata_abc_replicate1 with metadata.

metadata_abc_replicates2 <- metadata_abc_replicates1[-69:-70,] #remove samples that weren't collected.

#export metadata_abc_replicate2
#CONSIDER REMOVING write.csv(metadata_abc_replicates2, "metadata_abc_replicates2.csv") #CHANGE DIRECROREY FOR FINAL ONE

metadata_abc_replicate_tend <- subset(metadata_abc_replicates2, Timepoint_char=="Tend") #subset for just tend.

#-----------------------------------------------------------------------------------------------------------------------------
#Work up unifrac data.
generate.square.dist.matrix(unifrac, metadata_tech_replicates, 1) #work up tech_replicates data

replicate_unifrac.dist.matrix <- dist.matrix1 #store output as new data frame

replicate_unifrac.dist <- square.dist.matrix #store output as new data frame

generate.square.dist.matrix(unifrac, metadata_abc_replicates2, 2) #work up abc_replicate2 data

abc_replicate2_unifrac.dist.matrix <- dist.matrix1 #store output as new data frame

abc_replicate2_unifrac.dist <- square.dist.matrix #store output as new data frame

generate.square.dist.matrix(unifrac, metadata_abc_replicate_tend, 2) #work up abc_replicate_tend data

abc_replicate_tend_unifrac.dist.matrix <- dist.matrix1 #store output as new data frame

abc_replicate_tend_unifrac.dist <- square.dist.matrix #store output as new data frame

#-----------------------------------------------------------------------------------------------------------------------------
#Generate up NMDS objects from the unifrac square dist matrices.
nmds.replicate <- metaMDS(replicate_unifrac.dist, k=2, trymax=100)

nmds.abc.replicate2 <- metaMDS(abc_replicate2_unifrac.dist, k=2, trymax=100)

nmds.abc.tend <- metaMDS(abc_replicate_tend_unifrac.dist, k=2, trymax=100)

#-----------------------------------------------------------------------------------------------------------------------------
#Visualize NMDSes, coloring ordinations by sample type.
colvec <- c("blue","green","red","yellow","purple") #set colors for sample types

#plot
plot(nmds.replicate, type="n")
points(nmds.replicate, display="sites",
       col=colvec[as.integer(metadata_tech_replicates$Experiment)])
legend("topleft", legend=c(levels(metadata_tech_replicates$Experiment)), pch=c(16,16,16,16), col=colvec, bty="n", cex=.5, y.intersp = .5)
#PCR positive and negative controls cluster far away from abcDOM samples (green).

#Now, visualize just the abcDOM samples, colored by replicate with replicates connected by an arrow.
colvec2 <- c("blue","red") #set colors for replicates

#plot
metadata_abc_replicates2$Replicate <- as.factor(metadata_abc_replicates2$Replicate) #convert to factor

plot(nmds.abc.replicate2, type="n")
points(nmds.abc.replicate2, display="sites",
       col=colvec2[metadata_abc_replicates2$Replicate],
       pch=16)
ordiarrows(nmds.abc.replicate2,groups=metadata_abc_replicates2$Sample_Name)
#The left hand cluster corresponds to the T0 samples, which seem to group nicely. The right hand cluster corresponds to the Tfinal samples, for which some technical replicates do not group together.

#-----------------------------------------------------------------------------------------------------------------------------
#Test the effect of replicate on community structure.
permanova.replicate <- adonis2(abc_replicate2_unifrac.dist~Replicate, data=metadata_abc_replicates2, by="margin", permutations=999) #generate model

permanova.replicate #No significant effect and very low R2 for technical replicates.

#-----------------------------------------------------------------------------------------------------------------------------
#Visualize NMDS for just the tend samples.
cost.col.fill<-c("dodgerblue1","firebrick1", "white", "white", "grey70", "grey70") #set fill vec

cost.col.line<-c("dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3") #set outline vec

pointvec <- c(21,24) #set point shape vector

#metadata_abc_replicate_tend$Treatment <- as.factor(metadata_abc_replicate_tend$Treatment) #CONSIDER REMOVING
metadata_abc_replicate_tend$Treatment <- factor(metadata_abc_replicate_tend$Treatment, levels = c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", 'Bleached + Heated', "Ambient Water Control", "Heated Water Control")) #convert Treatment to factor

#plot
plot(nmds.abc.tend,type="n")
points(nmds.abc.tend,display="sites",
       col=cost.col.line[metadata_abc_replicate_tend$Treatment],
       bg=cost.col.fill[metadata_abc_replicate_tend$Treatment],
       cex=2,
       pch=21)
ordiarrows(nmds.abc.tend,groups=metadata_abc_replicate_tend$Sample_Name)
ordilabel(nmds.abc.tend,labels=metadata_abc_replicate_tend$Sample_Name_Unique,border=NA,fill=NA,cex=.5)
#Technical replicates for ABC_055 and ABC_059 did not cluster well together.

#-----------------------------------------------------------------------------------------------------------------------------
#Identify statistical outliers n tend samples.
lower_tri_abc_replicate_tend_unifrac.dist.matrix <- abc_replicate_tend_unifrac.dist.matrix[lower.tri(abc_replicate_tend_unifrac.dist.matrix, diag=T)] #extract lower triangle from unifrac data

hist(lower_tri_abc_replicate_tend_unifrac.dist.matrix) #visualize histogram. Two outliers identified on far right.

sd_lower_tri <- sd(lower_tri_abc_replicate_tend_unifrac.dist.matrix) #calculate SD of pairwise unifrac distances

mean_lower_tri <- mean(lower_tri_abc_replicate_tend_unifrac.dist.matrix) #calculate mean of pairwise unifrac distances

mean_lower_tri + 3*sd_lower_tri #cutoff will be a maximum unifrac dist between tech replicates of mean(pairwise unifrac distances) + 3*SD(pairwise unifrac distances). This cutoff is .7227101. ABC_055 and ABC_059 fall above this cutoff.
#rerun 16S pipeline merging all tech replicates except ABC_055 and ABC_059

#-----------------------------------------------------------------------------------------------------------------------------
#clean up
#CONSIDER REMOVINGrm(list=ls())


