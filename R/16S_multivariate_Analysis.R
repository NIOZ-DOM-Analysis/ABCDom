#load libraries
library(ggplot2)
library(vegan)
library(pairwiseAdonis)
library(stats)
library(dendextend)
library(scales)
library(FSA)
library(graphics)
library(grDevices)

#source custom functions
source(file="generate.square.dist.script.07.27.2020.R")

#read in metadata and data
#Read in unifrac.
unifrac <- read.csv(file.path(dirRAW, "16S", "otu_repr_100.tre1.wsummary.csv"))

#read in ASV table
abund <- read.csv(file.path(dirRAW, "16S", "abundance_table_100.database.csv"))

#Read in metadata. Unnecessary since metadata is included in global environment
#metadata <- read.csv(file="metadata1.csv")
#Read in the 16S metadata
metadata_16S <- read.csv(file.path(dirOutput, "metadata_16S_abcDOM.csv"))

#Work up metadata and unifrac data.
#subset metadata for only samples that made it through the 16S pipeline
samples_16s <- colnames(abund)[colnames(abund) %in% metadata_16S$Sample_Name_Unique] #create a vector of samples from metadata_16S that made it through the 16S pipeline.
metadata_16S_1 <- metadata_16S[metadata_16S$Sample_Name_Unique %in% samples_16s,] #subset metadata 16S for these samples.

#Remove outliers. Outlier tech replicates defined in 16S_multivariate_technical_replicate_Analysis.R. Include technical replicates that are closer to the centroid of their treatment level.
metadata.outliers <- metadata_16S_1[metadata_16S_1$Sample_Name_Unique!="ABC_055" & metadata_16S_1$Sample_Name_Unique!="ABC_059" & metadata_16S_1$Sample_Name_Unique!="ABC_MOC_B",] #Remove samples ABC_055 and ABC_059 because these tech reps were weird, include the v1 versions. Remove ABC_MOC_B beause I dont know what it is.

#generate a new factor that combines Treatment1 and Timepoint.
metadata.outliers$Treatment_Timepoint <- paste(metadata.outliers$Treatment, metadata.outliers$Timepoint, sep="_")

#Now define and remove potential biological replicate outliers.
#generate unifrac dist matrix from metadata.outliers.
generate.square.dist.matrix(unifrac,metadata.outliers,3)
unifrac.dist.matrix.outliers=dist.matrix1
unifrac.dist.outliers=as.dist(unifrac.dist.matrix.outliers)

#use betadisper to define within group variability.
dispersion <- betadisper(unifrac.dist.outliers, group=metadata.outliers$Treatment_Timepoint, type="centroid")
boxplot(dispersion)
dist_centroid <- as.data.frame(dispersion$distances) #store distance to centroid values in df
dist_centroid$sample <- rownames(dist_centroid) #add a sample column
dist_centroid <- merge(dist_centroid, metadata.outliers, by.x="sample", by.y="Sample_Name_Unique") #merge in metadata

#Assess distribution of distance from centroid values.
hist(log10(dist_centroid$`dispersion$distances`)) #needs to be log10 transformed to be normal

dist_centroid$log10_dist <- log10(dist_centroid$`dispersion$distances`) #add a log10 transformed dis centroid variable to df.

mean(log10(dist_centroid$`dispersion$distances`))+1.5*sd(log10(dist_centroid$`dispersion$distances`)) #Outliers defined as log10(dist) value greater than 1.5*sd(log10(dist))+mean(log10(dist)) = -.627
#This defines ABC072 and 064 as outliers.

#remove outliers
metadata1 <- metadata.outliers[c(-26,-34),]

#update treatment factor levels
metadata1$Treatment <- factor(metadata1$Treatment, levels=c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", "Bleached + Heated", "Ambient Water Control", "Heated Water Control"))

#subset for various treatments/timepoints.
metadata.t0=metadata1[metadata1$Timepoint_char=="T0",] #subset for only t0 samples
metadata.tend=metadata1[metadata1$Timepoint_char=="Tend",] #subset for only the tend samples
metadata.coral.tend=metadata.tend[metadata.tend$Origin_PlanC!="control",] #subset for only tend coral samples.

#Work up unifrac.
#work up unifrac for metadata1
generate.square.dist.matrix(unifrac,metadata1,3)
unifrac.dist.matrix1=dist.matrix1
unifrac.dist1=as.dist(unifrac.dist.matrix1)

#work up unifrac for metadata.t0
generate.square.dist.matrix(unifrac,metadata.t0,3)
unifrac.dist.matrix.t0=dist.matrix1
unifrac.dist.t0=as.dist(unifrac.dist.matrix.t0)

#work up unifrac for metadata.tend
generate.square.dist.matrix(unifrac,metadata.tend,3)
unifrac.dist.matrix.tend=dist.matrix1
unifrac.dist.tend=as.dist(unifrac.dist.matrix.tend)

#work up unifrac for metadata.coral.tend
generate.square.dist.matrix(unifrac,metadata.coral.tend,3)
unifrac.dist.matrix.coral.tend=dist.matrix1
unifrac.dist.coral.tend=as.dist(unifrac.dist.matrix.coral.tend)

#Visualize clustering dendogram of multivariate samples.
clust.coral.tend <- hclust(unifrac.dist.coral.tend,method="ward.D") #cluster tend coral samples
clust.tend <- hclust(unifrac.dist.tend,method="ward.D") #cluster tend samples

#rotate leaves of dendogram
clust.coral.tend1 <- rotate(clust.coral.tend, order=c("ABC_068","ABC_067","ABC_069","ABC_065","ABC_055_v1","ABC_057","ABC_056","ABC_060","ABC_059_v1","ABC_066"))
clust.tend1 <- rotate(clust.tend, order=c("ABC_061","ABC_062","ABC_063","ABC_070","ABC_071","ABC_068","ABC_067","ABC_069","ABC_065","ABC_055_v1","ABC_057","ABC_056","ABC_060","ABC_059_v1","ABC_066"))

#adjust labels
clust.tend1$labels <- metadata.tend$Treatment

#plot legend (seperately) for heirarchical clustering plot and export as jpeg
dev.off()
jpeg("../figures/16S_dendogram_legend.jpg",width=1500, height=2100, res=300)
plot(NULL, axes=F, bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(0, .5, legend=levels(metadata.tend$Treatment), col=c("#00BFC4","#F8766D","#00BFC4","#F8766D","black","gray"), pt.bg=c("dodgerblue3","firebrick3","NA","NA","NA","NA"), pch=21, pt.cex=2, pt.lwd=2.75, cex=.75, y.intersp=1.5, bty='n')
dev.off()

#generate outline and fill vectors for heirarchical clustering plot
fgvec <- c("#00BFC4","#F8766D","#00BFC4","#F8766D","black","gray")
bgvec <- c("dodgerblue3","firebrick3","NA","NA","NA","NA")

#plot heirarchical clustering plot, save as jpeg
jpeg("../figures/16S_dendogram.jpg",width=1500, height=2100, res=300)
plot(clust.tend1, xlab=NA, sub=NA, hang=-1,)#+
#legend("topleft",legend=levels(metadata.tend$Treatment1), col=c("#F8766D","#00BFC4","#F8766D","#00BFC4","black","gray"), pt.bg=c("firebrick3","dodgerblue3","NA","NA","NA","NA"), pch=21, pt.cex=.5, cex=.5, bty='n')
symbols(1:15, rep(0,15), circles=rep(1,15), add=T, fg=fgvec[metadata.tend$Treatment][clust.tend1$order], bg=bgvec[metadata.tend$Treatment][clust.tend1$order], inches=.09, xpd=T, lwd=3)
dev.off()

#Test statistical effects using permanova.
#first test for an effect of treatment at t0.
permanova.t0 <- adonis2(unifrac.dist.t0~Treatment, by="margin", permutations=999, data=metadata.t0)
permanova.t0 #Treatment1 is significant.

#now test tend data for all treatments.
permanova.tend=adonis2(unifrac.dist.tend~Treatment, by="margin", permutations=999, data=metadata.tend)
permanova.tend #treatment is significant, with very high R2 and medium F values.
pairwsie.adonis.tend <- pairwise.adonis(unifrac.dist.tend,factors=metadata.tend$Treatment) #no significant differences because sample size is too small

#export pairwise adonis tend for visualization in JMP
write.csv(pairwsie.adonis.tend, file.path(dirOutput, "pairwise.adonis.tend.csv"), )

#now test tend coral data for all treatments
permanova.coral.tend=adonis2(unifrac.dist.coral.tend~Treatment, by="margin", permutations=999, data=metadata.coral.tend)
permanova.coral.tend #treatment is significant.
pairwise.adonis(unifrac.dist.coral.tend,factors=metadata.coral.tend$Treatment) #no significant differences because sample size is too small

