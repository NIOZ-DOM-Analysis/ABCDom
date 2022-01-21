#Load libraries

library(ggplot2)
library(vegan)
library(pairwiseAdonis)
library(stats)
library(dendextend)
library(dplyr)
library(scales)
library(FSA)
library(phyloseq)
library(DESeq2)
library(pheatmap)

#load custom functions

source(file="generate.square.dist.script.07.27.2020.R")
source(file="generate.long.format.script.R")
source(file="cull.otu.script.R")
source(file="log.2.fold.change.calculation.05272021.R")

#Read in data/metadata.
metadata1_16S <- read.csv(file.path(dirOutput, "metadata1_16S.csv"))
abund <- read.csv(file.path(dirRAW, "16S", file="abundance_table_100.shared.csv"))
relabund <- read.csv(file.path(dirRAW, "16S", "abundance_table_100.database.csv"))

#Work up metadata.
metadata.tend <- subset(metadata1_16S, Timepoint_char=="Tend")
metadata.coral.tend <- subset(metadata.tend, Origin_PlanC!="control")

#Subset abund df to correspond to the associated metadata.
abund.coral.tend <- abund[abund$Group %in% metadata.coral.tend$Sample_Name,]

#Cull low abundance ASVs
#Format adbund df for culling.
rownames(abund.coral.tend) <- abund.coral.tend$Group #update rownames
abund.coral.tend.1 <- abund.coral.tend[,-1:-3] #remove unnecessary columns

#cull ASVs
cull.otu(abund.coral.tend.1,3,12,120) #minimum number of reads = 12 (corresponds to relabund .001) in 3 samples or 120 (relabund = .01) in 1 sample. 121 ASVs remain.
abund.coral.tend.1.cull <- relabund.df.cull #save as new df
abund.coral.tend.1.cull.t <- as.data.frame(t(abund.coral.tend.1.cull)) #transpose

#because we will be working with log 2 fold change relative to the control, we need to additionally cull ASVs that exhibit high degrees of variance in the control. To do this we will calculate SD for each ASV in the controls and normalize to the mean abund of that ASV in each control (CV).
sd.abund.coral.tend.1.cull <- as.data.frame(apply(abund.coral.tend.1.cull[8:10,], 2, FUN=sd)) #calculate stderror for each ASV in the controls, save as new df
colnames(sd.abund.coral.tend.1.cull) <- "SD"

sd.abund.coral.tend.1.cull$mean <- apply(abund.coral.tend.1.cull[8:10,], 2, FUN=mean) #calculate mean

sd.abund.coral.tend.1.cull$CV <- sd.abund.coral.tend.1.cull$SD/sd.abund.coral.tend.1.cull$mean #calculate CV

sd.abund.coral.tend.1.cull$CV[is.na(sd.abund.coral.tend.1.cull$CV)] <- 0 #manually replace NAs with 0s.

hist(sd.abund.coral.tend.1.cull$CV, breaks=seq(from=0,to=2,by=.1)) ##check distribution. for starters, .9 and 1 seem like a good threshold.

#determine appropriate CV cutoff threshold
#mean and sd of CV
mean(sd.abund.coral.tend.1.cull$CV) #.56
sd(sd.abund.coral.tend.1.cull$CV) #.37
mean(sd.abund.coral.tend.1.cull$CV) + sd(sd.abund.coral.tend.1.cull$CV) #threshold = mean CV + one SD = .93

#cull ASVS who's CV in control is above this threshold
abund.coral.tend.1.cull1 <- abund.coral.tend.1.cull[,sd.abund.coral.tend.1.cull$CV <= 0.9316542] #that leaves 103 ASVs.

abund.coral.tend.1.cull1.t <- as.data.frame(t(abund.coral.tend.1.cull1))#transpose

#Prep taxonomy data for analysis in DESEQ2
#extract taxonomy columns from relabund df.
taxonomy <- relabund[,234:236]
taxonomy1 <- cbind(taxonomy,t(as.data.frame(strsplit(as.character(taxonomy$OTUConTaxonomy),split=";")))) #split the OTU constaxonomy column and add to df.
colnames(taxonomy1)[c(1,4:9)] <- c("OTUNumber","Domain","Phylum","Class","Order","Family","Genus")

#subset taxonomy for only culled ASVs.
taxonomy1.cull=taxonomy1[taxonomy1$OTUNumber %in% colnames(abund.coral.tend.1.cull1),]
rownames(taxonomy1.cull)=taxonomy1.cull$OTUNumber #adjust rownames

#generat longformat of abund data
colnames(metadata.coral.tend)[3] <- "Sample_ID" #prep colnames
rownames(metadata.coral.tend)=metadata.coral.tend$Sample_Name #prep rownames
abund.raw.longformat <- generate.long.format(abund.coral.tend.1.cull1.t,metadata.coral.tend,taxonomy1.cull)

#Convert dfs to phyloseq object for use in DESeq2.
physeq.abund.coral.tend.1.cull1=otu_table(abund.coral.tend.1.cull1,taxa_are_rows=FALSE) #convert abund object

physeq.tax.abund.coral.tend.1.cull1=tax_table(as.matrix(taxonomy1.cull[,-1:-3])) #remove fasta, OTUNumber, and mothur taxonomy string from df, convert to matrix prior to converting to physeq object.

rownames(metadata.coral.tend)=metadata.coral.tend$Sample_ID #convert rownames to sample names
physeq.metadata.coral.tend=sample_data(metadata.coral.tend) #convert to physeq object

physeq.coral.tend.cull1=phyloseq(physeq.abund.coral.tend.1.cull1,physeq.tax.abund.coral.tend.1.cull1,physeq.metadata.coral.tend) #combine

#Run DESeq2
mod.deseq2 <- phyloseq_to_deseq2(physeq.coral.tend.cull1, ~ Treatment)

mod.deseq2 <- estimateSizeFactors(mod.deseq2,type="poscounts")

mod.deseq2 <- estimateDispersions(mod.deseq2)

mod.deseq2 <- nbinomWaldTest(mod.deseq2)

mod.deseq2.heated <- as.data.frame(results(mod.deseq2,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Non-bleached + Ambient", "Non-bleached + Heated")))
mod.deseq2.bleached <- as.data.frame(results(mod.deseq2,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Non-bleached + Ambient", "Bleached + Ambient")))
mod.deseq2.bleach.thermal=as.data.frame(results(mod.deseq2,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment","Non-bleached + Ambient", "Bleached + Heated")))

#check distributions of pvals
hist(mod.deseq2.heated$padj, breaks=seq(from=0,to=1,by=.05)) #2 sig ASVs
hist(mod.deseq2.bleached$padj, breaks=seq(from=0,to=1,by=.05)) #5 sig ASVs
hist(mod.deseq2.bleach.thermal$padj, breaks=seq(from=0,to=1,by=.05)) #0 sig ASV
