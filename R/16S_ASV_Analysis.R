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
abund.coral.tend <- abund[abund$Group %in% metadata.coral.tend$Sample_Name_Unique,]

#Cull low abundance ASVs
#Format adbund df for culling.
rownames(abund.coral.tend) <- abund.coral.tend$Group #update rownames
abund.coral.tend.1 <- abund.coral.tend[,-1:-3] #remove unnecessary columns

#cull ASVs
cull.otu(abund.coral.tend.1,3,12,120) #minimum number of reads = 12 (corresponds to relabund .001) in 3 samples or 120 (relabund = .01) in 1 sample. 117 ASVs remain.
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
mean(sd.abund.coral.tend.1.cull$CV) #.54
sd(sd.abund.coral.tend.1.cull$CV) #.36
mean(sd.abund.coral.tend.1.cull$CV) + sd(sd.abund.coral.tend.1.cull$CV) #threshold = mean CV + one SD = .9

#cull ASVS who's CV in control is above this threshold
abund.coral.tend.1.cull1 <- abund.coral.tend.1.cull[,sd.abund.coral.tend.1.cull$CV <= 0.9030586] #that leaves 100 ASVs.

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
colnames(metadata.coral.tend)[4] <- "Sample_ID" #prep colnames
rownames(metadata.coral.tend)=metadata.coral.tend$Sample_ID #prep rownames
abund.raw.longformat <- generate.long.format(abund.coral.tend.1.cull1.t,metadata.coral.tend,taxonomy1.cull)

#Convert dfs to phyloseq object for use in DESeq2.
physeq.abund.coral.tend.1.cull1=otu_table(abund.coral.tend.1.cull1,taxa_are_rows=FALSE) #convert abund object

physeq.tax.abund.coral.tend.1.cull1=tax_table(as.matrix(taxonomy1.cull[,-1:-3])) #remove fasta, OTUNumber, and mothur taxonomy string from df, convert to matrix prior to converting to physeq object.

physeq.metadata.coral.tend=sample_data(metadata.coral.tend) #convert to physeq object

physeq.coral.tend.cull1=phyloseq(physeq.abund.coral.tend.1.cull1,physeq.tax.abund.coral.tend.1.cull1,physeq.metadata.coral.tend) #combine

#Run DESeq2
mod.deseq2 <- phyloseq_to_deseq2(physeq.coral.tend.cull1,~ Treatment)

mod.deseq2 <- estimateSizeFactors(mod.deseq2,type="poscounts")

mod.deseq2 <- estimateDispersions(mod.deseq2)

mod.deseq2 <- nbinomWaldTest(mod.deseq2)

mod.deseq2.heated <- as.data.frame(results(mod.deseq2,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Non-bleached + Ambient", "Non-bleached + Heated")))
mod.deseq2.bleached <- as.data.frame(results(mod.deseq2,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Non-bleached + Ambient", "Bleached + Ambient")))
mod.deseq2.bleach.heated=as.data.frame(results(mod.deseq2,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment","Non-bleached + Ambient", "Bleached + Heated")))

#check distributions of pvals
hist(mod.deseq2.heated$padj, breaks=seq(from=0,to=1,by=.05)) #4 sig ASVs
hist(mod.deseq2.bleached$padj, breaks=seq(from=0,to=1,by=.05)) #many sig ASVs
hist(mod.deseq2.bleach.heated$padj, breaks=seq(from=0,to=1,by=.05)) #1 sig ASV

#extract pvalues and add to new df
sig.asvs <- as.data.frame(rownames(as.data.frame(mod.deseq2.bleached))) #add asv names
colnames(sig.asvs) <- "ASV" #change colnames
sig.asvs$heated.p <- mod.deseq2.heated$pvalue
sig.asvs$heated.padj <- mod.deseq2.heated$padj
sig.asvs$bleached.p <- mod.deseq2.bleached$pvalue
sig.asvs$bleached.padj <- mod.deseq2.bleached$padj
sig.asvs$bleached.heated.p <- mod.deseq2.bleach.heated$pvalue
sig.asvs$bleached.heated.padj <- mod.deseq2.bleach.heated$padj

#add new columns indicating significance Y/N
sig.asvs$heated.sig <- sig.asvs$heated.padj #duplicate padj
sig.asvs$heated.sig[sig.asvs$heated.sig >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs$heated.sig[sig.asvs$heated.sig <=.05] <- "Y"
sig.asvs$bleached.sig <- sig.asvs$bleached.padj #duplicate padj
sig.asvs$bleached.sig[sig.asvs$bleached.sig >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs$bleached.sig[sig.asvs$bleached.sig <=.05] <- "Y"
sig.asvs$bleached.heated.sig <- sig.asvs$bleached.heated.padj #duplicate padj
sig.asvs$bleached.heated.sig[sig.asvs$bleached.heated.sig >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs$bleached.heated.sig[sig.asvs$bleached.heated.sig <=.05] <- "Y"

#Manually calculate log2fold change and visualize.
#workup data for log2foldchange function.
abund.coral.tend.1.cull1$Treatment <- metadata.coral.tend$Treatment #add treatment column.

log2.fold.change(abund.coral.tend.1.cull1,101,"Non-bleached + Ambient",100) #run the function

hist(unlist(log2.fold.df1[,-101]))#check the distribution of log2foldchange data.

#Generate longformat of log2FC data.
log2.fold.df1.t=as.data.frame(t(log2.fold.df1[,-101])) #work up the data

abund.lfc.longformat=generate.long.format(log2.fold.df1.t,metadata.coral.tend,taxonomy1.cull) #generate longformat
colnames(abund.lfc.longformat)[3]="Log2FoldChange" #rename column

#Merge the two abund/lfc longformat dfs.
#work up the dfs.
abund.raw.longformat$Sample_OTU=paste(abund.raw.longformat$Sample,abund.raw.longformat$OTU)
abund.lfc.longformat$Sample_OTU=paste(abund.lfc.longformat$Sample,abund.lfc.longformat$OTU)

abund.longformat.merged=merge(abund.raw.longformat,abund.lfc.longformat,by.x="Sample_OTU",by.y="Sample_OTU") #merge

#Visualize
#First, perform heirarchical clustering on the l2fc values
lfc.clustering=pheatmap(log2.fold.df1[-8:-10,-101])

lfc.clustering1=pheatmap(log2.fold.df1[-8:-10,-101],cluster_rows = FALSE)

cluster.OTUs=lfc.clustering$tree_col$labels[lfc.clustering$tree_col$order]
cluster.samples=lfc.clustering$tree_row$labels[lfc.clustering$tree_row$order]
cluster.Genus=taxonomy1.cull$Genus[lfc.clustering$tree_col$order]
cluster.Genus_OTUs=paste(cluster.Genus,cluster.OTUs,sep="_")

#reorder OTUs in abund.longformat.merged
abund.longformat.merged$OTU.x=factor(abund.longformat.merged$OTU.x,levels=cluster.OTUs)

#add a new column correspondong to family_OTU
abund.longformat.merged$Genus_OTU=paste(abund.longformat.merged$Genus.x,abund.longformat.merged$OTU.x,sep="_")
abund.longformat.merged$Genus_OTU=factor(abund.longformat.merged$Genus_OTU,levels=cluster.Genus_OTUs)

#Visualize
ggplot(subset(abund.longformat.merged,Sample.x!="ABC_067" & Sample.x!="ABC_068" & Sample.x!="ABC_069"),aes(y=Genus_OTU,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
  geom_point()+
  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.1375,.55,.6625,1))+
  theme(axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
  xlab(label="Sample")+
  labs(size="Abundance")
#theme_bw()
ggsave('ASV_l2fc.jpg', path=dirFigs, width=6, height=16, dpi = 1200)

#Visualize, faceting by Family
ggplot(subset(abund.longformat.merged,Sample.x!="ABC_067" & Sample.x!="ABC_068" & Sample.x!="ABC_069"),aes(y=Genus_OTU,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
  geom_point()+
  facet_wrap(.~Family.x,scales="free")+
  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.1375,.55,.6625,1))+
  theme(axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
  xlab(label="Sample")+
  labs(size="Abundance")
#theme_bw()

ggsave('ASV_l2fc_family.jpg', path=dirFigs, width=28, height=26, dpi = 600)
