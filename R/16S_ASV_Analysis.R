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
library(cowplot)
library(ggforce)
library(grid)
library(gridExtra)
library(stringr)
library(cowplot)


#load custom functions

source(file="generate.square.dist.script.07.27.2020.R")
source(file="generate.long.format.script.R")
source(file="cull.otu.script.R")
source(file="log.2.fold.change.calculation.05272021.R")

#Read in data/metadata.
metadata1_16S <- read.csv(file.path(dirOutput, "metadata1_16S.csv"))
abund <- read.csv(file.path(dirRAW, "16S", file="abundance_table_100.shared.csv"))
abund.nosub <- read.csv(file.path(dirRAW, "16S", file="all_multipletonsFilter_100.count_table.csv")) #read in non-subsampled data
relabund <- read.csv(file.path(dirRAW, "16S", "abundance_table_100.database.csv"))
nosub.asv.list <- read.csv(file.path(dirRAW, "16S", "all_multipletonsFilter_100.list.csv")) #read in list converting esvs to OTUs

#Work up metadata.
metadata.tend <- subset(metadata1_16S, Timepoint_char=="Tend")
metadata.coral.tend <- subset(metadata.tend, Origin_PlanC!="control")

#Subset abund df to correspond to the associated metadata.
abund.coral.tend <- abund[abund$Group %in% metadata.coral.tend$Sample_Name_Unique,]

#RE DO CULLING WITH ABUND.NOSUB ASVS.
abund.nosub1 <- merge(abund.nosub, nosub.asv.list, by.x="Representative_Sequence", by.y="ESV") #merge abund.nosub and the nosub.asv.list
abund.nosub2 <- abund.nosub1[,-1:-2] #remove unnecessary columns

#Next, workup abund.nosub.cull so that it only contains the relevant samples.
rownames(abund.nosub2) <- abund.nosub2$OTU #update rownames
abund.nosub.coral.tend <- abund.nosub2[,colnames(abund.nosub2) %in% metadata.coral.tend$Sample_ID] #subset for relevant samples.

#Cull ow abundance ASVs
hist(unlist(abund.nosub.coral.tend)[unlist(abund.nosub.coral.tend)!=0], xlim=c(0,5000), breaks=10000) #visualize distribution of abundance values. After removing zeros, it looks like most ASVs have an abundance < 1000.
abund.nosub.coral.tend.t <- as.data.frame(t(abund.nosub.coral.tend)) #transpose

#cull ASVs
cull.otu(abund.nosub.coral.tend.t,3,20,1000) #minimum number of reads = 50 in 3 samples or 1000 in 1 sample. 301 ASVs remain.
abund.nosub.coral.tend.t.cull <- relabund.df.cull #save as new df
abund.nosub.coral.tend.cull <- as.data.frame(t(abund.nosub.coral.tend.t.cull)) #transpose







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
rownames(metadata.coral.tend) <- metadata.coral.tend$Sample_ID #prep rownames
abund.raw.longformat <- generate.long.format(abund.coral.tend.1.cull1.t,metadata.coral.tend,taxonomy1.cull)

#Convert dfs to phyloseq object for use in DESeq2.
physeq.abund.coral.tend.1.cull1=otu_table(abund.coral.tend.1.cull1[,-101],taxa_are_rows=FALSE) #convert abund object, be sure to remove treatment column

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

#extract l2fc values and add to sig.asvs
sig.asvs$heated.l2fc <- mod.deseq2.heated$log2FoldChange
sig.asvs$bleached.l2fc <- mod.deseq2.bleached$log2FoldChange
sig.asvs$bleached.heated.l2fc <- mod.deseq2.bleach.heated$log2FoldChange

#add new columns indicating significance Y/N
sig.asvs$heated.sig <- sig.asvs$heated.padj #duplicate padj
sig.asvs$heated.sig[sig.asvs$heated.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs$heated.sig[sig.asvs$heated.padj <=.05] <- "Y"
sig.asvs$bleached.sig <- sig.asvs$bleached.padj #duplicate padj
sig.asvs$bleached.sig[sig.asvs$bleached.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs$bleached.sig[sig.asvs$bleached.padj <=.05] <- "Y"
sig.asvs$bleached.heated.sig <- sig.asvs$bleached.heated.padj #duplicate padj
sig.asvs$bleached.heated.sig[sig.asvs$bleached.heated.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs$bleached.heated.sig[sig.asvs$bleached.heated.padj <=.05] <- "Y"

#merge in taxonomy
sig.asvs1 <- cbind(sig.asvs, taxonomy1.cull)

#add in mean l2fc values for each ASV in each treatment
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

abund.longformat.merged$Treatment_OTU <- paste(abund.longformat.merged$Treatment.x, abund.longformat.merged$OTU.x, sep="_") #add a treatment_OTU column

#calculate mean relabund and l2fc for each treatment_OTU
mean.abund.l2fc <- as.data.frame(aggregate(abund.longformat.merged[,4], by=list(abund.longformat.merged$Treatment_OTU), FUN=mean)) #calculate mean abund
mean.abund.l2fc1 <- cbind(mean.abund.l2fc, as.data.frame(aggregate(abund.longformat.merged[,36], by=list(abund.longformat.merged$Treatment_OTU), FUN=mean))) #calculate mean l2fc, add to df
mean.abund.l2fc2 <- mean.abund.l2fc1[,-3] #remove unnecessary columns
colnames(mean.abund.l2fc2) <- c("Treatment_OTU", "Mean Abundance", "Mean Log2 Fold Change") #rename columns
mean.abund.l2fc3 <- cbind(mean.abund.l2fc2, t(as.data.frame(strsplit(mean.abund.l2fc2$Treatment_OTU, split="_")))) #split Treatment_OTU
colnames(mean.abund.l2fc3)[4:5] <- c("Treatment", "OTU") #update colnames

#combine mean.abund.l2fc3 values with sig.asvs1
#first for abundance
sig.asvs1$mean.abund.heated <- mean.abund.l2fc3$`Mean Abundance`[mean.abund.l2fc3$Treatment == "Non-bleached + Heated"]
sig.asvs1$mean.abund.bleached <- mean.abund.l2fc3$`Mean Abundance`[mean.abund.l2fc3$Treatment == "Bleached + Ambient"]
sig.asvs1$mean.abund.bleached.heated <- mean.abund.l2fc3$`Mean Abundance`[mean.abund.l2fc3$Treatment == "Bleached + Heated"]
sig.asvs1$mean.abund.ambient <- mean.abund.l2fc3$`Mean Abundance`[mean.abund.l2fc3$Treatment == "Non-bleached + Ambient"]

#then for l2fc
sig.asvs1$mean.l2fc.heated <- mean.abund.l2fc3$`Mean Log2 Fold Change`[mean.abund.l2fc3$Treatment == "Non-bleached + Heated"]
sig.asvs1$mean.l2fc.bleached <- mean.abund.l2fc3$`Mean Log2 Fold Change`[mean.abund.l2fc3$Treatment == "Bleached + Ambient"]
sig.asvs1$mean.l2fc.bleached.heated <- mean.abund.l2fc3$`Mean Log2 Fold Change`[mean.abund.l2fc3$Treatment == "Bleached + Heated"]

#compare the l2fc means that I calculated from the ones DESeq2 calcualted.
#First multiply the deseq2 l2fc values by -1 to get the correct sign for comparison.
sig.asvs1$heated.l2fc.inv <- -1*sig.asvs1$heated.l2fc
sig.asvs1$bleached.l2fc.inv <- -1*sig.asvs1$bleached.l2fc
sig.asvs1$bleached.heated.l2fc.inv <- -1*sig.asvs1$bleached.heated.l2fc

#calculate the difference.
sig.asvs1$l2fc.diff.heated <- sig.asvs1$mean.l2fc.heated - sig.asvs1$heated.l2fc.inv
sig.asvs1$l2fc.diff.bleached <- sig.asvs1$mean.l2fc.bleached - sig.asvs1$bleached.l2fc.inv
sig.asvs1$l2fc.diff.bleached.heated <- sig.asvs1$mean.l2fc.bleached.heated - sig.asvs1$bleached.heated.l2fc.inv

#plot a comparison of deseq2 vs. personally calculated l2fc
bleached.l2fc.comparison.plot <- ggplot(sig.asvs1, aes(x=l2fc.diff.bleached, y=ASV, fill=mean.l2fc.bleached))+
  geom_bar(stat="identity", color="black")+
  scale_fill_gradient2(low="blue",mid="white",high="red")+
  xlab("Caclulated l2fc - DEseq2 l2fc")+
  ggtitle("Bleached")+
  theme(axis.text.y=element_blank(), axis.title.y=element_blank())

heated.l2fc.comparison.plot <- ggplot(sig.asvs1, aes(x=l2fc.diff.heated, y=ASV, fill=mean.l2fc.heated))+
  geom_bar(stat="identity", color="black")+
  scale_fill_gradient2(low="blue",mid="white",high="red")+
  xlab("Caclulated l2fc - DEseq2 l2fc")+
  ggtitle("Heated")+
  theme(axis.text.y=element_blank(), axis.title.y=element_blank())

bleached.heated.l2fc.comparison.plot <- ggplot(sig.asvs1, aes(x=l2fc.diff.bleached.heated, y=ASV, fill=mean.l2fc.bleached.heated))+
  geom_bar(stat="identity", color="black")+
  scale_fill_gradient2(low="blue",mid="white",high="red")+
  xlab("Caclulated l2fc - DEseq2 l2fc")+
  ggtitle("Bleached+Heated")

png(filename="../figures/ASV l2fc comparison.png", width=10000, height=10000, res=600)
plot_grid(bleached.heated.l2fc.comparison.plot, heated.l2fc.comparison.plot, bleached.l2fc.comparison.plot, nrow=1, rel_widths = c(1.3,1,1))
dev.off()

#Next, try DESEq2 on raw, unsubsampled, un-lulu'ed data.

#First, work up abund.nosub so it only contains the ASVs from the prior, culled, abund df.


abund.nosub.cull <- abund.nosub2[abund.nosub2$OTU %in% taxonomy1.cull$OTUNumber,] #subset rows for only relevant OTUs


#work up in phyloseq for DESEq2
abund.nosub.cull1.t <- as.data.frame(t(abund.nosub.cull1)) #transpose
physeq.abund.nosub.cull1=otu_table(abund.nosub.cull1.t,taxa_are_rows=FALSE) #convert abund object

physeq.coral.nosub.cull1=phyloseq(physeq.abund.nosub.cull1,physeq.tax.abund.coral.tend.1.cull1,physeq.metadata.coral.tend) #combine

#Run DESEq2 on this data.
mod.deseq3 <- phyloseq_to_deseq2(physeq.coral.nosub.cull1,~ Treatment)

mod.deseq3 <- estimateSizeFactors(mod.deseq3,type="poscounts")

mod.deseq3 <- estimateDispersions(mod.deseq3)

mod.deseq3 <- nbinomWaldTest(mod.deseq3)

mod.deseq3.heated <- as.data.frame(results(mod.deseq3,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Non-bleached + Ambient", "Non-bleached + Heated")))
mod.deseq3.bleached <- as.data.frame(results(mod.deseq3,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Non-bleached + Ambient", "Bleached + Ambient")))
mod.deseq3.bleach.heated <- as.data.frame(results(mod.deseq3,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment","Non-bleached + Ambient", "Bleached + Heated")))

#visualize pvalues
hist(mod.deseq3.heated$padj, breaks=seq(from=0,to=1,by=.05)) #6 sig ASVs
hist(mod.deseq3.bleached$padj, breaks=seq(from=0,to=1,by=.05)) #13? sig ASVs
hist(mod.deseq3.bleach.heated$padj, breaks=seq(from=0,to=1,by=.05)) #5 sig ASVs

#extract pvalues and add to new df
sig.asvs2 <- as.data.frame(rownames(as.data.frame(mod.deseq3.bleached))) #add asv names
colnames(sig.asvs2) <- "ASV" #change colnames
sig.asvs2$heated.p <- mod.deseq3.heated$pvalue
sig.asvs2$heated.padj <- mod.deseq3.heated$padj
sig.asvs2$bleached.p <- mod.deseq3.bleached$pvalue
sig.asvs2$bleached.padj <- mod.deseq3.bleached$padj
sig.asvs2$bleached.heated.p <- mod.deseq3.bleach.heated$pvalue
sig.asvs2$bleached.heated.padj <- mod.deseq3.bleach.heated$padj

#extract l2fc values and add to sig.asvs2
sig.asvs2$heated.l2fc <- mod.deseq3.heated$log2FoldChange
sig.asvs2$bleached.l2fc <- mod.deseq3.bleached$log2FoldChange
sig.asvs2$bleached.heated.l2fc <- mod.deseq3.bleach.heated$log2FoldChange

#add new columns indicating significance Y/N
sig.asvs2$heated.sig <- sig.asvs2$heated.padj #duplicate padj
sig.asvs2$heated.sig[sig.asvs2$heated.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs2$heated.sig[sig.asvs2$heated.padj <=.05] <- "Y"
sig.asvs2$bleached.sig <- sig.asvs2$bleached.padj #duplicate padj
sig.asvs2$bleached.sig[sig.asvs2$bleached.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs2$bleached.sig[sig.asvs2$bleached.padj <=.05] <- "Y"
sig.asvs2$bleached.heated.sig <- sig.asvs2$bleached.heated.padj #duplicate padj
sig.asvs2$bleached.heated.sig[sig.asvs2$bleached.heated.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs2$bleached.heated.sig[sig.asvs2$bleached.heated.padj <=.05] <- "Y"

#export
#write.csv(sig.asvs2, file.path(dirOutput, "sig.asvs2.csv"))

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

#remove control samples for visualization
abund.longformat.merged1 <- subset(abund.longformat.merged,Sample.x!="ABC_067" & Sample.x!="ABC_068" & Sample.x!="ABC_069")

#pad Genus_OTU strings to 42 characters. Need to do this for alignment of faceted plots to work
abund.longformat.merged1$Genus_OTU1 <- str_pad(abund.longformat.merged1$Genus_OTU, width=42, "left")

#make all lowercase and then pad to see if this resolves the alignment problem
abund.longformat.merged1$Genus_OTU2 <- tolower(abund.longformat.merged1$Genus_OTU)
abund.longformat.merged1$Genus_OTU3 <- str_pad(abund.longformat.merged1$Genus_OTU2, width=42, "left", pad="_")

#Visualize
ggplot(abund.longformat.merged1,aes(y=Genus_OTU,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
  geom_point()+
  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.1375,.55,.6625,1))+
  theme(axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
  xlab(label="Sample")+
  labs(size="Abundance")
#theme_bw()
ggsave('ASV_l2fc.jpg', path=dirFigs, width=6, height=16, dpi = 1200)

#Visualize, faceting by Order
ggplot(abund.longformat.merged1,aes(y=Genus_OTU,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
  geom_point()+
  facet_wrap(.~Order.x,scales="free")+
  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.1375,.55,.6625,1))+
  theme(axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
  xlab(label="Sample")+
  labs(size="Abundance")
#theme_bw()
ggsave('ASV_l2fc_order.jpg', path=dirFigs, width=28, height=26, dpi = 600)

#Visualize, faceting by class
ggplot(abund.longformat.merged1,aes(y=Genus_OTU,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
  geom_point()+
  facet_grid(row=vars(Class.x), scales="free_y", space="fre")+
  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.1375,.55,.6625,1))+
  theme(axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
  xlab(label="Sample")+
  labs(size="Abundance")
#theme_bw()
ggsave('ASV_l2fc_class.jpg', path=dirFigs, width=20, height=26, dpi = 600)

#for class, facet out the abundant classes (aprot, gprot, bacteroi) and the not abundantclasses (oxy, prot_unclass, dprot, bact_unclass)
bubbleplot_class1 <- ggplot(subset(abund.longformat.merged1, Class.x=="Alphaproteobacteria" | Class.x=="Gammaproteobacteria" | Class.x=="Bacteroidia"),aes(y=Genus_OTU1,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
  geom_point()+
  facet_wrap(.~Class.x,scales="free_y")+
  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.325,.5,.675,1), limits=c(-10,10))+
  theme(axis.text.y=element_text(family="mono"), axis.text.x=element_blank(), axis.title.x=element_blank())+
  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
  labs(size="Abundance")

bubbleplot_class2 <- ggplot(subset(abund.longformat.merged1, Class.x=="Deltaproteobacteria" | Class.x=="Oxyphotobacteria" | Class.x=="Proteobacteria_unclassified" | Class.x=="Bacteria_unclassified"),aes(y=Genus_OTU1,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
  geom_point()+
  facet_wrap(.~Class.x,scales="free_y", nrow=2, ncol=3)+
  scale_size_continuous(limits=c(0,5000))+
  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.325,.5,.675,1), limits=c(-10,10))+
  theme(legend.position="none",axis.text.y=element_text(family="mono"),axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
  xlab(label="Sample")+
  labs(size="Abundance")

#combine plots with plot.grid. PROBLEMS WITH ALIGNMENT CUZ OF YAXIS LABELS NEED TO PAD STRINGS
plot_grid(bubbleplot_class1, bubbleplot_class2, align="v", ncol=1, axis="lr", rel_heights = c(1,.7))
#8x17 landscape

#different method
bubbleplot_class3 <- ggplot(subset(abund.longformat.merged1, Class.x=="Alphaproteobacteria" | Class.x=="Gammaproteobacteria" | Class.x=="Bacteroidia"),aes(y=Genus_OTU,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
  geom_point()+
  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.325,.5,.675,1), limits=c(-10,10))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())+
  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
  labs(size="Abundance")+
  facet_col(vars(Class.x), scales="free_y", space="free")

bubbleplot_class3_v1 <- ggplot(subset(abund.longformat.merged1, Class.x=="Alphaproteobacteria" | Class.x=="Gammaproteobacteria" | Class.x=="Bacteroidia"),aes(y=Genus_OTU,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
  geom_point()+
  scale_size_continuous(limits=c(0,5000))+
  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.325,.5,.675,1), limits=c(-10,10))+
  theme(legend.position="none",axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
  xlab(label="Sample")+
  labs(size="Abundance")+
  facet_col(vars(Class.x), scales="free_y", space="free")

bubbleplot_class4 <- ggplot(subset(abund.longformat.merged1, Class.x=="Deltaproteobacteria" | Class.x=="Oxyphotobacteria" | Class.x=="Proteobacteria_unclassified" | Class.x=="Bacteria_unclassified"),aes(y=Genus_OTU,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
  geom_point()+
  scale_size_continuous(limits=c(0,5000))+
  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.325,.5,.675,1), limits=c(-10,10))+
  theme(legend.position="none",axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
  xlab(label="Sample")+
  labs(size="Abundance")+
  facet_col(vars(Class.x), scales="free_y", space="free")

#plotgrid
plot_grid(bubbleplot_class3, bubbleplot_class4, align="v", ncol=1, axis="lr", rel_heights = c(1,.25))

#Next, visualize boxplots of sig ASVs from mod.deseq2 and mod.deseq3
sig.asvs3 <- sig.asvs2[is.na(sig.asvs2$heated.p)==FALSE,] #remove NAs

abund.longformat.merged.sig <- abund.longformat.merged[abund.longformat.merged$OTU.x %in% sig.asvs3$ASV[sig.asvs3$heated.sig == "Y" | sig.asvs3$bleached.sig == "Y" | sig.asvs3$bleached.heated.sig == "Y"],] #subset abund.longformat.merged1 for just mod.deseq3 sig asvs.
abund.longformat.merged.sig$Treatment.x <- factor(abund.longformat.merged.sig$Treatment.x, levels=levels(fact.all.treat)) #adjust factor levels

ggplot(abund.longformat.merged.sig, aes(y=abund, x=Treatment.x, color=Treatment.x, fill=Treatment.x))+
  geom_boxplot()+
  facet_wrap(.~Genus_OTU, scales="free")+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
ggsave('Sig ASVs mod.deseq3.jpeg', path = dirFigs, width = 24, height = 18, dpi = 600)
