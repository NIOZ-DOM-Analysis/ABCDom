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
library(MetBrewer)

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
nosub.taxonomy <- read.csv(file.path(dirRAW, "16s", "all_multipletonsFilter_100.taxonomy.csv")) #read in taxonomy of nosub asvs

#Work up metadata.
metadata.tend <- subset(metadata1_16S, Timepoint_char=="Tend")
metadata.coral.tend <- subset(metadata.tend, Origin_PlanC!="control")

#Subset abund df to correspond to the associated metadata.
abund.coral.tend <- abund[abund$Group %in% metadata.coral.tend$Sample_Name_Unique,]

#Cull abund.nosub asvs
nosub.asv.info <- merge(nosub.asv.list, nosub.taxonomy, by.x="ESV", by.y="ESV") #merge asv list and taxonomy
abund.nosub1 <- merge(abund.nosub, nosub.asv.info, by.x="Representative_Sequence", by.y="ESV") #merge abund.nosub and the nosub.asv.list
abund.nosub2 <- abund.nosub1[,-1:-2] #remove unnecessary columns

#Next, workup abund.nosub.cull so that it only contains the relevant samples.
rownames(abund.nosub2) <- abund.nosub2$OTU #update rownames
abund.nosub.coral.tend <- abund.nosub2[,colnames(abund.nosub2) %in% metadata.coral.tend$Sample_Name_Unique] #subset for relevant samples.

#Next, convert to relabund df for later boxplots
#create relabund function
relabund.calculate <- function(x) { #where x is a vector of numbers

  relabund.vec <- c() #create blank vector

  for (i in 1:length(x)) { #for each position (numeric value) in vector x
    relabund.vec[i] <- x[i]/sum(x) #take it's value and divide it by the total sum of vector x, and put that number in the corresponding position of the new vector.
  }

  relabund.vec1 <<- relabund.vec #store output as new vector
}

relabund.nosub.coral.tend <- as.data.frame(apply(abund.nosub.coral.tend, 2, FUN=relabund.calculate)) #calculate relabund
rownames(relabund.nosub.coral.tend) <- rownames(abund.nosub.coral.tend) #update rownames

#Cull low abundance ASVs
hist(unlist(abund.nosub.coral.tend)[unlist(abund.nosub.coral.tend)!=0], xlim=c(0,5000), breaks=10000) #visualize distribution of abundance values. After removing zeros, it looks like most ASVs have an abundance < 1000.
abund.nosub.coral.tend.t <- as.data.frame(t(abund.nosub.coral.tend)) #transpose. I checked and there are a total of 1068 ASVs in this set (the rest=0 in all samples).

#cull ASVs
cull.otu(abund.nosub.coral.tend.t,3,50,1000) #minimum number of reads = 50 in 3 samples or 1000 in 1 sample. 187 ASVs remain.
abund.nosub.coral.tend.t.cull <- relabund.df.cull #save as new df
abund.nosub.coral.tend.cull <- as.data.frame(t(abund.nosub.coral.tend.t.cull)) #transpose

#because we will be working with log 2 fold change relative to the control, we need to additionally cull ASVs that exhibit high degrees of variance in the control. To do this we will calculate SD for each ASV in the controls and normalize to the mean abund of that ASV in each control (CV).
sd.abund.nosub.coral.tend.t.cull <- as.data.frame(apply(abund.nosub.coral.tend.t.cull[8:10,], 2, FUN=sd)) #calculate stderror for each ASV in the controls, save as new df
colnames(sd.abund.nosub.coral.tend.t.cull) <- "SD"

sd.abund.nosub.coral.tend.t.cull$mean <- apply(abund.nosub.coral.tend.t.cull[8:10,], 2, FUN=mean) #calculate mean

sd.abund.nosub.coral.tend.t.cull$CV <- sd.abund.nosub.coral.tend.t.cull$SD/sd.abund.nosub.coral.tend.t.cull$mean #calculate CV

sd.abund.nosub.coral.tend.t.cull$CV[is.na(sd.abund.nosub.coral.tend.t.cull$CV)] <- 0 #manually replace NAs with 0s.

hist(sd.abund.nosub.coral.tend.t.cull$CV, breaks=seq(from=0,to=2,by=.1)) ##check distribution. for starters, .8, .9, and 1 seem like a good threshold.

#determine appropriate CV cutoff threshold
#mean and sd of CV
mean(sd.abund.nosub.coral.tend.t.cull$CV) #.505
sd(sd.abund.nosub.coral.tend.t.cull$CV) #.40
mean(sd.abund.nosub.coral.tend.t.cull$CV) + sd(sd.abund.nosub.coral.tend.t.cull$CV) #threshold = mean CV + one SD =  0.9111981

#cull ASVS who's CV in control is above this threshold
abund.nosub.coral.tend.t.cull1 <- abund.nosub.coral.tend.t.cull[,sd.abund.nosub.coral.tend.t.cull$CV <=  0.9111981] #that leaves 159
abund.nosub.coral.tend.cull1 <- as.data.frame(t(abund.nosub.coral.tend.t.cull1)) #transpose

#work up data for deseq2

#Prep taxonomy data for analysis in DESEQ2
#extract taxonomy data
taxonomy <- nosub.asv.info

#subset taxonomy for only culled ASVs.
taxonomy.cull=taxonomy[taxonomy$OTU %in% colnames(abund.nosub.coral.tend.t.cull1),]
rownames(taxonomy.cull)=taxonomy.cull$OTU #adjust rownames
colnames(taxonomy.cull)[2] <- "OTUNumber" #adjust colnames

#generat longformat of abund data
colnames(metadata.coral.tend)[4] <- "Sample_ID" #prep colnames
rownames(metadata.coral.tend) <- metadata.coral.tend$Sample_ID #prep rownames
abund.nosub.longformat <- generate.long.format(abund.nosub.coral.tend.cull1,metadata.coral.tend,taxonomy.cull)
relabund.nosub.longformat <- generate.long.format(relabund.nosub.coral.tend,metadata.coral.tend,taxonomy.cull) #generate longformat for relabund data as well.

#Convert dfs to phyloseq object for use in DESeq2.
physeq.abund.nosub <- otu_table(abund.nosub.coral.tend.cull1,taxa_are_rows=TRUE) #convert abund object, be sure to remove treatment column

physeq.tax.abund.nosub <- tax_table(as.matrix(taxonomy.cull)) #remove fasta, OTUNumber, and mothur taxonomy string from df, convert to matrix prior to converting to physeq object.

physeq.metadata.coral.tend <- sample_data(metadata.coral.tend) #convert to physeq object

physeq.nosub <- phyloseq(physeq.abund.nosub, physeq.tax.abund.nosub, physeq.metadata.coral.tend) #combine

#Run DESEq2 on this data.
mod.deseq4 <- phyloseq_to_deseq2(physeq.nosub,~ Treatment)

mod.deseq4 <- estimateSizeFactors(mod.deseq4,type="poscounts")

mod.deseq4 <- estimateDispersions(mod.deseq4)

mod.deseq4 <- nbinomWaldTest(mod.deseq4)

mod.deseq4.heated <- as.data.frame(results(mod.deseq4,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Non-bleached + Heated", "Non-bleached + Ambient")))
mod.deseq4.bleached <- as.data.frame(results(mod.deseq4,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Bleached + Ambient", "Non-bleached + Ambient")))
mod.deseq4.bleach.heated <- as.data.frame(results(mod.deseq4,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Bleached + Heated","Non-bleached + Ambient")))

#visualize pvalues
hist(mod.deseq4.heated$padj, breaks=seq(from=0,to=1,by=.05)) #12 sig ASVs
hist(mod.deseq4.bleached$padj, breaks=seq(from=0,to=1,by=.05)) #21 sig ASVs
hist(mod.deseq4.bleach.heated$padj, breaks=seq(from=0,to=1,by=.05)) #19 sig ASVs

#extract pvalues and add to new df
sig.asvs.v1 <- as.data.frame(rownames(as.data.frame(mod.deseq4.bleached))) #add asv names
colnames(sig.asvs.v1) <- "ASV" #change colnames
sig.asvs.v1$heated.p <- mod.deseq4.heated$pvalue
sig.asvs.v1$heated.padj <- mod.deseq4.heated$padj
sig.asvs.v1$bleached.p <- mod.deseq4.bleached$pvalue
sig.asvs.v1$bleached.padj <- mod.deseq4.bleached$padj
sig.asvs.v1$bleached.heated.p <- mod.deseq4.bleach.heated$pvalue
sig.asvs.v1$bleached.heated.padj <- mod.deseq4.bleach.heated$padj

#extract l2fc values and add to sig.asvs.v1
sig.asvs.v1$heated.l2fc <- mod.deseq4.heated$log2FoldChange
sig.asvs.v1$bleached.l2fc <- mod.deseq4.bleached$log2FoldChange
sig.asvs.v1$bleached.heated.l2fc <- mod.deseq4.bleach.heated$log2FoldChange

#add new columns indicating significance Y/N
sig.asvs.v1$heated.sig <- sig.asvs.v1$heated.padj #duplicate padj
sig.asvs.v1$heated.sig[sig.asvs.v1$heated.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs.v1$heated.sig[sig.asvs.v1$heated.padj <=.05] <- "Y"
sig.asvs.v1$bleached.sig <- sig.asvs.v1$bleached.padj #duplicate padj
sig.asvs.v1$bleached.sig[sig.asvs.v1$bleached.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs.v1$bleached.sig[sig.asvs.v1$bleached.padj <=.05] <- "Y"
sig.asvs.v1$bleached.heated.sig <- sig.asvs.v1$bleached.heated.padj #duplicate padj
sig.asvs.v1$bleached.heated.sig[sig.asvs.v1$bleached.heated.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs.v1$bleached.heated.sig[sig.asvs.v1$bleached.heated.padj <=.05] <- "Y"

#merge in taxonomy
sig.asvs.v2 <- cbind(sig.asvs.v1, taxonomy.cull)

#add a "differentially abundant in" column
sig.asvs.v2$DA_in <- paste(sig.asvs.v2$heated.sig, sig.asvs.v2$bleached.sig, sig.asvs.v2$bleached.heated.sig, sep="") #paste the individual significance Y/N characters
#replace three character string with full description
sig.asvs.v2$DA_in[sig.asvs.v2$DA_in=="YNN"] <- "Non-bleached + Heated"
sig.asvs.v2$DA_in[sig.asvs.v2$DA_in=="NYN" | sig.asvs.v2$DA_in=="NAYN"] <- "Bleached + Ambient"
sig.asvs.v2$DA_in[sig.asvs.v2$DA_in=="NNY" | sig.asvs.v2$DA_in=="NANY"] <- "Bleached + Heated"
sig.asvs.v2$DA_in[sig.asvs.v2$DA_in=="YYN"] <- "Non-bleached + Heated and Bleached + Ambient"
sig.asvs.v2$DA_in[sig.asvs.v2$DA_in=="YNY"] <- "Non-bleached + Heated and Bleached + Heated"
sig.asvs.v2$DA_in[sig.asvs.v2$DA_in=="NYY" | sig.asvs.v2$DA_in=="NAYY"] <- "Bleached + Ambient and Bleached + Heated"
sig.asvs.v2$DA_in[sig.asvs.v2$DA_in=="YYY"] <- "Non-bleached + Heated and Bleached + Ambient and Bleached + Heated"

#export.
write.csv(sig.asvs.v2, file.path(dirOutput, "sig.asvs.v2.csv"), )

#calculate mean abundance for each treatment
abund.nosub.coral.tend.t.cull2 <- abund.nosub.coral.tend.t.cull1 #duplicate df
abund.nosub.coral.tend.t.cull2$SampleID <- rownames(abund.nosub.coral.tend.t.cull2) #add sample ID column
abund.nosub.coral.tend.t.cull3 <- merge(abund.nosub.coral.tend.t.cull2, metadata.coral.tend, by.x="SampleID", by.y="Sample_ID") #add in metadata
abund.nosub.coral.tend.mean <- as.data.frame(aggregate(abund.nosub.coral.tend.t.cull3[,c(2:160,174)], by=list(abund.nosub.coral.tend.t.cull3$Treatment), FUN=mean)) #calculate mean

#work up for combining with sig.asvs.v2
abund.nosub.coral.tend.mean1 <- abund.nosub.coral.tend.mean[,-1] #remove first column
rownames(abund.nosub.coral.tend.mean1) <- c("Mean Bleached + Ambient Abundance", "Mean Bleached + Heated Abundance", "Mean Non-bleached + Ambient Abundance", "Mean Non-bleached + Heated Abundance") #adjust rownames
abund.nosub.coral.tend.mean.t <- as.data.frame(t(abund.nosub.coral.tend.mean1)) #transpose
abund.nosub.coral.tend.mean.t$OTU <- rownames(abund.nosub.coral.tend.mean.t) #add in OT column

#combine with sig.asvs.v2 data
abund.nosub.asv <- merge(sig.asvs.v2, abund.nosub.coral.tend.mean.t, by.x="ASV", by.y="OTU") #merge

#add a fontface columns
abund.nosub.asv$sig <- "N" #add new sig column
abund.nosub.asv$sig[abund.nosub.asv$heated.sig=="Y" | abund.nosub.asv$bleached.sig=="Y" | abund.nosub.asv$bleached.heated.sig=="Y"] <- "Y" #subset for all sig asvs in any treatment, replace with "Y"
abund.nosub.asv$fontface <- "plain" #add fontface column
abund.nosub.asv$fontface[abund.nosub.asv$sig=="Y"] <- "bold" #replace plain with bold fontface for sig asvs

#manually make into longformat.
abund.nosub.asv.longformat <- as.data.frame(rep(abund.nosub.asv$ASV, times=3))
abund.nosub.asv.longformat$adjusted_pval <-c(abund.nosub.asv$heated.padj, abund.nosub.asv$bleached.padj, abund.nosub.asv$bleached.heated.padj) #add padjust column
abund.nosub.asv.longformat$Treatment <- c(rep("Non-bleached + Heated", times=159), rep("Bleached + Ambient", times=159), rep("Bleached + Heated", times=159)) #add treatment column
abund.nosub.asv.longformat$l2fc <- c(abund.nosub.asv$heated.l2fc, abund.nosub.asv$bleached.l2fc, abund.nosub.asv$bleached.heated.l2fc) #add l2fc column
abund.nosub.asv.longformat$Mean_Abundance <- c(abund.nosub.asv$`Mean Non-bleached + Heated Abundance`, abund.nosub.asv$`Mean Bleached + Ambient Abundance`, abund.nosub.asv$`Mean Bleached + Heated Abundance`) #add a mean abundance column
abund.nosub.asv.longformat$significant <- c(abund.nosub.asv$heated.sig, abund.nosub.asv$bleached.sig, abund.nosub.asv$bleached.heated.sig) #add significance column
abund.nosub.asv.longformat1 <- cbind(abund.nosub.asv.longformat, rbind(abund.nosub.asv[,14:21], abund.nosub.asv[,14:21], abund.nosub.asv[,14:21])) #add taxonomy columns
abund.nosub.asv.longformat1$Genus_OTU <- paste(abund.nosub.asv.longformat1$Genus, abund.nosub.asv.longformat1$OTUNumber, sep="_")
abund.nosub.asv.longformat1$Family_Genus_OTU <- paste(abund.nosub.asv.longformat1$Family, abund.nosub.asv.longformat1$Genus, abund.nosub.asv.longformat1$OTUNumber, sep="_")
abund.nosub.asv.longformat1$Family_OTU <- paste(abund.nosub.asv.longformat1$Family, abund.nosub.asv.longformat1$OTUNumber, sep="_")
abund.nosub.asv.longformat1$significant1 <- rep(abund.nosub.asv$sig, times=3)
abund.nosub.asv.longformat1$fontface <- rep(abund.nosub.asv$fontface, times=3)

#Next, perform heirarchical clustering on the l2fc values
rownames(abund.nosub.asv) <- abund.nosub.asv$OTUNumber
lfc.clustering <- pheatmap(abund.nosub.asv[,8:10])

cluster.OTUs <- lfc.clustering$tree_row$labels[lfc.clustering$tree_row$order]
cluster.Genus <- abund.nosub.asv$Genus[lfc.clustering$tree_row$order]
cluster.Family <- abund.nosub.asv$Family[lfc.clustering$tree_row$order]
cluster.Genus_OTUs <- paste(cluster.Genus,cluster.OTUs,sep="_")
cluster.Family_Genus_OTUs <- paste(cluster.Family, cluster.Genus, cluster.OTUs,sep="_")
cluster.fontface <- abund.nosub.asv$fontface[lfc.clustering$tree_row$order]

#reorder OTUs in abund.longformat.merged
abund.nosub.asv.longformat1$Genus_OTU <- factor(abund.nosub.asv.longformat1$Genus_OTU, levels=cluster.Genus_OTUs)

#reorder OTUs in abund.nosub.asv
abund.nosub.asv1 <- abund.nosub.asv[lfc.clustering$tree_row$order,]

#Visualize significant ASVs
ggplot(abund.nosub.asv.longformat1, aes(y=Genus_OTU, x=Treatment, size=Mean_Abundance, color=l2fc, group=Genus_OTU))+
  geom_point()+
  scale_size_continuous(limits=c(0,5000))+
  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.325,.5,.675,1), limits=c(-10,10))+
  theme(legend.position="none",axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
  labs(size="Abundance")+
  facet_wrap(.~Class, scales="free_y")

#subset abund.nosub.asv.longformat1 to remove unnecessary Classes and lump "Other taxa
abund.nosub.asv.longformat2 <- subset(abund.nosub.asv.longformat1, Class!="Bacteria_unclassified" & Class!="Acidimicrobiia" & Class!="Proteobacteria_unclassified" & Class!="Oxyphotobacteria")
abund.nosub.asv.longformat2$Class1 <- factor(abund.nosub.asv.longformat2$Class, levels=c(levels(abund.nosub.asv.longformat2$Class), "Other")) #duplicate calss, add new "Other" Class
abund.nosub.asv.longformat2$Class1[abund.nosub.asv.longformat2$Class1=="Deltaproteobacteria" | abund.nosub.asv.longformat2$Class1=="Thermoplasmata"] = "Other" #rename Thermoplasmata an deltaprot as "other"

#visualize again
#Visualize significant ASVs
test=ggplot(abund.nosub.asv.longformat2, aes(y=Family_Genus_OTU, x=Treatment, size=Mean_Abundance, color=l2fc, group=Genus_OTU))+
  geom_point()+
  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.4,.5,.6,1), limits=c(-26,26))+
  theme(axis.text.x=element_text(angle=90,vjust=.5,size = 15), axis.text.y=element_text(face=abund.nosub.asv.longformat2$fontface))+
  labs(size="Abundance", color="Log2 Fold Change")+
  facet_wrap(.~Class1, scales="free_y")
ggsave('ASV bubbpleplot class.png', path=dirFigs, units="in", width=20, height=20, dpi=600)

#visualize again without "Other" AKA just bacterial ASVs no Archeae
ggplot(subset(abund.nosub.asv.longformat2, Class1!="Other"), aes(y=Family_Genus_OTU, x=Treatment, size=Mean_Abundance, color=l2fc, group=Genus_OTU))+
  geom_point()+
  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.4,.5,.6,1), limits=c(-26,26))+
  theme(axis.text.x=element_text(angle=90,vjust=.5,size = 15), axis.text.y=element_text(face=c("bold","plain")))+
  labs(size="Abundance", color="Log2 Fold Change")+
  facet_wrap(.~Class1, scales="free_y")

#visualize again, without other, generating 3 seperate plots by Class and combining with plot_grid
gamma.bubbleplot <- ggplot(subset(abund.nosub.asv.longformat2, Class1=="Gammaproteobacteria"), aes(y=Family_Genus_OTU, x=Treatment, size=Mean_Abundance, fill=l2fc, group=Genus_OTU))+
  geom_point(shape=21)+
  scale_size_continuous(range=c(2,10))+
  scale_fill_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.4,.5,.6,1), limits=c(-26,26))+
  theme_bw()+
  theme(legend.position="none",axis.text.x=element_text(angle=90,vjust=.5,size = 15), axis.text.y=element_text(face=subset(abund.nosub.asv1, Class=="Gammaproteobacteria")$fontface), plot.title=element_text(hjust = 0.5), axis.title.x=element_blank())+
  labs(size="Abundance", color="Log2 Fold Change")+
  ggtitle(label="Gammaproteobacteria")

#visualize again, without other, generating 3 seperate plots by Class and combining with plot_grid
alpha.bubbleplot <- ggplot(subset(abund.nosub.asv.longformat2, Class1=="Alphaproteobacteria"), aes(y=Family_Genus_OTU, x=Treatment, size=Mean_Abundance, fill=l2fc, group=Genus_OTU))+
  geom_point(shape=21)+
  scale_size_continuous(range=c(2,10))+
  scale_fill_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.4,.5,.6,1), limits=c(-26,26))+
  theme_bw()+
  theme(legend.position="none",axis.text.x=element_text(angle=90,vjust=.5,size = 15), axis.text.y=element_text(face=subset(abund.nosub.asv1, Class=="Alphaproteobacteria")$fontface), plot.title=element_text(hjust = 0.5), axis.title.x=element_blank())+
  labs(size="Abundance", color="Log2 Fold Change")+
  ggtitle(label="Alphaproteobacteria")

#visualize again, without other, generating 3 seperate plots by Class and combining with plot_grid
bact.bubbleplot <- ggplot(subset(abund.nosub.asv.longformat2, Class1=="Bacteroidia"), aes(y=Family_Genus_OTU, x=Treatment, size=Mean_Abundance, fill=l2fc, group=Genus_OTU))+
  geom_point(shape=21)+
  scale_size_continuous(range=c(2,10))+
  scale_fill_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.4,.5,.6,1), limits=c(-26,26))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=.5,size = 15), axis.text.y=element_text(face=subset(abund.nosub.asv1, Class=="Bacteroidia")$fontface), plot.title=element_text(hjust = 0.5), axis.title.x=element_blank())+
  labs(size="Abundance", fill="Log2 Fold Change")+
  ggtitle(label="Bacteroidia")

#plot and export
png("../figures/ASV bubbleplot class v2.png", width=21, height=19, units="in", res=600)
plot_grid(alpha.bubbleplot, gamma.bubbleplot, bact.bubbleplot, nrow=1, rel_widths = c(1, 1, 1.2))
dev.off()

#generate volcanoe plots
ggplot(subset(abund.nosub.asv.longformat2, Treatment=="Non-bleached + Heated"), aes(x=l2fc, y=-log10(adjusted_pval)))+
  geom_point()+
  geom_hline(yintercept=1.30103)+
  theme_bw()

ggplot(subset(abund.nosub.asv.longformat2, Treatment=="Bleached + Heated"), aes(x=l2fc, y=-log10(adjusted_pval)))+
  geom_point()+
  geom_hline(yintercept=1.30103)+
  theme_bw()

ggplot(subset(abund.nosub.asv.longformat2, Treatment=="Bleached + Ambient"), aes(x=l2fc, y=-log10(adjusted_pval)))+
  geom_point()+
  geom_hline(yintercept=1.30103)+
  theme_bw()

#visualize significant ASVs in other ways
ggplot(subset(abund.nosub.asv.longformat2, Treatment=="Non-bleached + Heated" & significant1=="Y"), aes(x=Family_OTU, y=l2fc, fill=l2fc))+
  geom_bar(stat="identity")+
  scale_fill_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.4,.5,.6,1), limits=c(-26,26))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90))

ggplot(subset(abund.nosub.asv.longformat2, Treatment=="Bleached + Heated" & significant1=="Y"), aes(x=Family_OTU, y=l2fc, fill=l2fc))+
  geom_bar(stat="identity")+
  scale_fill_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.4,.5,.6,1), limits=c(-26,26))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90))

ggplot(subset(abund.nosub.asv.longformat2, Treatment=="Bleached + Ambient" & significant1=="Y"), aes(x=Family_OTU, y=l2fc, fill=l2fc))+
  geom_bar(stat="identity")+
  scale_fill_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.4,.5,.6,1), limits=c(-26,26))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90))

ggplot(subset(abund.nosub.asv.longformat2, significant1=="Y"), aes(x=Family_OTU, y=l2fc, fill=l2fc))+
  geom_bar(stat="identity", color="black")+
  facet_grid(rows=vars(Treatment))+
  scale_fill_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.4,.5,.6,1), limits=c(-26,26))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90), text=element_text(size=24))+
  ylab(label="log2 fold change")
  #scale_x_discrete(labels=sort(subset(abund.nosub.asv.longformat2, Treatment=="Non-bleached + Heated")$Family[subset(abund.nosub.asv.longformat2, Treatment=="Non-bleached + Heated")$significant1=="Y"], ascending=T))
ggsave("ASV l2fc barplot v1.jpeg", path = dirFigs, width = 20, height = 15, dpi = 600)


ggplot(subset(abund.nosub.asv.longformat2, significant1=="Y"), aes(x=Family_OTU, y=l2fc, fill=l2fc, size=Mean_Abundance))+
  geom_point(color="black", shape=21)+
  scale_size_continuous(range=c(5,20))+
  facet_grid(rows=vars(Treatment))+
  scale_fill_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.4,.5,.6,1), limits=c(-26,26))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90), text=element_text(size=24))+
  ylab(label="log2 fold change")
#scale_x_discrete(labels=sort(subset(abund.nosub.asv.longformat2, Treatment=="Non-bleached + Heated")$Family[subset(abund.nosub.asv.longformat2, Treatment=="Non-bleached + Heated")$significant1=="Y"], ascending=T))
ggsave("ASV l2fc dotplot.jpeg", path = dirFigs, width = 20, height = 15, dpi = 600)

#next, plot sig asvs as boxplots.
relabund.nosub.longformat.sig <- relabund.nosub.longformat[relabund.nosub.longformat$OTU %in% abund.nosub.asv1$ASV[abund.nosub.asv1$sig=="Y"],] #subset for just sig asvs
relabund.nosub.longformat.sig$Family_Genus_OTU <- paste(relabund.nosub.longformat.sig$Family, relabund.nosub.longformat.sig$Genus, relabund.nosub.longformat.sig$OTU, sep="_") #generate new column
relabund.nosub.longformat.sig$Treatment <- factor(relabund.nosub.longformat.sig$Treatment, levels=levels(fact.all.treat)) #adjust treatment factor levels
relabund.nosub.longformat.sig.v1 <- merge(relabund.nosub.longformat.sig, sig.asvs.v2[,c(1,27)], by.x="OTU", by.y="ASV", all.x=T, all.y=F)

ggplot(relabund.nosub.longformat.sig, aes(y=abund, x=Treatment, color=Treatment, fill=Treatment))+
  geom_boxplot()+
  facet_wrap(.~Family_Genus_OTU, scales="free")+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.key.size=unit(2, "cm"), legend.text=element_text(size=20))+
  ylab("Relative Abundance")
#ggsave('Sig ASVs mod.deseq4.jpeg', path = dirFigs, width = 30, height = 22, dpi = 600)

sig.BA.boxplot <- ggplot(subset(relabund.nosub.longformat.sig.v1, DA_in=="Bleached + Ambient"), aes(y=abund, x=Treatment, color=Treatment, fill=Treatment))+
  geom_boxplot()+
  facet_wrap(.~Family_Genus_OTU, scales="free", ncol=2)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.key.size=unit(2, "cm"), legend.text=element_text(size=20), legend.position="none", plot.title=element_text(size=24))
  #ylab("Relative Abundance")
  #ggtitle("Bleached + Ambient")
ggsave('sig.BA.boxplot.jpeg', path=dirFigs, width=9, height=11, dpi=300)

sig.NbH.boxplot <- ggplot(subset(relabund.nosub.longformat.sig.v1, DA_in=="Non-bleached + Heated"), aes(y=abund, x=Treatment, color=Treatment, fill=Treatment))+
  geom_boxplot()+
  facet_wrap(.~Family_Genus_OTU, scales="free", ncol=1)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.key.size=unit(2, "cm"), legend.text=element_text(size=20), legend.position="none", axis.title.y=element_blank(), plot.title=element_text(size=24))
  #ylab("Relative Abundance")+
  #ggtitle("Non-bleached + Heated")
ggsave('sig.NbH.boxplot.jpeg', path=dirFigs, width=4.5, height=1.833333, dpi=300)

sig.BH.boxplot <- ggplot(subset(relabund.nosub.longformat.sig.v1, DA_in=="Bleached + Heated"), aes(y=abund, x=Treatment, color=Treatment, fill=Treatment))+
  geom_boxplot()+
  facet_wrap(.~Family_Genus_OTU, scales="free", ncol=1)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.key.size=unit(2, "cm"), legend.text=element_text(size=20), legend.position="none")
  #ylab("Relative Abundance")+
  #ggtitle("Bleached + Heated")
ggsave('sig.BH.boxplot.jpeg', path=dirFigs, width=4.5, height=9.166667, dpi=300)

sig.NbH.and.BA.boxplot <- ggplot(subset(relabund.nosub.longformat.sig.v1, DA_in=="Non-bleached + Heated and Bleached + Ambient"), aes(y=abund, x=Treatment, color=Treatment, fill=Treatment))+
  geom_boxplot()+
  facet_wrap(.~Family_Genus_OTU, scales="free", ncol=1)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.key.size=unit(2, "cm"), legend.text=element_text(size=20), legend.position="none")
  #ylab("Relative Abundance")+
  #ggtitle("Bleached + Ambient\nNon-bleached + Heated")
  ggsave('sig.NbH.and.BA.boxplot.jpeg', path=dirFigs, width=4.5, height=5.5, dpi=300)

sig.BA.and.BH.boxplot <- ggplot(subset(relabund.nosub.longformat.sig.v1, DA_in=="Bleached + Ambient and Bleached + Heated"), aes(y=abund, x=Treatment, color=Treatment, fill=Treatment))+
  geom_boxplot()+
  facet_wrap(.~Family_Genus_OTU, scales="free", ncol=1)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.key.size=unit(2, "cm"), legend.text=element_text(size=20), legend.position="none")
  #ylab("Relative Abundance")+
  #ggtitle("Bleached + Ambient\nBleached + Heated")
ggsave('sig.BA.and.BH.boxplot.jpeg', path=dirFigs, width=4.5, height=3.666667, dpi=300)

sig.NbH.and.BA.and.BH.boxplot <- ggplot(subset(relabund.nosub.longformat.sig.v1, DA_in=="Non-bleached + Heated and Bleached + Ambient and Bleached + Heated"), aes(y=abund, x=Treatment, color=Treatment, fill=Treatment))+
  geom_boxplot()+
  facet_wrap(.~Family_Genus_OTU, scales="free", ncol=1)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),legend.key.size=unit(2, "cm"), legend.text=element_text(size=20), legend.position="none")
  #ylab("Relative Abundance")+
  #ggtitle("Bleached + Ambient\nBleached + Heated\nNon-bleached + Heated")
ggsave('sig.NbH.and.BA.and.BH.boxplot.jpeg', path=dirFigs, width=4.5, height=7.333333, dpi=300)

sig.NbH.and.BH.boxplot <- ggplot(subset(relabund.nosub.longformat.sig.v1, DA_in=="Non-bleached + Heated and Bleached + Heated"), aes(y=abund, x=Treatment, color=Treatment, fill=Treatment))+
  geom_boxplot()+
  facet_wrap(.~Family_Genus_OTU, scales="free", ncol=1)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.key.size=unit(2, "cm"), legend.text=element_text(size=20), legend.position="none")
  #ylab("Relative Abundance")+
  #ggtitle("Bleached + Heated\nNon-bleached + Heated")
ggsave('sig.NbH.and.BH.boxplot.jpeg', path=dirFigs, width=4.5, height=7.333333, dpi=300)

plot_grid(sig.BA.boxplot, sig.BH.boxplot, sig.NbH.boxplot, sig.BA.and.BH.boxplot, sig.NbH.and.BA.boxplot, sig.NbH.and.BH.boxplot, sig.NbH.and.BA.and.BH.boxplot)

#visualize relabund data as stacked barcharts

#work up relabund data
rownames(relabund) <- relabund$OTUNumber #update rownames
relabund.tend <- relabund[,colnames(relabund) %in% metadata.tend$Sample_Name_Unique] #subset for only abcdom tend 16s samples

#merge in taxonomy data.
relabund.tend$OTU <- rownames(relabund.tend) #add otu column
relabund.tend1 <- merge(relabund.tend, taxonomy1, by.x="OTU", by.y="OTUNumber")

#aggregate by family
relabund.tend1.family <- as.data.frame(aggregate(relabund.tend1[,c(2:16)], by=list(relabund.tend1$Family), FUN=sum)) #aggregate by family
rownames(relabund.tend1.family) <- relabund.tend1.family$Group.1 #update rownames
relabund.tend2.family <- relabund.tend1.family[,-1] #remove family column
relabund.tend2.family.t <- as.data.frame(t(relabund.tend2.family)) #transpose

#cull low abundance familes.
cull.otu(relabund.tend2.family.t, 1, .03, .03) #cull
relabund.tend.family.cull.t <- relabund.df.cull #save output as new df
relabund.tend.family.cull <- as.data.frame(t(relabund.tend.family.cull.t)) #transpose

#work up metadata.tend
metadata.tend$Sample_ID <- metadata.tend$Sample_Name_Unique #add Sample_ID column

#work up taxonomy.
taxonomy.tend.family.cull <- as.data.frame(rownames(relabund.tend.family.cull)) #extract culed fmaily names
colnames(taxonomy.tend.family.cull) <- "Family" #update colnames
taxonomy.tend.family.cull1 <- merge(taxonomy.tend.family.cull, taxonomy1[,4:8], by.x="Family", by.y="Family", all.x=T, all.y=F)
taxonomy.tend.family.cull1$OTUNumber <- taxonomy.tend.family.cull1$Family #create otunumber column for long format generation
taxonomy.tend.family.cull2 <- taxonomy.tend.family.cull1[duplicated(taxonomy.tend.family.cull1)==FALSE,]

#generate longformat of
generate.long.format(relabund.tend.family.cull, metadata.tend, taxonomy.tend.family.cull2)
family.relabund.tend.longformat <- abund.longformat #save as new output

#order samples according to multivariate clustering dendrogram
family.relabund.tend.longformat$Sample <- factor(family.relabund.tend.longformat$Sample, levels=levels(as.factor(family.relabund.tend.longformat$Sample))[clust.tend1$order]) #convert sample to factor and reorder leels to match clustering dendrogram

#visualize
ggplot(family.relabund.tend.longformat, aes(x=Sample, y=abund, fill=Family))+
  geom_bar(stat="identity", color="black")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90, vjust=.5))+
  scale_x_discrete(labels=clust.tend1$labels[clust.tend1$order])+
  scale_fill_manual(values=met.brewer(name="Signac",n=13,type="discrete",direction=-1))+
  ylab("Relative Abundance")
ggsave("16S_Stackedbar_Family_v1.jpg", path=dirFigs, width=10.5, height=10.5, dpi=600)

#work up relabund.tend.family.cull.t to use as a supplemental table in MS
relabund.tend.family.cull.t1 <- relabund.tend.family.cull.t #duplicate df
relabund.tend.family.cull.t1$Sample <- rownames(relabund.tend.family.cull.t)
relabund.tend.family.cull.t2 <- merge(relabund.tend.family.cull.t1, metadata.tend, by.x="Sample", by.y="Sample_Name_Unique") #merge in metadata
relabund.tend.family.cull.t2.mean <- as.data.frame(aggregate(relabund.tend.family.cull.t2[,c(2:14,28)], by=list(relabund.tend.family.cull.t2$Treatment), FUN=mean))
colnames(relabund.tend.family.cull.t2.mean)[1] <- "Treatment" #update colnames
#write.csv(relabund.tend.family.cull.t2.mean, file.path(dirOutput, "family.mean.tend.csv"), ) #export

#generate longformat.
generate.long.format(relabund.tend, metadata.tend, taxonomy1)
relabund.longformat.tend <- abund.longformat #save output as new df

#order samples according to multivariate clustering dendrogram
relabund.longformat.tend$Sample <- factor(relabund.longformat.tend$Sample, levels=levels(as.factor(relabund.longformat.tend$Sample))[clust.tend1$order]) #convert sample to factor and reorder leels to match clustering dendrogram

#visualize
ggplot(relabund.longformat.tend, aes(x=Sample, y=abund, fill=Family))+
  geom_bar(stat="identity", color="black")+
  theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=.5))+
  scale_x_discrete(labels=clust.tend1$labels[clust.tend1$order])+
  xlab("")+
  ylab("Relative Abundance")







#Run DESeq2
#mod.deseq2 <- phyloseq_to_deseq2(physeq.coral.tend.cull1,~ Treatment)
#mod.deseq2 <- estimateSizeFactors(mod.deseq2,type="poscounts")
#mod.deseq2 <- estimateDispersions(mod.deseq2)
#mod.deseq2 <- nbinomWaldTest(mod.deseq2)
#mod.deseq2.heated <- as.data.frame(results(mod.deseq2,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Non-bleached + Ambient", "Non-bleached + Heated")))
#mod.deseq2.bleached <- as.data.frame(results(mod.deseq2,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Non-bleached + Ambient", "Bleached + Ambient")))
#mod.deseq2.bleach.heated=as.data.frame(results(mod.deseq2,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment","Non-bleached + Ambient", "Bleached + Heated")))
#check distributions of pvals
#hist(mod.deseq2.heated$padj, breaks=seq(from=0,to=1,by=.05)) #4 sig ASVs
#hist(mod.deseq2.bleached$padj, breaks=seq(from=0,to=1,by=.05)) #many sig ASVs
#hist(mod.deseq2.bleach.heated$padj, breaks=seq(from=0,to=1,by=.05)) #1 sig ASV

#extract pvalues and add to new df
#sig.asvs <- as.data.frame(rownames(as.data.frame(mod.deseq2.bleached))) #add asv names
#colnames(sig.asvs) <- "ASV" #change colnames
#sig.asvs$heated.p <- mod.deseq2.heated$pvalue
#sig.asvs$heated.padj <- mod.deseq2.heated$padj
#sig.asvs$bleached.p <- mod.deseq2.bleached$pvalue
#sig.asvs$bleached.padj <- mod.deseq2.bleached$padj
#sig.asvs$bleached.heated.p <- mod.deseq2.bleach.heated$pvalue
#sig.asvs$bleached.heated.padj <- mod.deseq2.bleach.heated$padj

#extract l2fc values and add to sig.asvs
#sig.asvs$heated.l2fc <- mod.deseq2.heated$log2FoldChange
#sig.asvs$bleached.l2fc <- mod.deseq2.bleached$log2FoldChange
#sig.asvs$bleached.heated.l2fc <- mod.deseq2.bleach.heated$log2FoldChange

#add new columns indicating significance Y/N
#sig.asvs$heated.sig <- sig.asvs$heated.padj #duplicate padj
#sig.asvs$heated.sig[sig.asvs$heated.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
#sig.asvs$heated.sig[sig.asvs$heated.padj <=.05] <- "Y"
#sig.asvs$bleached.sig <- sig.asvs$bleached.padj #duplicate padj
#sig.asvs$bleached.sig[sig.asvs$bleached.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
#sig.asvs$bleached.sig[sig.asvs$bleached.padj <=.05] <- "Y"
#sig.asvs$bleached.heated.sig <- sig.asvs$bleached.heated.padj #duplicate padj
#sig.asvs$bleached.heated.sig[sig.asvs$bleached.heated.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
#sig.asvs$bleached.heated.sig[sig.asvs$bleached.heated.padj <=.05] <- "Y"

#merge in taxonomy
#sig.asvs1 <- cbind(sig.asvs, taxonomy1.cull)

#add in mean l2fc values for each ASV in each treatment
#Manually calculate log2fold change and visualize.
#workup data for log2foldchange function.
#abund.coral.tend.1.cull1$Treatment <- metadata.coral.tend$Treatment #add treatment column.

#og2.fold.change(abund.coral.tend.1.cull1,101,"Non-bleached + Ambient",100) #run the function

#hist(unlist(log2.fold.df1[,-101]))#check the distribution of log2foldchange data.

#Generate longformat of log2FC data.
#log2.fold.df1.t=as.data.frame(t(log2.fold.df1[,-101])) #work up the data

#abund.lfc.longformat=generate.long.format(log2.fold.df1.t,metadata.coral.tend,taxonomy1.cull) #generate longformat
#colnames(abund.lfc.longformat)[3]="Log2FoldChange" #rename column

#Merge the two abund/lfc longformat dfs.
#work up the dfs.
#abund.raw.longformat$Sample_OTU=paste(abund.raw.longformat$Sample,abund.raw.longformat$OTU)
#abund.lfc.longformat$Sample_OTU=paste(abund.lfc.longformat$Sample,abund.lfc.longformat$OTU)

#abund.longformat.merged=merge(abund.raw.longformat,abund.lfc.longformat,by.x="Sample_OTU",by.y="Sample_OTU") #merge

#abund.longformat.merged$Treatment_OTU <- paste(abund.longformat.merged$Treatment.x, abund.longformat.merged$OTU.x, sep="_") #add a treatment_OTU column

#calculate mean relabund and l2fc for each treatment_OTU
#mean.abund.l2fc <- as.data.frame(aggregate(abund.longformat.merged[,4], by=list(abund.longformat.merged$Treatment_OTU), FUN=mean)) #calculate mean abund
#mean.abund.l2fc1 <- cbind(mean.abund.l2fc, as.data.frame(aggregate(abund.longformat.merged[,36], by=list(abund.longformat.merged$Treatment_OTU), FUN=mean))) #calculate mean l2fc, add to df
#mean.abund.l2fc2 <- mean.abund.l2fc1[,-3] #remove unnecessary columns
#colnames(mean.abund.l2fc2) <- c("Treatment_OTU", "Mean Abundance", "Mean Log2 Fold Change") #rename columns
#mean.abund.l2fc3 <- cbind(mean.abund.l2fc2, t(as.data.frame(strsplit(mean.abund.l2fc2$Treatment_OTU, split="_")))) #split Treatment_OTU
#colnames(mean.abund.l2fc3)[4:5] <- c("Treatment", "OTU") #update colnames

#combine mean.abund.l2fc3 values with sig.asvs1
#first for abundance
#sig.asvs1$mean.abund.heated <- mean.abund.l2fc3$`Mean Abundance`[mean.abund.l2fc3$Treatment == "Non-bleached + Heated"]
#sig.asvs1$mean.abund.bleached <- mean.abund.l2fc3$`Mean Abundance`[mean.abund.l2fc3$Treatment == "Bleached + Ambient"]
#sig.asvs1$mean.abund.bleached.heated <- mean.abund.l2fc3$`Mean Abundance`[mean.abund.l2fc3$Treatment == "Bleached + Heated"]
#sig.asvs1$mean.abund.ambient <- mean.abund.l2fc3$`Mean Abundance`[mean.abund.l2fc3$Treatment == "Non-bleached + Ambient"]

#then for l2fc
#sig.asvs1$mean.l2fc.heated <- mean.abund.l2fc3$`Mean Log2 Fold Change`[mean.abund.l2fc3$Treatment == "Non-bleached + Heated"]
#sig.asvs1$mean.l2fc.bleached <- mean.abund.l2fc3$`Mean Log2 Fold Change`[mean.abund.l2fc3$Treatment == "Bleached + Ambient"]
#sig.asvs1$mean.l2fc.bleached.heated <- mean.abund.l2fc3$`Mean Log2 Fold Change`[mean.abund.l2fc3$Treatment == "Bleached + Heated"]

#compare the l2fc means that I calculated from the ones DESeq2 calcualted.
#First multiply the deseq2 l2fc values by -1 to get the correct sign for comparison.
#sig.asvs1$heated.l2fc.inv <- -1*sig.asvs1$heated.l2fc
#sig.asvs1$bleached.l2fc.inv <- -1*sig.asvs1$bleached.l2fc
#sig.asvs1$bleached.heated.l2fc.inv <- -1*sig.asvs1$bleached.heated.l2fc

#calculate the difference.
#sig.asvs1$l2fc.diff.heated <- sig.asvs1$mean.l2fc.heated - sig.asvs1$heated.l2fc.inv
#sig.asvs1$l2fc.diff.bleached <- sig.asvs1$mean.l2fc.bleached - sig.asvs1$bleached.l2fc.inv
#sig.asvs1$l2fc.diff.bleached.heated <- sig.asvs1$mean.l2fc.bleached.heated - sig.asvs1$bleached.heated.l2fc.inv

#plot a comparison of deseq2 vs. personally calculated l2fc
#bleached.l2fc.comparison.plot <- ggplot(sig.asvs1, aes(x=l2fc.diff.bleached, y=ASV, fill=mean.l2fc.bleached))+
#  geom_bar(stat="identity", color="black")+
#  scale_fill_gradient2(low="blue",mid="white",high="red")+
#  xlab("Caclulated l2fc - DEseq2 l2fc")+
#  ggtitle("Bleached")+
#  theme(axis.text.y=element_blank(), axis.title.y=element_blank())

#heated.l2fc.comparison.plot <- ggplot(sig.asvs1, aes(x=l2fc.diff.heated, y=ASV, fill=mean.l2fc.heated))+
#  geom_bar(stat="identity", color="black")+
#  scale_fill_gradient2(low="blue",mid="white",high="red")+
#  xlab("Caclulated l2fc - DEseq2 l2fc")+
#  ggtitle("Heated")+
#  theme(axis.text.y=element_blank(), axis.title.y=element_blank())

#bleached.heated.l2fc.comparison.plot <- ggplot(sig.asvs1, aes(x=l2fc.diff.bleached.heated, y=ASV, fill=mean.l2fc.bleached.heated))+
#  geom_bar(stat="identity", color="black")+
#  scale_fill_gradient2(low="blue",mid="white",high="red")+
#  xlab("Caclulated l2fc - DEseq2 l2fc")+
#  ggtitle("Bleached+Heated")

#png(filename="../figures/ASV l2fc comparison.png", width=10000, height=10000, res=600)
#plot_grid(bleached.heated.l2fc.comparison.plot, heated.l2fc.comparison.plot, bleached.l2fc.comparison.plot, nrow=1, rel_widths = c(1.3,1,1))
#dev.off()

#Next, try DESEq2 on raw, unsubsampled, un-lulu'ed data.

#First, work up abund.nosub so it only contains the ASVs from the prior, culled, abund df.


#abund.nosub.cull <- abund.nosub2[abund.nosub2$OTU %in% taxonomy1.cull$OTUNumber,] #subset rows for only relevant OTUs


##work up in phyloseq for DESEq2
#abund.nosub.cull1.t <- as.data.frame(t(abund.nosub.cull1)) #transpose
#physeq.abund.nosub.cull1=otu_table(abund.nosub.cull1.t,taxa_are_rows=FALSE) #convert abund object

#physeq.coral.nosub.cull1=phyloseq(physeq.abund.nosub.cull1,physeq.tax.abund.coral.tend.1.cull1,physeq.metadata.coral.tend) #combine

#Run DESEq2 on this data.
#mod.deseq3 <- phyloseq_to_deseq2(physeq.coral.nosub.cull1,~ Treatment)

#mod.deseq3 <- estimateSizeFactors(mod.deseq3,type="poscounts")

#mod.deseq3 <- estimateDispersions(mod.deseq3)

#mod.deseq3 <- nbinomWaldTest(mod.deseq3)

#mod.deseq3.heated <- as.data.frame(results(mod.deseq3,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Non-bleached + Ambient", "Non-bleached + Heated")))
#mod.deseq3.bleached <- as.data.frame(results(mod.deseq3,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Non-bleached + Ambient", "Bleached + Ambient")))
#mod.deseq3.bleach.heated <- as.data.frame(results(mod.deseq3,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment","Non-bleached + Ambient", "Bleached + Heated")))

#visualize pvalues
#hist(mod.deseq3.heated$padj, breaks=seq(from=0,to=1,by=.05)) #6 sig ASVs
#hist(mod.deseq3.bleached$padj, breaks=seq(from=0,to=1,by=.05)) #13? sig ASVs
#hist(mod.deseq3.bleach.heated$padj, breaks=seq(from=0,to=1,by=.05)) #5 sig ASVs

##extract pvalues and add to new df
#sig.asvs2 <- as.data.frame(rownames(as.data.frame(mod.deseq3.bleached))) #add asv names
#colnames(sig.asvs2) <- "ASV" #change colnames
#sig.asvs2$heated.p <- mod.deseq3.heated$pvalue
#sig.asvs2$heated.padj <- mod.deseq3.heated$padj
#sig.asvs2$bleached.p <- mod.deseq3.bleached$pvalue
#sig.asvs2$bleached.padj <- mod.deseq3.bleached$padj
#sig.asvs2$bleached.heated.p <- mod.deseq3.bleach.heated$pvalue
#sig.asvs2$bleached.heated.padj <- mod.deseq3.bleach.heated$padj

#extract l2fc values and add to sig.asvs2
#sig.asvs2$heated.l2fc <- mod.deseq3.heated$log2FoldChange
#sig.asvs2$bleached.l2fc <- mod.deseq3.bleached$log2FoldChange
#sig.asvs2$bleached.heated.l2fc <- mod.deseq3.bleach.heated$log2FoldChange

#add new columns indicating significance Y/N
#sig.asvs2$heated.sig <- sig.asvs2$heated.padj #duplicate padj
#sig.asvs2$heated.sig[sig.asvs2$heated.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
#sig.asvs2$heated.sig[sig.asvs2$heated.padj <=.05] <- "Y"
#sig.asvs2$bleached.sig <- sig.asvs2$bleached.padj #duplicate padj
#sig.asvs2$bleached.sig[sig.asvs2$bleached.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
#sig.asvs2$bleached.sig[sig.asvs2$bleached.padj <=.05] <- "Y"
#sig.asvs2$bleached.heated.sig <- sig.asvs2$bleached.heated.padj #duplicate padj
#sig.asvs2$bleached.heated.sig[sig.asvs2$bleached.heated.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
#sig.asvs2$bleached.heated.sig[sig.asvs2$bleached.heated.padj <=.05] <- "Y"

#export
#write.csv(sig.asvs2, file.path(dirOutput, "sig.asvs2.csv"))

#Visualize
#First, perform heirarchical clustering on the l2fc values
#lfc.clustering=pheatmap(log2.fold.df1[-8:-10,-101])

#lfc.clustering1=pheatmap(log2.fold.df1[-8:-10,-101],cluster_rows = FALSE)

#cluster.OTUs=lfc.clustering$tree_col$labels[lfc.clustering$tree_col$order]
#cluster.samples=lfc.clustering$tree_row$labels[lfc.clustering$tree_row$order]
#cluster.Genus=taxonomy1.cull$Genus[lfc.clustering$tree_col$order]
#cluster.Genus_OTUs=paste(cluster.Genus,cluster.OTUs,sep="_")

#reorder OTUs in abund.longformat.merged
#abund.longformat.merged$OTU.x=factor(abund.longformat.merged$OTU.x,levels=cluster.OTUs)

#add a new column correspondong to family_OTU
#abund.longformat.merged$Genus_OTU=paste(abund.longformat.merged$Genus.x,abund.longformat.merged$OTU.x,sep="_")
#abund.longformat.merged$Genus_OTU=factor(abund.longformat.merged$Genus_OTU,levels=cluster.Genus_OTUs)

#remove control samples for visualization
#abund.longformat.merged1 <- subset(abund.longformat.merged,Sample.x!="ABC_067" & Sample.x!="ABC_068" & Sample.x!="ABC_069")

#pad Genus_OTU strings to 42 characters. Need to do this for alignment of faceted plots to work
#abund.longformat.merged1$Genus_OTU1 <- str_pad(abund.longformat.merged1$Genus_OTU, width=42, "left")

#make all lowercase and then pad to see if this resolves the alignment problem
#abund.longformat.merged1$Genus_OTU2 <- tolower(abund.longformat.merged1$Genus_OTU)
#abund.longformat.merged1$Genus_OTU3 <- str_pad(abund.longformat.merged1$Genus_OTU2, width=42, "left", pad="_")

#Visualize
#ggplot(abund.longformat.merged1,aes(y=Genus_OTU,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
#  geom_point()+
#  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.1375,.55,.6625,1))+
#  theme(axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
#  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
#  xlab(label="Sample")+
#  labs(size="Abundance")
#theme_bw()
#ggsave('ASV_l2fc.jpg', path=dirFigs, width=6, height=16, dpi = 1200)

#Visualize, faceting by Order
#ggplot(abund.longformat.merged1,aes(y=Genus_OTU,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
#  geom_point()+
#  facet_wrap(.~Order.x,scales="free")+
#  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.1375,.55,.6625,1))+
#  theme(axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
#  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
#  xlab(label="Sample")+
#  labs(size="Abundance")
#theme_bw()
#ggsave('ASV_l2fc_order.jpg', path=dirFigs, width=28, height=26, dpi = 600)

#Visualize, faceting by class
#ggplot(abund.longformat.merged1,aes(y=Genus_OTU,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
#  geom_point()+
#  facet_grid(row=vars(Class.x), scales="free_y", space="fre")+
#  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.1375,.55,.6625,1))+
#  theme(axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
#  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
#  xlab(label="Sample")+
#  labs(size="Abundance")
#theme_bw()
#ggsave('ASV_l2fc_class.jpg', path=dirFigs, width=20, height=26, dpi = 600)

#for class, facet out the abundant classes (aprot, gprot, bacteroi) and the not abundantclasses (oxy, prot_unclass, dprot, bact_unclass)
#bubbleplot_class1 <- ggplot(subset(abund.longformat.merged1, Class.x=="Alphaproteobacteria" | Class.x=="Gammaproteobacteria" | Class.x=="Bacteroidia"),aes(y=Genus_OTU1,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
#  geom_point()+
#  facet_wrap(.~Class.x,scales="free_y")+
#  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.325,.5,.675,1), limits=c(-10,10))+
#  theme(axis.text.y=element_text(family="mono"), axis.text.x=element_blank(), axis.title.x=element_blank())+
#  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
#  labs(size="Abundance")

#bubbleplot_class2 <- ggplot(subset(abund.longformat.merged1, Class.x=="Deltaproteobacteria" | Class.x=="Oxyphotobacteria" | Class.x=="Proteobacteria_unclassified" | Class.x=="Bacteria_unclassified"),aes(y=Genus_OTU1,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
#  geom_point()+
#  facet_wrap(.~Class.x,scales="free_y", nrow=2, ncol=3)+
#  scale_size_continuous(limits=c(0,5000))+
#  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.325,.5,.675,1), limits=c(-10,10))+
#  theme(legend.position="none",axis.text.y=element_text(family="mono"),axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
#  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
#  xlab(label="Sample")+
#  labs(size="Abundance")

#combine plots with plot.grid. PROBLEMS WITH ALIGNMENT CUZ OF YAXIS LABELS NEED TO PAD STRINGS
#plot_grid(bubbleplot_class1, bubbleplot_class2, align="v", ncol=1, axis="lr", rel_heights = c(1,.7))
#8x17 landscape

#different method
#bubbleplot_class3 <- ggplot(subset(abund.longformat.merged1, Class.x=="Alphaproteobacteria" | Class.x=="Gammaproteobacteria" | Class.x=="Bacteroidia"),aes(y=Genus_OTU,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
#  geom_point()+
#  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.325,.5,.675,1), limits=c(-10,10))+
#  theme(axis.text.x=element_blank(), axis.title.x=element_blank())+
#  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
#  labs(size="Abundance")+
#  facet_col(vars(Class.x), scales="free_y", space="free")

#bubbleplot_class3_v1 <- ggplot(subset(abund.longformat.merged1, Class.x=="Alphaproteobacteria" | Class.x=="Gammaproteobacteria" | Class.x=="Bacteroidia"),aes(y=Genus_OTU,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
#  geom_point()+
#  scale_size_continuous(limits=c(0,5000))+
#  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.325,.5,.675,1), limits=c(-10,10))+
#  theme(legend.position="none",axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
#  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
#  xlab(label="Sample")+
#  labs(size="Abundance")+
#  facet_col(vars(Class.x), scales="free_y", space="free")

#bubbleplot_class4 <- ggplot(subset(abund.longformat.merged1, Class.x=="Deltaproteobacteria" | Class.x=="Oxyphotobacteria" | Class.x=="Proteobacteria_unclassified" | Class.x=="Bacteria_unclassified"),aes(y=Genus_OTU,x=Sample.x,size=abund,color=Log2FoldChange,group=Genus_OTU))+
#  geom_point()+
#  scale_size_continuous(limits=c(0,5000))+
#  scale_color_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.325,.5,.675,1), limits=c(-10,10))+
#  theme(legend.position="none",axis.text.x=element_text(angle=90,vjust=.5,size = 5))+
#  scale_x_discrete(labels=paste(metadata.coral.tend$Treatment, metadata.coral.tend$Sample_ID, sep=" "))+
#  xlab(label="Sample")+
#  labs(size="Abundance")+
#  facet_col(vars(Class.x), scales="free_y", space="free")

#plotgrid
#plot_grid(bubbleplot_class3, bubbleplot_class4, align="v", ncol=1, axis="lr", rel_heights = c(1,.25))

#Next, visualize boxplots of sig ASVs from mod.deseq2 and mod.deseq3
#sig.asvs3 <- sig.asvs2[is.na(sig.asvs2$heated.p)==FALSE,] #remove NAs

#abund.longformat.merged.sig <- abund.longformat.merged[abund.longformat.merged$OTU.x %in% sig.asvs3$ASV[sig.asvs3$heated.sig == "Y" | sig.asvs3$bleached.sig == "Y" | sig.asvs3$bleached.heated.sig == "Y"],] #subset abund.longformat.merged1 for just mod.deseq3 sig asvs.
#abund.longformat.merged.sig$Treatment.x <- factor(abund.longformat.merged.sig$Treatment.x, levels=levels(fact.all.treat)) #adjust factor levels

#ggplot(abund.longformat.merged.sig, aes(y=abund, x=Treatment.x, color=Treatment.x, fill=Treatment.x))+
#  geom_boxplot()+
#  facet_wrap(.~Genus_OTU, scales="free")+
#  scale_color_manual(values=cost.col.line)+
#  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
#  scale_x_discrete(guide = guide_axis(n.dodge = 2))
#ggsave('Sig ASVs mod.deseq3.jpeg', path = dirFigs, width = 24, height = 18, dpi = 600)






#Cull low abundance ASVs
#Format adbund df for culling.
#rownames(abund.coral.tend) <- abund.coral.tend$Group #update rownames
#abund.coral.tend.1 <- abund.coral.tend[,-1:-3] #remove unnecessary columns

#cull ASVs
#cull.otu(abund.coral.tend.1,3,12,120) #minimum number of reads = 12 (corresponds to relabund .001) in 3 samples or 120 (relabund = .01) in 1 sample. 117 ASVs remain.
#abund.coral.tend.1.cull <- relabund.df.cull #save as new df
#abund.coral.tend.1.cull.t <- as.data.frame(t(abund.coral.tend.1.cull)) #transpose

#because we will be working with log 2 fold change relative to the control, we need to additionally cull ASVs that exhibit high degrees of variance in the control. To do this we will calculate SD for each ASV in the controls and normalize to the mean abund of that ASV in each control (CV).
#sd.abund.coral.tend.1.cull <- as.data.frame(apply(abund.coral.tend.1.cull[8:10,], 2, FUN=sd)) #calculate stderror for each ASV in the controls, save as new df
#colnames(sd.abund.coral.tend.1.cull) <- "SD"

#sd.abund.coral.tend.1.cull$mean <- apply(abund.coral.tend.1.cull[8:10,], 2, FUN=mean) #calculate mean

#sd.abund.coral.tend.1.cull$CV <- sd.abund.coral.tend.1.cull$SD/sd.abund.coral.tend.1.cull$mean #calculate CV

#sd.abund.coral.tend.1.cull$CV[is.na(sd.abund.coral.tend.1.cull$CV)] <- 0 #manually replace NAs with 0s.

#hist(sd.abund.coral.tend.1.cull$CV, breaks=seq(from=0,to=2,by=.1)) ##check distribution. for starters, .9 and 1 seem like a good threshold.

#determine appropriate CV cutoff threshold
#mean and sd of CV
#mean(sd.abund.coral.tend.1.cull$CV) #.54
#sd(sd.abund.coral.tend.1.cull$CV) #.36
#mean(sd.abund.coral.tend.1.cull$CV) + sd(sd.abund.coral.tend.1.cull$CV) #threshold = mean CV + one SD = .9

#cull ASVS who's CV in control is above this threshold
#abund.coral.tend.1.cull1 <- abund.coral.tend.1.cull[,sd.abund.coral.tend.1.cull$CV <= 0.9030586] #that leaves 100 ASVs.

#abund.coral.tend.1.cull1.t <- as.data.frame(t(abund.coral.tend.1.cull1))#transpose
