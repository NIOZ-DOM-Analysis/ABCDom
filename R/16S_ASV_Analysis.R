#install deseq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

#Load libraries
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(cowplot)
library(MetBrewer)
library(igraph)
library(phyloseq)

#-----------------------------------------------------------------------------------------------------------------------------
#load custom functions
source(file="generate.long.format.OTU.R")
source(file="cull.otu.script.R")
source(file="log.2.fold.change.calculation.05272021.R")
source(file="zscore.calculation.script.R")

#-----------------------------------------------------------------------------------------------------------------------------
#read in metadata and data

#read in non-subsampled data
abund.nosub <- read.csv(file.path(dirRAW, "16s", "all_multipletonsFilter_100.count_table.csv"))

relabund <- read.csv(file.path(dirRAW, "16S", "abundance_table_100.database.csv"))

nosub.asv.list <- read.csv(file.path(dirRAW, "16S", "all_multipletonsFilter_100.list.csv")) #read in list converting esvs to OTUs

nosub.taxonomy <- read.csv(file.path(dirRAW, "16S", "all_multipletonsFilter_100.taxonomy.csv")) #read in taxonomy of nosub asvs

#-----------------------------------------------------------------------------------------------------------------------------
#work up the relabund data
rownames(relabund) <- relabund$OTUNumber #add rownames corresponding to OTUs

#subset relabund data for only desired samples
relabund.tend <- relabund[, colnames(relabund) %in% metadata.tend$Sample_Name_Unique] #for tend samples
relabund.coral.tend <- relabund[, colnames(relabund) %in% metadata.coral.tend$Sample_Name_Unique] #for tend samples

#cull relabund data
relabund.tend.t <- as.data.frame(t(relabund.tend)) #transpose
relabund.tend.t.cull <- cull.otu(relabund.tend.t, 3, .001, .01)
relabund.tend.cull <- as.data.frame(t(relabund.tend.t.cull)) #153 ASVs remain

#-----------------------------------------------------------------------------------------------------------------------------
#Work up raw abund data
rownames(abund) <- abund$OTUNumber #add ASV names for rownames

abund.coral.tend <- abund[, colnames(abund) %in% metadata.coral.tend$Sample_Name_Unique] #Subset abund df to correspond to the associated metadata.

#-----------------------------------------------------------------------------------------------------------------------------
#work up and cull raw abundance ASV data
#work up raw abund ASV data
nosub.asv.info <- merge(nosub.asv.list, nosub.taxonomy, by.x="ESV", by.y="ESV") #merge asv list and taxonomy
abund.nosub1 <- merge(abund.nosub, nosub.asv.info, by.x="Representative_Sequence", by.y="ESV") #merge abund.nosub and the nosub.asv.list
abund.nosub2 <- abund.nosub1[,-1:-2] #remove unnecessary columns

#Next, workup abund.nosub.cull so that it only contains the relevant samples.
rownames(abund.nosub2) <- abund.nosub2$OTU #update rownames
abund.nosub.coral.tend <- abund.nosub2[,colnames(abund.nosub2) %in% metadata.coral.tend$Sample_Name_Unique] #subset for relevant samples.

#Next, convert to relabund for later boxplots
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
abund.nosub.coral.tend.t <- as.data.frame(t(abund.nosub.coral.tend)) #transpose.

cull.otu(abund.nosub.coral.tend.t,3,50,1000) #minimum number of reads = 50 in 3 samples or 1000 in 1 sample. 255 ASVs remain.
abund.nosub.coral.tend.t.cull <- relabund.df.cull #save as new df
abund.nosub.coral.tend.cull <- as.data.frame(t(abund.nosub.coral.tend.t.cull)) #transpose

#-----------------------------------------------------------------------------------------------------------------------------
#We will be analyzing log 2 fold change of ASVs relative to the control, so we will need to additionally cull ASVs that exhibit high degrees of variance in the control. To do this we will calculate SD for each ASV in the controls and normalize to the mean abund of that ASV in each control (CV).
sd.abund.nosub.coral.tend.t.cull <- as.data.frame(apply(abund.nosub.coral.tend.t.cull[8:10,], 2, FUN=sd)) #calculate stderror for each ASV in the controls, save as new df
colnames(sd.abund.nosub.coral.tend.t.cull) <- "SD"

sd.abund.nosub.coral.tend.t.cull$mean <- apply(abund.nosub.coral.tend.t.cull[8:10,], 2, FUN=mean) #calculate mean

sd.abund.nosub.coral.tend.t.cull$CV <- sd.abund.nosub.coral.tend.t.cull$SD/sd.abund.nosub.coral.tend.t.cull$mean #calculate CV

sd.abund.nosub.coral.tend.t.cull$CV[is.na(sd.abund.nosub.coral.tend.t.cull$CV)] <- 0 #manually replace NAs with 0s.

#determine appropriate CV cutoff threshold

#mean and sd of CV
mean(sd.abund.nosub.coral.tend.t.cull$CV) #.505
sd(sd.abund.nosub.coral.tend.t.cull$CV) #.405
mean(sd.abund.nosub.coral.tend.t.cull$CV) + sd(sd.abund.nosub.coral.tend.t.cull$CV) #threshold = mean CV + one SD =  0.9111981

#cull ASVS who's CV in control is above this threshold
abund.nosub.coral.tend.t.cull1 <- abund.nosub.coral.tend.t.cull[,sd.abund.nosub.coral.tend.t.cull$CV <=  0.9111981] #cull ASVs, that leaves 159

abund.nosub.coral.tend.cull1 <- as.data.frame(t(abund.nosub.coral.tend.t.cull1)) #transpose

#-----------------------------------------------------------------------------------------------------------------------------
#work up data for deseq2
#Prep taxonomy data for analysis in DESEQ2
taxonomy.cull <- nosub.asv.info[nosub.asv.info$OTU %in% colnames(abund.nosub.coral.tend.t.cull1),] #subset taxonomy for only culled ASVs.
rownames(taxonomy.cull) <- taxonomy.cull$OTU #adjust rownames
colnames(taxonomy.cull)[2] <- "OTUNumber" #adjust colnames
#write.csv(taxonomy.cull, file="taxonomy.cull.csv") #export

#generat longformat of abund data
colnames(metadata.coral.tend)[4] <- "Sample_ID" #prep colnames

rownames(metadata.coral.tend) <- metadata.coral.tend$Sample_Name_Unique #prep rownames

abund.nosub.longformat <- generate.long.format.OTUs(abund.nosub.coral.tend.cull1,metadata.coral.tend,taxonomy.cull) #generate longformat for abund data

relabund.nosub.longformat <- generate.long.format.OTUs(relabund.nosub.coral.tend,metadata.coral.tend,taxonomy.cull) #generate longformat for relabund data

#Convert to phyloseq object for use in DESeq2.
physeq.abund.nosub <- otu_table(abund.nosub.coral.tend.cull1,taxa_are_rows=TRUE) #convert abund df to phyloseq object

physeq.tax.abund.nosub <- tax_table(as.matrix(taxonomy.cull[],)) #convert taxonomy df to phyloseq object #remove fasta, OTUNumber, and mothur taxonomy string from df, convert to matrix prior to converting to physeq object.

physeq.metadata.coral.tend <- sample_data(metadata.coral.tend) #convert to physeq object

physeq.nosub <- phyloseq(physeq.abund.nosub, physeq.tax.abund.nosub, physeq.metadata.coral.tend) #combine

#-----------------------------------------------------------------------------------------------------------------------------
#Run DESEq2 on this data.
mod.deseq4 <- phyloseq_to_deseq2(physeq.nosub,~ Treatment)

mod.deseq4 <- estimateSizeFactors(mod.deseq4,type="poscounts")

mod.deseq4 <- estimateDispersions(mod.deseq4)

mod.deseq4 <- nbinomWaldTest(mod.deseq4)

#extract pairwise model results
mod.deseq4.heated <- as.data.frame(results(mod.deseq4,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Heated", "Control")))
mod.deseq4.bleached <- as.data.frame(results(mod.deseq4,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Bleached", "Control")))
mod.deseq4.bleach.heated <- as.data.frame(results(mod.deseq4,pAdjustMethod="BH",alpha=0.05,contrast=c("Treatment", "Bleached + Heated","Control")))

#extract pvalues and add to new df
sig.asvs.v1 <- as.data.frame(rownames(as.data.frame(mod.deseq4.bleached))) #add asv names

colnames(sig.asvs.v1) <- "ASV" #change colnames

#add p and adjusted p values
sig.asvs.v1$heated.p <- mod.deseq4.heated$pvalue
sig.asvs.v1$heated.padj <- mod.deseq4.heated$padj
sig.asvs.v1$bleached.p <- mod.deseq4.bleached$pvalue
sig.asvs.v1$bleached.padj <- mod.deseq4.bleached$padj
sig.asvs.v1$bleached.heated.p <- mod.deseq4.bleach.heated$pvalue
sig.asvs.v1$bleached.heated.padj <- mod.deseq4.bleach.heated$padj

#extract l2fc values and l2fcSE values and add to sig.asvs.v1
sig.asvs.v1$heated.l2fc <- mod.deseq4.heated$log2FoldChange
sig.asvs.v1$bleached.l2fc <- mod.deseq4.bleached$log2FoldChange
sig.asvs.v1$bleached.heated.l2fc <- mod.deseq4.bleach.heated$log2FoldChange
sig.asvs.v1$heated.l2fcSE <- mod.deseq4.heated$lfcSE
sig.asvs.v1$bleached.l2fcSE <- mod.deseq4.bleached$lfcSE
sig.asvs.v1$bleached.heated.l2fcSE <- mod.deseq4.bleach.heated$lfcSE

#add new columns indicating significance Y/N, based on a padj <= .05 threshold
sig.asvs.v1$heated.sig <- sig.asvs.v1$heated.padj #duplicate padj
sig.asvs.v1$heated.sig[sig.asvs.v1$heated.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs.v1$heated.sig[sig.asvs.v1$heated.padj <=.05] <- "Y"
sig.asvs.v1$bleached.sig <- sig.asvs.v1$bleached.padj #duplicate padj
sig.asvs.v1$bleached.sig[sig.asvs.v1$bleached.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs.v1$bleached.sig[sig.asvs.v1$bleached.padj <=.05] <- "Y"
sig.asvs.v1$bleached.heated.sig <- sig.asvs.v1$bleached.heated.padj #duplicate padj
sig.asvs.v1$bleached.heated.sig[sig.asvs.v1$bleached.heated.padj >= .05] <- "N" #replace padj values that are greater than .05 w N
sig.asvs.v1$bleached.heated.sig[sig.asvs.v1$bleached.heated.padj <=.05] <- "Y"

#add in taxonomy
sig.asvs.v2 <- cbind(sig.asvs.v1, taxonomy.cull[,1:8])

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

#-----------------------------------------------------------------------------------------------------------------------------
#work up abund data for visualization with deseq2 results
#calculate mean abundance for each treatment
abund.nosub.coral.tend.t.cull2 <- abund.nosub.coral.tend.t.cull1 #duplicate df
abund.nosub.coral.tend.t.cull2$SampleID <- rownames(abund.nosub.coral.tend.t.cull2) #add sample ID column
abund.nosub.coral.tend.t.cull3 <- merge(abund.nosub.coral.tend.t.cull2, metadata.coral.tend, by.x="SampleID", by.y="Sample_Name_Unique") #add in metadata
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
abund.nosub.asv.longformat1 <- cbind(abund.nosub.asv.longformat, rbind(abund.nosub.asv[,18:24], abund.nosub.asv[,18:24], abund.nosub.asv[,18:24])) #add taxonomy columns
abund.nosub.asv.longformat1$Genus_OTU <- paste(abund.nosub.asv.longformat1$Genus, abund.nosub.asv.longformat1$OTUNumber, sep="_")
abund.nosub.asv.longformat1$Family_Genus_OTU <- paste(abund.nosub.asv.longformat1$Family, abund.nosub.asv.longformat1$Genus, abund.nosub.asv.longformat1$OTUNumber, sep="_")
abund.nosub.asv.longformat1$Family_OTU <- paste(abund.nosub.asv.longformat1$Family, abund.nosub.asv.longformat1$OTUNumber, sep="_")
abund.nosub.asv.longformat1$significant1 <- rep(abund.nosub.asv$sig, times=3)
abund.nosub.asv.longformat1$fontface <- rep(abund.nosub.asv$fontface, times=3)
abund.nosub.asv.longformat1$l2fcSE <- c(abund.nosub.asv$heated.l2fcSE, abund.nosub.asv$bleached.l2fcSE, abund.nosub.asv$bleached.heated.l2fcSE) #add l2fcSE columns

#Next, perform heirarchical clustering on the l2fc values
rownames(abund.nosub.asv) <- abund.nosub.asv$OTUNumber

lfc.clustering <- pheatmap(abund.nosub.asv[,8:10]) #cluster

#work up clustering data
cluster.OTUs <- lfc.clustering$tree_row$labels[lfc.clustering$tree_row$order]
cluster.Genus <- abund.nosub.asv$Genus[lfc.clustering$tree_row$order]
cluster.Family <- abund.nosub.asv$Family[lfc.clustering$tree_row$order]
cluster.Genus_OTUs <- paste(cluster.Genus,cluster.OTUs,sep="_")
cluster.Family_Genus_OTUs <- paste(cluster.Family, cluster.Genus, cluster.OTUs,sep="_")
cluster.fontface <- abund.nosub.asv$fontface[lfc.clustering$tree_row$order]

#reorder OTUs in abund.longformat.merged
abund.nosub.asv.longformat1$Family_Genus_OTU <- factor(abund.nosub.asv.longformat1$Family_Genus_OTU, levels=cluster.Family_Genus_OTUs)

#reorder OTUs in abund.nosub.asv
abund.nosub.asv1 <- abund.nosub.asv[lfc.clustering$tree_row$order,]

#-----------------------------------------------------------------------------------------------------------------------------
#Visualize significant ASVs
#subset abund.nosub.asv.longformat1 to remove unnecessary Classes and lump "Other taxa
abund.nosub.asv.longformat2 <- subset(abund.nosub.asv.longformat1, Class!="Bacteria_unclassified" & Class!="Acidimicrobiia" & Class!="Proteobacteria_unclassified" & Class!="Oxyphotobacteria")
#abund.nosub.asv.longformat2$Class1 <- factor(abund.nosub.asv.longformat2$Class, levels=c(levels(abund.nosub.asv.longformat2$Class), "Other")) #duplicate calss, add new "Other" Class
#abund.nosub.asv.longformat2$Class1[abund.nosub.asv.longformat2$Class1=="Deltaproteobacteria" | abund.nosub.asv.longformat2$Class1=="Thermoplasmata"] = "Other" #rename Thermoplasmata an deltaprot as "other"
abund.nosub.asv.longformat2$Treatment[abund.nosub.asv.longformat2$Treatment=="Bleached + Ambient"] <- "Bleached"
abund.nosub.asv.longformat2$Treatment[abund.nosub.asv.longformat2$Treatment=="Non-bleached + Heated"] <- "Heated"

#visualize again, without other, generating 3 seperate plots by Class and combining with plot_grid
gamma.bubbleplot <- ggplot(subset(abund.nosub.asv.longformat2, Class=="Gammaproteobacteria"), aes(y=Family_Genus_OTU, x=Treatment, size=Mean_Abundance, fill=l2fc, group=Genus_OTU))+
  geom_point(shape=21)+
  scale_size_continuous(range=c(2,10))+
  scale_fill_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.4,.5,.6,1), limits=c(-26,26))+
  theme_bw()+
  theme(legend.position="none",axis.text.x=element_text(angle=90,vjust=.5,size = 15), axis.text.y=element_text(face=subset(abund.nosub.asv1, Class=="Gammaproteobacteria")$fontface), plot.title=element_text(hjust = 0.5), axis.title.x=element_blank())+
  labs(size="Abundance", color="Log2 Fold Change")+
  ggtitle(label="Gammaproteobacteria")

#visualize again, without other, generating 3 seperate plots by Class and combining with plot_grid
alpha.bubbleplot <- ggplot(subset(abund.nosub.asv.longformat2, Class=="Alphaproteobacteria"), aes(y=Family_Genus_OTU, x=Treatment, size=Mean_Abundance, fill=l2fc, group=Genus_OTU))+
  geom_point(shape=21)+
  scale_size_continuous(range=c(2,10))+
  scale_fill_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.4,.5,.6,1), limits=c(-26,26))+
  theme_bw()+
  theme(legend.position="none",axis.text.x=element_text(angle=90,vjust=.5,size = 15), axis.text.y=element_text(face=subset(abund.nosub.asv1, Class=="Alphaproteobacteria")$fontface), plot.title=element_text(hjust = 0.5), axis.title.x=element_blank())+
  labs(size="Abundance", color="Log2 Fold Change")+
  ggtitle(label="Alphaproteobacteria")

#visualize again, without other, generating 3 seperate plots by Class and combining with plot_grid
bact.bubbleplot <- ggplot(subset(abund.nosub.asv.longformat2, Class=="Bacteroidia"), aes(y=Family_Genus_OTU, x=Treatment, size=Mean_Abundance, fill=l2fc, group=Genus_OTU))+
  geom_point(shape=21)+
  scale_size_continuous(range=c(2,10))+
  scale_fill_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.4,.5,.6,1), limits=c(-26,26))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=.5,size = 15), axis.text.y=element_text(face=subset(abund.nosub.asv1, Class=="Bacteroidia")$fontface), plot.title=element_text(hjust = 0.5), axis.title.x=element_blank())+
  labs(size="Abundance", fill="Log2 Fold Change")+
  ggtitle(label="Bacteroidia")

#plot and export, use for supplemental figure 6
#png("ASV bubbleplot class v2.png", width=21, height=19, units="in", res=600)
#plot_grid(alpha.bubbleplot, gamma.bubbleplot, bact.bubbleplot, nrow=1, rel_widths = c(1, 1, 1.2))
#dev.off()

#subset for just significant asvs
abund.nosub.asv.longformat2.sig <- subset(abund.nosub.asv.longformat2, significant1=="Y")
abund.nosub.asv.longformat2.sig$Genus_OTU <- factor(abund.nosub.asv.longformat2.sig$Genus_OTU, levels=abund.nosub.asv.longformat2.sig$Genus_OTU[1:31][c(20,19,11,21,17,6,13,9,14,25,5,10,16,24,28,29,1,23,12,22,2,4,18,3,7,8,26,27,30,31,15)])

#generate dotplto for figure 3C
ggplot(abund.nosub.asv.longformat2.sig, aes(x=Genus_OTU, y=l2fc, fill=l2fc, size=Mean_Abundance))+
  geom_point(color="black", shape=21)+
  scale_size_continuous(range=c(5,20))+
  geom_linerange(color="black", size=.5, aes(ymin=l2fc-l2fcSE, ymax=l2fc+l2fcSE))+
  scale_size_continuous(range=c(5,20))+
  facet_grid(rows=vars(Treatment))+
  scale_fill_gradientn(colours=c("blue4","blue","white","red","red4"),values=c(0,.4,.5,.6,1), limits=c(-26,26))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1), text=element_text(size=24))+
  ylab(label="log2 fold change")+
  xlab("")+
  geom_hline(yintercept=0, linetype="dashed")
#ggsave("Figure 3C.jpeg", path = dirFigs, width = 20, height = 15, dpi = 600)

#export data used to make figure 3C
#write.csv(abund.nosub.asv.longformat2.sig, file="figure.3c.data")

#visualize relabund data as 2 way heatmaps at the family level
#work up taxonomy data from relabund data frame
taxonomy <- as.data.frame(t(as.data.frame(strsplit(relabund$OTUConTaxonomy, split=";")))) #extract taxonomy string and separate taxonomy levels
taxonomy$OTU <- relabund$OTUNumber #add OTU names
rownames(taxonomy) <- taxonomy$OTU #update rownames
colnames(taxonomy) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "OTUNumber") #update colnames

#merge in taxonomy data with relabund data
relabund.tend$OTU <- rownames(relabund.tend) #add otu column
relabund.tend1 <- merge(relabund.tend, taxonomy, by.x="OTU", by.y="OTUNumber")

#aggregate by family
relabund.tend1.family <- as.data.frame(aggregate(relabund.tend1[,c(2:16)], by=list(relabund.tend1$Family), FUN=sum)) #aggregate by family
rownames(relabund.tend1.family) <- relabund.tend1.family$Group.1 #update rownames
relabund.tend2.family <- relabund.tend1.family[,-1] #remove family column
relabund.tend2.family.t <- as.data.frame(t(relabund.tend2.family)) #transpose

#cull low abundance familes.
cull.otu(relabund.tend2.family.t, 2, .005, .1) #cull
relabund.tend.family.cull.t <- relabund.df.cull #save output as new df
relabund.tend.family.cull <- as.data.frame(t(relabund.tend.family.cull.t)) #transpose

#add a column corresponding to culled families called "other".
calculate.rare <- function(x) { #generate a function to calculate rare, "other" portions of community

  #rare.vec <- c()

  rare.vec <- 1-sum(x)

  rare.vec <<- rare.vec

}

relabund.tend.family.cull.rare <- apply(relabund.tend.family.cull, 2, FUN=calculate.rare) #apply function to df
relabund.tend.family.cull1 <- rbind(relabund.tend.family.cull, relabund.tend.family.cull.rare) #merge rare row back with df.
rownames(relabund.tend.family.cull1)[24] <- "Other" #rename rare families as "other"

#work up metadata.tend for longformat
metadata.tend$Sample_ID <- metadata.tend$Sample_Name_Unique #add Sample_ID column

#work up taxonomy for longformat
taxonomy.tend.family.cull <- as.data.frame(rownames(relabund.tend.family.cull1)) #extract culled family names
#colnames(taxonomy.tend.family.cull) <- "Family" #update colnames
#taxonomy.tend.family.cull1 <- merge(taxonomy.tend.family.cull, taxonomy[,1:5], by.x="Family", by.y="Family", all.x=T, all.y=F)
taxonomy.tend.family.cull$OTUNumber <- taxonomy.tend.family.cull$Family #create otunumber column for long format generation
#taxonomy.tend.family.cull2 <- taxonomy.tend.family.cull1[duplicated(taxonomy.tend.family.cull1)==FALSE,]

#generate summed family clustering heatmap, mean for each treatment
#work up relabund data
relabund.tend.family.cull1.t <- as.data.frame(t(relabund.tend.family.cull1)) #transpose
relabund.tend.family.cull1.t$Treatment <- metadata.tend$Treatment #add treatment column

relabund.tend.family.cull1.t.mean <- aggregate(relabund.tend.family.cull1.t[,-25], by=list(relabund.tend.family.cull1.t$Treatment), FUN=mean) #calculate mean relabund for each family in each treatment
colnames(relabund.tend.family.cull1.t.mean)[1] <- "Treatment" #update colnames
rownames(relabund.tend.family.cull1.t.mean) <- relabund.tend.family.cull1.t.mean$Treatment #update rownames
rownames(relabund.tend.family.cull1.t.mean)[c(1:2,4:6)] <- c("Negative Control", "Bleached", "Negative Control + Heated", "Control", "Heated") #update rownames

relabund.tend.family.cull1.t.mean1 <- relabund.tend.family.cull1.t.mean[,-1] #remove treatment column

#prepare the data for 2 way heatmap
relabund.tend.family.cull1.t.mean.trans <- asin(sqrt(relabund.tend.family.cull1.t.mean1)) #arcsine sqrt transform

relabund.tend.family.cull1.t.mean.zscore <- apply(relabund.tend.family.cull1.t.mean.trans, 2, FUN=zscore.calculation) #zscore

#visualize and save as figure 3B
family.cull.heatmap <- pheatmap(relabund.tend.family.cull1.t.mean.zscore)
#ggsave("Figure 3B.png", family.cull.heatmap, dpi=600, width=6, height=4) #UPDATE DIRECTORY IN FINAL REPO

#generate longformat of the family aggregated relabund data
colnames(taxonomy.tend.family.cull) <- "OTUNumber" #work up taxonomy data

generate.long.format.OTUs(relabund.tend.family.cull1, metadata.tend, taxonomy.tend.family.cull)
family.relabund.tend.longformat <- abund.longformat #save as new output
colnames(family.relabund.tend.longformat)[2] <- "Family" #change colnames

#generate pie charts for figure 5
ggplot(subset(family.relabund.tend.longformat, Treatment=="Non-bleached + Ambient"), aes(x="", y=abund, fill=Family))+
  geom_bar(stat="summary", fun.y="mean", color="black")+
  theme_void()+
  theme(legend.position="none")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("gold", met.brewer(name="Signac",n=12,type="discrete",direction=-1), met.brewer(name="Renoir",n=4,type="discrete",direction=-1), "gray", "orange", "firebrick1", "darkolivegreen3", met.brewer(name="Juarez",n=6,type="discrete",direction=1)))
#FIX DIRECTORIES FOR FINAL SCRIPTggsave("Figure 5 piechart control.png", path=dirFigs, height=6, width=6, dpi=600)

ggplot(subset(family.relabund.tend.longformat, Treatment=="Non-bleached + Heated"), aes(x="", y=abund, fill=Family))+
  geom_bar(stat="summary", fun.y="mean", color="black")+
  theme_void()+
  theme(legend.position="none")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("gold", met.brewer(name="Signac",n=12,type="discrete",direction=-1), met.brewer(name="Renoir",n=4,type="discrete",direction=-1), "gray", "orange", "firebrick1", "darkolivegreen3", met.brewer(name="Juarez",n=6,type="discrete",direction=1)))
#FIX DIRECTORIES FOR FINAL SCRIPTggsave("Figure 5 piechart heated.png", path=dirFigs, height=6, width=6, dpi=600)

ggplot(subset(family.relabund.tend.longformat, Treatment=="Bleached + Ambient"), aes(x="", y=abund, fill=Family))+
  geom_bar(stat="summary", fun.y="mean", color="black")+
  theme_void()+
  theme(legend.position="none")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("gold", met.brewer(name="Signac",n=12,type="discrete",direction=-1), met.brewer(name="Renoir",n=4,type="discrete",direction=-1), "gray", "orange", "firebrick1", "darkolivegreen3", met.brewer(name="Juarez",n=6,type="discrete",direction=1)))
#FIX DIRECTORIES FOR FINAL SCRIPTggsave("Figure 5 piechart bleached.png", path=dirFigs, height=6, width=6, dpi=600)

ggplot(subset(family.relabund.tend.longformat, Treatment=="Bleached + Heated"), aes(x="", y=abund, fill=Family))+
  geom_bar(stat="summary", fun.y="mean", color="black")+
  theme_void()+
  theme(legend.position="none")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("gold", met.brewer(name="Signac",n=12,type="discrete",direction=-1), met.brewer(name="Renoir",n=4,type="discrete",direction=-1), "gray", "orange", "firebrick1", "darkolivegreen3", met.brewer(name="Juarez",n=6,type="discrete",direction=1)))
#FIX DIRECTORIES FOR FINAL SCRIPTggsave("Figure 5 piechart bleached and heated.png", path=dirFigs, height=6, width=6, dpi=600)

#export data
#FIX OUTPUT DIRECTORY write.csv(family.relabund.tend.longformat, file="figure.5.relativeabundance.data.csv")

#work up supplemental figure 7
#add a column to sig asvs that includes if it's enriched or depleted.
sig.asvs.v2$heated.enrich <- rep("", times=159) #create dummy column
sig.asvs.v2$heated.enrich[sig.asvs.v2$heated.l2fc>0] <- "enriched" #for all heated l2fc values greater than 0, code as "enriched"
sig.asvs.v2$heated.enrich[sig.asvs.v2$heated.l2fc<0] <- "depleted" #do the opposite for depleted
#repeat for bleached
sig.asvs.v2$bleached.enrich <- rep("", times=159) #create dummy column
sig.asvs.v2$bleached.enrich[sig.asvs.v2$bleached.l2fc>0] <- "enriched"
sig.asvs.v2$bleached.enrich[sig.asvs.v2$bleached.l2fc<0] <- "depleted"
#repeat for bleached.heated
sig.asvs.v2$bleached.heated.enrich <- rep("", times=159) #create dummy column
sig.asvs.v2$bleached.heated.enrich[sig.asvs.v2$bleached.heated.l2fc>0] <- "enriched"
sig.asvs.v2$bleached.heated.enrich[sig.asvs.v2$bleached.heated.l2fc<0] <- "depleted"

#merge sig.asvs.v2 with abund.nosub.asv1
sig.asvs.v3 <- merge(sig.asvs.v2, abund.nosub.asv1[,c(1,26:30)], by.x="ASV", by.y="ASV", all.x=T, all.y=T)

#
sig.asvs.v3$sig.enrichment <- rep("", times=159) #create dummy variable
for (i in 1:nrow(sig.asvs.v3)) {

  sig.asvs.v3$sig.enrichment[i] <- paste(ifelse(sig.asvs.v3$heated.sig[i]=="Y", sig.asvs.v3$heated.enrich[i], ""), ifelse(sig.asvs.v3$bleached.sig[i]=="Y", sig.asvs.v3$bleached.enrich[i], ""), ifelse(sig.asvs.v3$bleached.heated.sig[i]=="Y", sig.asvs.v3$bleached.heated.enrich[i], ""), sep="")

}

#replace duplicates with just one
sig.asvs.v3$sig.enrichment[sig.asvs.v3$sig.enrichment=="NAenrichedenriched" | sig.asvs.v3$sig.enrichment=="NAenriched" | sig.asvs.v3$sig.enrichment=="enrichedenrichedenriched" | sig.asvs.v3$sig.enrichment=="enrichedenriched" | sig.asvs.v3$sig.enrichment=="enriched"] <- "enriched"
sig.asvs.v3$sig.enrichment[sig.asvs.v3$sig.enrichment=="NAdepleted" | sig.asvs.v3$sig.enrichment=="depleteddepleted" | sig.asvs.v3$sig.enrichment=="depleted"] <- "depleted"

#merge sig ASVs with relabund.nosub.longformat
relabund.nosub.longformat.sig <- merge(relabund.nosub.longformat, sig.asvs.v3, by.x="OTU", by.y="ASV", all.x=T, all.y=F)
colnames(relabund.nosub.longformat.sig)[9] <- "Family" #update colnames for Family

#add DA_in names vector
DA_names <- c(
  `Non-bleached + Heated` = "Heated",
  `Bleached + Heated` = "Bleached + Heated",
  `Bleached + Ambient` = "Bleached",
  `Bleached + Ambient and Bleached + Heated` = "Bleached, Bleached + Heated",
  `Non-bleached + Heated and Bleached + Ambient` = "Heated, Bleached",
  `Non-bleached + Heated and Bleached + Ambient and Bleached + Heated` = "Heated, Bleached, Bleached + Heated",
  `Non-bleached + Heated and Bleached + Heated` = "Heated, Bleached + Heated"
)

ggplot(subset(relabund.nosub.longformat.sig, sig=="Y"), aes(y=abund, x=Treatment, fill=Family, group=OTU))+
  facet_grid(cols=vars(sig.enrichment), rows=vars(DA_in), labeller=as_labeller(DA_names), scale="free")+
  geom_bar(stat="summary", color="black", fun.y="mean")+
  scale_fill_manual(values=met.brewer(name="Signac",n=13,type="discrete",direction=1))+
  ylab("Relative Abundance")+
  theme(strip.text.y = element_text(size = 8))+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))
#ggsave("Figure S7.png", path=dirFigs, width=9, height=15, dpi=600)

#-----------------------------------------------------------------------------------------------------------------------------
#work up relabund.tend.family.cull.t to use as a supplemental table S1
relabund.tend.family.cull.t1 <- relabund.tend.family.cull.t #duplicate df

relabund.tend.family.cull.t1$Sample <- rownames(relabund.tend.family.cull.t) #add in sample column

relabund.tend.family.cull.t2 <- merge(relabund.tend.family.cull.t1, metadata.tend, by.x="Sample", by.y="Sample_Name_Unique") #merge in metadata

relabund.tend.family.cull.t2.mean <- as.data.frame(aggregate(relabund.tend.family.cull.t2[,2:24], by=list(relabund.tend.family.cull.t2$Treatment), FUN=mean)) #calculate mean

colnames(relabund.tend.family.cull.t2.mean)[1] <- "Treatment" #update colnames

#write.csv(relabund.tend.family.cull.t2.mean, "table.s1.dat.csv") #export FIX DIRECTORY

#-----------------------------------------------------------------------------------------------------------------------------
#generate network for supplement
#Work up abund.nosub.coral.tend.cull1 to include water samples for network visualization.
abund.nosub.tend <- abund.nosub2[,colnames(abund.nosub2) %in% metadata.tend$Sample_Name_Unique] #subset for just tend samples




#DO THIS IN SEPERATE R SCRIPT, EXPORT ABUND NOSUB 2 AND WORK UP.






abund.nosub.tend.cull <- abund.nosub.tend[rownames(abund.nosub.tend) %in% rownames(abund.nosub.coral.tend.cull1),] #subset for just 159 ASVs run thru DESEQ2

#Construct network using spiec.easi on a phyloseq object.
taxonomy.cull.physeq <- tax_table(as.matrix(taxonomy.cull))
abund.nosub.tend.cull.physeq <- otu_table(abund.nosub.tend.cull, taxa_are_rows=T)
rownames(metadata.tend) <- metadata.tend$Sample_Name_Unique
metadata.tend.physeq <- sample_data(metadata.tend)
physeq.cull.tend <- phyloseq(taxonomy.cull.physeq, abund.nosub.tend.cull.physeq, metadata.tend.physeq)

#Construct network using spiec.easi.
#net.c <- spiec.easi(as.matrix(t(abund.nosub.tend.cull)), method='mb') #construct the network
net.c <- spiec.easi(physeq.cull.tend, method='mb')
net.c1 <- spiec.easi(as.matrix(t(abund.nosub.tend.cull)), method='mb')

net.c.ig <- adj2igraph(getRefit(net.c), vertex.attr=list(name=taxa_names(physeq.cull.tend))) #fit the network for igraph
net.c.ig <- adj2igraph(net.c) #fit the network for igraph







net.c <- spiec.easi(as.matrix(t(abund.nosub.tend.cull)), method='mb')
>>>>>>> 1fb7be5ce0ab88fa4b40924424051f048e0d757a

net.c$lambda
net.c$select$stars$summary

<<<<<<< HEAD
net.c.ig <- adj2igraph(getRefit(net.c), vertex.attr=list(name=taxa_names(physeq.cull.tend)))

plot(net.c.ig, vertex.size=8, type="taxa", main="MB")

#export
write.graph(net.c.ig, file=file.path(dirOutput,"net.c.ig.txt"),format="ncol")
















#clear workspace
rm(list=ls())
