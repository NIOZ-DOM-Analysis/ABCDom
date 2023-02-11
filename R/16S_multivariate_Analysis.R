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
library(pheatmap)
library(viridis)
library(ggforce)

#source custom functions
#source(file="generate.square.dist.script.07.27.2020.R")

#read in metadata and data
#Read in unifrac.
#unifrac <- read.csv(file.path(dirRAW, "16S", "otu_repr_100.tre1.wsummary.csv"))

#read in ASV table
#abund <- read.csv(file.path(dirRAW, "16S", "abundance_table_100.database.csv"))

#read in alpha diversity data
#adiv <- read.csv(file.path(dirRAW, "16S", "all_alphaDiversity_100.summary.csv"))

#Read in metadata. Unnecessary since metadata is included in global environment
#metadata <- read.csv(file="metadata1.csv")
#Read in the 16S metadata
#metadata_16S <- read.csv(file.path(dirOutput, "metadata_16S_abcDOM.csv"))

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

#export 16s metadata
#write.csv(metadata1, file.path(dirOutput, "metadata1_16S.csv"))

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
legend(0, .5, legend=levels(metadata.tend$Treatment), col=cost.col.line, pt.bg=cost.col.fill, pch=21, pt.cex=2, pt.lwd=2.75, cex=.75, y.intersp=1.5, bty='n')
dev.off()

#plot heirarchical clustering plot, save as jpeg
jpeg("../figures/16S_dendogram_v1.jpg",width=2100, height=2100, res=300)
clust.tend2 <- as.dendrogram(clust.tend1, hang=-1) #convert to dendrogram
plot(clust.tend2, type="rectangle", dLeaf=.1, ylim=c(-.5,1))
symbols(1:15, rep(-.05,15), circles=rep(1,15), add=T, fg=cost.col.line[metadata.tend$Treatment][clust.tend1$order], bg=cost.col.fill[metadata.tend$Treatment][clust.tend1$order], inches=.09, xpd=T, lwd=3)
dev.off()

#plot nmds of t0 and tend samples for supplement
#First workup nmds
nmds.16S <- metaMDS(unifrac.dist1, k=2, trymax=100) #generate nmds. stress=.04
nmds.16S.scores <- as.data.frame(scores(nmds.16S)) #extract scores
nmds.16S.scores$Sample <- rownames(nmds.16S.scores) #add sample names column
nmds.16S.scores1 <- merge(nmds.16S.scores, metadata1, by.x="Sample", by.y="Sample_Name_Unique", all.x=T, all.y=F) #merge with metadata

#plot nmds.
ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data=nmds.16S.scores1,
             aes(x=NMDS1, y=NMDS2, colour=Timepoint_char),
             size=4, stroke=1.5, shape=16) +
  scale_color_manual(values = c("black", "orange"), name = "Timepoint")+
  ggtitle("abcDOM 16S T0-Tfinal")+
  coord_fixed(ratio=1)
ggsave("NMDS_16S_all_data.jpg", path = dirFigs, width = 6.75, height = 5, units = "in", dpi = 320)

#plot nmds of tend samples for supplement
#First workup nmds
nmds.16S.tend <- metaMDS(unifrac.dist.tend, k=2, trymax=100) #generate nmds. stress=.09
nmds.16S.tend.scores <- as.data.frame(scores(nmds.16S.tend)) #extract scores
nmds.16S.tend.scores$Sample <- rownames(nmds.16S.tend.scores) #add sample names column
nmds.16S.tend.scores1 <- merge(nmds.16S.tend.scores, metadata1, by.x="Sample", by.y="Sample_Name_Unique", all.x=T, all.y=F) #merge with metadata
nmds.16S.tend.scores1$Stress_status_v1 <- nmds.16S.tend.scores1$Stress_status #duplicate stress status column
nmds.16S.tend.scores1$Stress_status_v1 <- factor(nmds.16S.tend.scores1$Stress_status_v1, levels=c(levels(nmds.16S.tend.scores1$Stress_status_v1), "Stressed")) #update levels
nmds.16S.tend.scores1$Stress_status_v1[nmds.16S.tend.scores1$Origin_PlanC=="control"] <- NA
nmds.16S.tend.scores1$Stress_status_v1[nmds.16S.tend.scores1$Treatment=="Bleached + Heated" | nmds.16S.tend.scores1$Treatment=="Bleached + Ambient" | nmds.16S.tend.scores1$Treatment=="Non-bleached + Heated"] <- "Stressed" #assign a new value for all stress treatments

#plot nmds.
ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data=nmds.16S.tend.scores1,
             aes(x=NMDS1, y=NMDS2, fill=Treatment, color=Treatment),
             size=5, stroke=1.5, shape=21) +
  scale_color_manual(values = cost.col.line, name = "Treatment")+
  scale_fill_manual(values = cost.col.fill, name = "Treatment")+
  #ggtitle("Microbial Communities")+
  coord_fixed(ratio=1.9)+
  theme_bw()+
  ggforce::geom_mark_ellipse(data=nmds.16S.tend.scores1, aes(x=NMDS1, y=NMDS2, linetype=Stress_status_v1, label=Stress_status_v1), con.type="none", label.buffer=unit(4,'mm'), show.legend=F)+
  annotate("text", label="p > 0.001 \n stress = 0.082", x=-.225, y=-.125)
ggsave("NMDS_16S_tfinal.jpg", path = dirFigs, width = 7.5, height = 5.37, units = "in", dpi = 320)

#Test statistical effects using permanova.
#first test the effect of timepoint on all samples.
permanova.timepoint <- adonis2(unifrac.dist1~Timepoint, by="margin", permutations=999, data=metadata1)
permanova.timepoint

#then test for an effect of treatment at t0.
permanova.t0 <- adonis2(unifrac.dist.t0~Treatment, by="margin", permutations=999, data=metadata.t0)
permanova.t0 #Treatment1 is significant.

#now test tend data for all treatments.
permanova.tend <- adonis2(unifrac.dist.tend~Treatment, by="margin", permutations=999, data=metadata.tend)
permanova.tend #treatment is significant, with very high R2 and medium F values.
pairwise.adonis.tend <- pairwise.adonis(unifrac.dist.tend,factors=metadata.tend$Treatment) #no significant differences because sample size is too small

#now test for the effect of coral vs. water
metadata.tend$organism <- as.character(metadata.tend$Origin_PlanC) #duplicate
metadata.tend$organism[metadata.tend$organism!="control"] <- "coral"
permanova.tend.organism <- adonis2(unifrac.dist.tend ~ organism, by="margin", permutations=999, data=metadata.tend)
permanova.tend.organism #significant

#work up pairwise.adonis.tend for export
pairwise.adonis.tend1 <- cbind(pairwise.adonis.tend, t(as.data.frame(strsplit(as.character(pairwise.adonis.tend$pairs), "vs")))) #split the "pairs" into the two seperate treatments, add as new columns
colnames(pairwise.adonis.tend1)[9] <- "Treatment1" #update colnames
colnames(pairwise.adonis.tend1)[10] <- "Treatment2" #update colnames

#export pairwise adonis tend for visualization in JMP
write.csv(pairwise.adonis.tend1, file.path(dirOutput, "pairwise.adonis.tend1.csv"), )

#read in and work up pairwise.adonis.tend1.square.matrix.csv for visualization
pairwise.adonis.tend1.square.matrix <- read.csv(file.path(dirOutput, "pairwise.adonis.tend1.square.matrix.csv"))
rownames(pairwise.adonis.tend1.square.matrix) <- pairwise.adonis.tend1.square.matrix$X #update rownames
colnames(pairwise.adonis.tend1.square.matrix)[2:7] <- as.vector(pairwise.adonis.tend1.square.matrix$X) #update colnames
pairwise.adonis.tend1.square.matrix1 <- pairwise.adonis.tend1.square.matrix[,-1] #remove extra column

#visualize with pheatmap
jpeg("../figures/16S_pairwise_adonis_clust.jpg",width=2100, height=2000, res=300)
pheatmap(pairwise.adonis.tend1.square.matrix1, color=viridis(n=256, alpha = 1, begin = 0, end = 1, direction = 1, option="B")
)
dev.off()

#now test tend coral data for all treatments
permanova.coral.tend=adonis2(unifrac.dist.coral.tend~Treatment, by="margin", permutations=999, data=metadata.coral.tend)
permanova.coral.tend #treatment is significant.
pairwise.adonis(unifrac.dist.coral.tend,factors=metadata.coral.tend$Treatment) #no significant differences because sample size is too small

#Work up and analyze the tfinal adiv data.
adiv.tend <- merge(adiv, metadata.tend, by.x="group", by.y="Sample_Name_Unique") #merge with metadata
levels(adiv.tend$Treatment) <- levels(fact.all.treat) #reorder levels

#visualize
adiv.tend.sobs <- ggplot(adiv.tend, aes(x=Treatment, y=sobs, fill=Treatment, col=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  theme_classic()+
  theme(legend.position="none")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylab("Observed ASVs")+
  xlab("")+
  ggtitle("Observed ASVs")

adiv.tend.chao <- ggplot(adiv.tend, aes(x=Treatment, y=chao, fill=Treatment, col=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  theme_classic()+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylab("Chao Diversity")+
  xlab("")+
  ggtitle("Chao Diversity")

adiv.tend.shannon <- ggplot(adiv.tend, aes(x=Treatment, y=shannon, fill=Treatment, col=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  theme_classic()+
  theme(legend.position="none")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylab("Shannon Diversity")+
  xlab("")+
  ggtitle("Shannon Diversity")

adiv.tend.shannoneven <- ggplot(adiv.tend, aes(x=Treatment, y=shannoneven, fill=Treatment, col=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  theme_classic()+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylab("Shannon's Eveness")+
  xlab("")+
  ggtitle("Shannon's Eveness")

#export
jpeg("../figures/16S_tend_alpha_diversity.jpg",width=4900, height=2800, res=300)
plot_grid(adiv.tend.sobs, adiv.tend.chao, adiv.tend.shannon, adiv.tend.shannoneven, nrow=2, ncol=2, rel_widths=c(1,1.3))
dev.off()

#run stats on adiv data
#first visualize distributions
hist(adiv.tend$sobs) #normal ish
hist(adiv.tend$chao) #more normal
hist(adiv.tend$shannon) #more normal
hist(adiv.tend$shannoneven) #normal

#run both anova and k-w test on sobs data
mod.sobs <- aov(sobs ~ Treatment, data=adiv.tend) #run anova
summary(mod.sobs) #not significant
mod.sobs.kw <- kruskal.test(sobs ~ Treatment, data=adiv.tend)
mod.sobs.kw #not significant

#run anova on chao data
mod.chao <- aov(chao ~ Treatment, data=adiv.tend) #run anova
summary(mod.chao) #not significant
mod.kw.chao <- kruskal.test(chao ~ Treatment, data=adiv.tend)
mod.kw.chao #not significant

#run anova on shannon data
mod.shannon <- aov(shannon ~ Treatment, data=adiv.tend) #run anova
summary(mod.shannon) #significant
TukeyHSD(mod.shannon, "Treatment") #run post hoc. no significant pairwise differences.

#run anova on shannoneven data
mod.shannoneven <- aov(shannoneven ~ Treatment, data=adiv.tend) #run anova
summary(mod.shannoneven) #not significant
mod.shannoneven.kw <- kruskal.test(shannoneven ~ Treatment, data=adiv.tend)
mod.shannoneven.kw #not significant
