#load libraries
library(ggplot2)
library(vegan)
library(stats)
library(cowplot)
library(dendextend)
library(scales)
library(FSA)
library(graphics)
library(grDevices)
library(pheatmap)
library(viridis)
library(ggforce)

#-----------------------------------------------------------------------------------------------------------------------------
#read in metadata and data
unifrac <- read.csv(file.path(dirRAW, "16S", "otu_repr_100.tre1.wsummary.csv")) #Read in unifrac.

#read in ASV table
abund <- read.csv(file.path(dirRAW, "16S", "abundance_table_100.database.csv"))

adiv <- read.csv(file.path(dirRAW, "16S", "all_alphaDiversity_100.summary.csv")) #read in alpha diversity data

#-----------------------------------------------------------------------------------------------------------------------------
#set global settings for visualization
cost.col.fill<-c("dodgerblue1","firebrick1", "white", "white", "grey70", "grey70") #set fill vec

cost.col.line<-c("dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3") #set outline vec

#-----------------------------------------------------------------------------------------------------------------------------
#work up metadata
metadata_16S_1 <- metadata_abc_replicates2[metadata_abc_replicates2$Sample_Name_Unique %in% colnames(abund),] #subset metadata 16S for these samples.

#Remove outliers. Outlier tech replicates defined in 16S_multivariate_technical_replicate_Analysis.R. Include the technical replicates that are closer to the centroid of their treatment.
metadata_16S_outliers <- metadata_16S_1[metadata_16S_1$Sample_Name_Unique!="ABC_055" & metadata_16S_1$Sample_Name_Unique!="ABC_059" & metadata_16S_1$Sample_Name_Unique!="ABC_MOC_B",] #Remove samples ABC_055 and ABC_059, include the v1 versions. Remove ABC_MOC_B because that is not a part of this experiment.

metadata_16S_outliers$Treatment_Timepoint <- paste(metadata_16S_outliers$Treatment, metadata_16S_outliers$Timepoint, sep="_") #generate a new factor that combines Treatment1 and Timepoint.

#-----------------------------------------------------------------------------------------------------------------------------
#Now define and remove potential biological  outliers.
generate.square.dist.matrix(unifrac, metadata_16S_outliers, 2) #generate unifrac dist matrix
unifrac.dist.matrix.outliers <- dist.matrix1 #store output as new data frame
unifrac.dist.outliers <- as.dist(unifrac.dist.matrix.outliers) #store output as new dist object

dispersion <- betadisper(unifrac.dist.outliers, group=metadata_16S_outliers$Treatment_Timepoint, type="centroid") #define within group variability.

boxplot(dispersion) #visualize dispersion

dist_centroid <- as.data.frame(dispersion$distances) #store distance to centroid values in dataframe

dist_centroid$sample <- rownames(dist_centroid) #add a sample column

dist_centroid1 <- merge(dist_centroid, metadata_16S_outliers, by.x="sample", by.y="Sample_Name_Unique") #merge in metadata

dist_centroid1$log10_dist <- log10(dist_centroid1$`dispersion$distances`) #add a log10 transformed dis centroid variable to df.

mean(log10(dist_centroid$`dispersion$distances`))+1.5*sd(log10(dist_centroid$`dispersion$distances`)) #Outliers defined as log10(dist) value greater than 1.5*sd(log10(dist))+mean(log10(dist)) = -.627
#This defines ABC072 and 064 as outliers.

metadata_16S_2 <- metadata_16S_outliers[c(-26,-34),] #remove outliers

metadata_16S_2$Treatment <- factor(metadata_16S_2$Treatment, levels=c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated")) #make treatment a factor

#generate a new stress status column
metadata_16S_2$stress_status_v1 <- as.character(metadata_16S_2$Treatment) #duplicate treatment column
metadata_16S_2$stress_status_v1[metadata_16S_2$stress_status_v1=="Control"] <- "Ambient"
metadata_16S_2$stress_status_v1[metadata_16S_2$stress_status_v1=="Negative Control" | metadata_16S_2$stress_status_v1=="Negative Control + Heated"] <- NA
metadata_16S_2$stress_status_v1[metadata_16S_2$stress_status_v1!="Ambient" & is.na(metadata_16S_2$stress_status_v1)!=T] <- "Stressed"

#-----------------------------------------------------------------------------------------------------------------------------
#subset metadata for various treatments/timepoints.
metadata.t0 <- metadata_16S_2[metadata_16S_2$Timepoint_char=="T0",] #subset for only t0 samples

metadata.tend <- metadata_16S_2[metadata_16S_2$Timepoint_char=="Tend",] #subset for only the tend samples

metadata.coral.tend <- subset(metadata.tend, Treatment!="Negative Control" & Treatment!="Negative Control + Heated") #subset for only tend coral samples.

#-----------------------------------------------------------------------------------------------------------------------------
#Work up unifrac.
generate.square.dist.matrix(unifrac, metadata_16S_2, 2) #work up unifrac for metadata_16S_2
unifrac.dist.matrix1 <- dist.matrix1 #store output as new data frame
unifrac.dist1 <- as.dist(unifrac.dist.matrix1) #store output as new dist object

generate.square.dist.matrix(unifrac, metadata.t0, 2) #work up unifrac for metadata.t0
unifrac.dist.matrix.t0 <- dist.matrix1 #store output as new data frame
unifrac.dist.t0 <- as.dist(unifrac.dist.matrix.t0) #store output as new dist object

generate.square.dist.matrix(unifrac,metadata.tend,2) #work up unifrac for metadata.tend
unifrac.dist.matrix.tend <- dist.matrix1
unifrac.dist.tend <- as.dist(unifrac.dist.matrix.tend)

generate.square.dist.matrix(unifrac, metadata.coral.tend, 2) #work up unifrac for metadata.coral.tend
unifrac.dist.matrix.coral.tend <- dist.matrix1
unifrac.dist.coral.tend <- as.dist(unifrac.dist.matrix.coral.tend)

#-----------------------------------------------------------------------------------------------------------------------------
#plot nmds of t0 and tend samples for supplemental figure 4
#First workup nmds
nmds.16S <- metaMDS(unifrac.dist1, k=2, trymax=100) #generate nmds. stress=.04

nmds.16S.scores <- as.data.frame(scores(nmds.16S)) #extract scores

nmds.16S.scores$Sample <- rownames(nmds.16S.scores) #add sample names column

nmds.16S.scores1 <- merge(nmds.16S.scores, metadata_16S_2, by.x="Sample", by.y="Sample_Name_Unique", all.x=T, all.y=F) #merge with metadata

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
ggsave("FigS4.jpg", path = dirFigs, width = 6.75, height = 5, units = "in", dpi = 320)

write.csv(nmds.16S.scores1, file=file.path(dirOutput, "FigS4_data.csv"))

#-----------------------------------------------------------------------------------------------------------------------------
#plot nmds of tend samples for supplement
#First workup nmds
nmds.16S.tend <- metaMDS(unifrac.dist.tend, k=2, trymax=100) #generate nmds. stress=.09

nmds.16S.tend.scores <- as.data.frame(scores(nmds.16S.tend)) #extract scores

nmds.16S.tend.scores$Sample <- rownames(nmds.16S.tend.scores) #add sample names column

nmds.16S.tend.scores1 <- merge(nmds.16S.tend.scores, metadata_16S_2, by.x="Sample", by.y="Sample_Name_Unique", all.x=T, all.y=F) #merge with metadata

#plot nmds. Store to be used as figure 3A.
fig3A<-ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data=nmds.16S.tend.scores1,
             aes(x=NMDS1, y=NMDS2, fill=Treatment, color=Treatment),
             size=5, stroke=1.5, shape=21) +
  scale_color_manual(values = cost.col.line, name = "Treatment", labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  scale_fill_manual(values = cost.col.fill, name = "Treatment", labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  #ggtitle("Microbial Communities")+
  coord_fixed(ratio=1.9)+
  theme_bw()+
  ggforce::geom_mark_ellipse(data=nmds.16S.tend.scores1, aes(x=NMDS1, y=NMDS2, linetype=stress_status_v1, label=stress_status_v1), con.type="none", label.buffer=unit(4,'mm'), show.legend=F)+
  annotate("text", label="p > 0.001 \n stress = 0.082", x=-.225, y=-.125)
ggsave("fig3A.jpg", plot = fig3A, path = dirFigs, dpi  = 300, width = 8, height = 6)

#NEED TO FIGURE OUT DIMENSIONSggsave("figure_3A.png", )

#export data used in figure 3a
write_csv(nmds.16S.tend.scores1, file=file.path(dirOutput, "Fig3a_data.csv"))

#-----------------------------------------------------------------------------------------------------------------------------
#Test statistical effects using permanova.
permanova.timepoint <- adonis2(unifrac.dist1~Timepoint, by="margin", permutations=999, data=metadata_16S_2) #first test the effect of timepoint on all samples.
permanova.timepoint

permanova.tend <- adonis2(unifrac.dist.tend~Treatment, by="margin", permutations=999, data=metadata.tend) #now test the effect of all treatments on the tend data.
permanova.tend

#now test for the effect of coral vs. water on the tend data.
metadata.tend$organism <- as.character(metadata.tend$Treatment) #duplicate

metadata.tend$organism[metadata.tend$organism=="Non-bleached + Ambient" | metadata.tend$organism=="Non-bleached + Heated" | metadata.tend$organism=="Bleached + Ambient" | metadata.tend$organism=="Bleached + Heated"] <- "coral" #replace coral treatments with "coral" character

metadata.tend$organism[metadata.tend$organism!="coral"] <- "water" #replace non-coral treatments with "water" character

permanova.tend.organism <- adonis2(unifrac.dist.tend ~ organism, by="margin", permutations=999, data=metadata.tend) #run the permanova
permanova.tend.organism

permanova.coral.tend <- adonis2(unifrac.dist.coral.tend~Stress_status, by="margin", permutations=999, data=metadata.coral.tend) #now test the effect of "stress status" on just the tend coral data
permanova.coral.tend

#-----------------------------------------------------------------------------------------------------------------------------
#Work up and analyze the tfinal adiv data.
adiv.tend <- merge(adiv, metadata.tend, by.x="group", by.y="Sample_Name_Unique") #merge with metadata

levels(adiv.tend$Treatment) <- levels(fact.all.treat) #reorder levels

#visualize
adiv.tend.sobs <- ggplot(adiv.tend, aes(x=Treatment, y=sobs, fill=Treatment, col=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  scale_color_manual(values=cost.col.line, labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control\n+ Heated"))+
  scale_fill_manual(values=cost.col.fill, labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control\n+ Heated"))+
  theme_classic()+
  theme(legend.position="none")+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control\n+ Heated"))+
  ylab("Observed ASVs")+
  xlab("")+
  ggtitle("Observed ASVs")

adiv.tend.chao <- ggplot(adiv.tend, aes(x=Treatment, y=chao, fill=Treatment, col=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  scale_color_manual(values=cost.col.line, labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control\n+ Heated"))+
  scale_fill_manual(values=cost.col.fill, labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control\n+ Heated"))+
  theme_classic()+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control\n+ Heated"))+
  ylab("Chao Diversity")+
  xlab("")+
  ggtitle("Chao Diversity")

adiv.tend.shannon <- ggplot(adiv.tend, aes(x=Treatment, y=shannon, fill=Treatment, col=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  scale_color_manual(values=cost.col.line, labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control\n+ Heated"))+
  scale_fill_manual(values=cost.col.fill, labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control\n+ Heated"))+
  theme_classic()+
  theme(legend.position="none")+
  scale_x_discrete( labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control\n+ Heated"))+
  ylab("Shannon Diversity")+
  xlab("")+
  ggtitle("Shannon Diversity")

adiv.tend.shannoneven <- ggplot(adiv.tend, aes(x=Treatment, y=shannoneven, fill=Treatment, col=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  scale_color_manual(values=cost.col.line, labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control\n+ Heated"))+
  scale_fill_manual(values=cost.col.fill, labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control\n+ Heated"))+
  theme_classic()+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control\n+ Heated"))+
  ylab("Shannon's Eveness")+
  xlab("")+
  ggtitle("Shannon's Eveness")

#export for supplemental figure 5
#jpeg(file=file.path(dirFigs, "Fig_S5.jpg"),width=5000, height=2800, res=300)
#plot_grid(adiv.tend.sobs, adiv.tend.chao, adiv.tend.shannon, adiv.tend.shannoneven, nrow=2, ncol=2, rel_widths=c(1,1.3), align = c("hv"), axis = "lt" )
#dev.off()

#export supplemental figure 5 data
#write.csv(adiv.tend, file=file.path(dirOutput, "FigS5_data.csv"))

#-----------------------------------------------------------------------------------------------------------------------------
#run stats on adiv data
#first visualize distributions
hist(adiv.tend$sobs)
hist(adiv.tend$chao)
hist(adiv.tend$shannon)
hist(adiv.tend$shannoneven)

#run both anova and k-w test on sobs data
mod.sobs <- aov(sobs ~ Treatment, data=adiv.tend) #run anova
summary(mod.sobs) #not significant

mod.sobs.kw <- kruskal.test(sobs ~ Treatment, data=adiv.tend) #run k-w
mod.sobs.kw #not significant

#run anova on chao data
mod.chao <- aov(chao ~ Treatment, data=adiv.tend) #run anova
summary(mod.chao) #not significant

#run anova on shannon data
mod.shannon <- aov(shannon ~ Treatment, data=adiv.tend) #run anova
summary(mod.shannon) #significant

TukeyHSD(mod.shannon, "Treatment") #run post hoc. no significant pairwise differences.

#run anova on shannoneven data
mod.shannoneven <- aov(shannoneven ~ Treatment, data=adiv.tend) #run anova
summary(mod.shannoneven) #not significant
