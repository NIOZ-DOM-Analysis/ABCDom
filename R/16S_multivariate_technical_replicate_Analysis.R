#load libraries
library(vegan)
library(ggplot2)

#load custom functions
source(file="generate.square.dist.script.07.27.2020.R") #need to update this

#load metadata and data
unifrac <- read.csv(file.path(dirRAW, "16S/ASRAAMP_MCR2019_16S_Mar_2021_2k_subsample_04212021", "otu_repr_100.tre1.wsummary.csv"))
metadata <- read.csv(file.path(dirRAW, "16S/ASRAAMP_MCR2019_16S_Mar_2021_2k_subsample_04212021", "metadata_tech_replicates.csv"))
metadata1 <- read.csv(file.path(dirRAW, "16S/ASRAAMP_MCR2019_16S_Mar_2021_2k_subsample_04212021", "metadata1.csv"))

#Work up the abcDOM tech replicate metadata.
#subset for just abc and pcr pos/neg samples
metadata_abc_replicate=metadata[metadata$Experiment!="ROBO" & metadata$Experiment!="ASRAMP" & metadata$Experiment!="UNKNOWN",]

#subset for just abc samples
metadata_abc_replicate1=metadata_abc_replicate[metadata_abc_replicate$Experiment=="ABC",]

#merge metadata_abc_replicate1 with metadata1
metadata_abc_replicate2=merge(metadata_abc_replicate1,metadata1,by.x="Sample_Name",by.y="Sample_NR",all.x=T,all.y=F)

#subset for just tend
metadata_abc_replicate_tend=metadata_abc_replicate2[metadata_abc_replicate2$Timepoint.x=="Tend",]
metadata_abc_replicate_tend=metadata_abc_replicate_tend[-35:-36,]

#Work up unifrac data.
#work up abc_replicate data
generate.square.dist.matrix(unifrac,metadata_abc_replicate,1)

#store outputs in new dfs
abc_replicate_unifrac.dist.matrix=dist.matrix1
abc_replicate_unifrac.dist=square.dist.matrix

#work up abc_replicate1 data
generate.square.dist.matrix(unifrac,metadata_abc_replicate1,1)

#store outputs in new dfs
abc_replicate1_unifrac.dist.matrix=dist.matrix1
abc_replicate1_unifrac.dist=square.dist.matrix

#work up abc_replicate_tend data
generate.square.dist.matrix(unifrac,metadata_abc_replicate_tend,2)
abc_replicate_tend_unifrac.dist.matrix=dist.matrix1
abc_replicate_tend_unifrac.dist=square.dist.matrix

#Work up NMDS objects from the unifrac square dist matrices.
nmds.abc.replicate=metaMDS(abc_replicate_unifrac.dist,k=2,trymax=100) #very low stress, probably cause PCR controls are leveraging the ordination

nmds.abc.replicate1=metaMDS(abc_replicate1_unifrac.dist,k=2,trymax=100) #still pretty low stress

nmds.abc.tend=metaMDS(abc_replicate_tend_unifrac.dist,k=2,trymax=100)

#Visualize, coloring ordinations by sample type.
colvec=c("blue","green","red","yellow","purple")

plot(nmds.abc.replicate,type="n")
points(nmds.abc.replicate,display="sites",
       col=colvec[metadata_abc_replicate$Experiment])
legend("topleft",legend=c(levels(metadata_abc_replicate$Experiment)),pch=c(16,16,16,16),col=colvec,bty="n",cex=.5,y.intersp = .5)
#PCR positive and negative controls (blue) cluster far away from abcDOM samples (green).

#Now, visualize just the abcDOM samples, colored by replicate with replicates connected by an arrow.
colvec2=c("blue","red")

plot(nmds.abc.replicate1,type="n")
points(nmds.abc.replicate1,display="sites",
       col=colvec2[metadata_abc_replicate1$Replicate],
       pch=16)
ordiarrows(nmds.abc.replicate1,groups=metadata_abc_replicate1$Sample_Name)
#The left hand cluster corresponds to the T0 samples, which seem to group nicely. The right hand cluster corresponds to the Tfinal samples, for which some technical replicates do not group together.

#Test the effect of replicate on community structure.
permanova.replicate=adonis2(abc_replicate1_unifrac.dist~Replicate,data=metadata_abc_replicate1,by="margin",permutations=999)
permanova.replicate #No effect and very low R2 for technical replicates.

#Now visualize just the tend samples. Save output for potential supplement.
colvec3=c("orange","red","yellow","green","gray","black")
pointvec=c(16,17)


jpeg("../figures/16S_technical_replicate_tfinal_nmds.jpg",width=2100, height=1500, res=300)
plot(nmds.abc.tend,type="n")
points(nmds.abc.tend,display="sites",
       col=colvec3[metadata_abc_replicate_tend$Treatment],
       pch=pointvec[metadata_abc_replicate_tend$PCR_setting],
       cex=2)
ordiarrows(nmds.abc.tend,groups=metadata_abc_replicate_tend$Sample_Name)
ordilabel(nmds.abc.tend,labels=metadata_abc_replicate_tend$Sample_Name_Unique,border=NA,fill=NA,cex=.5)
dev.off()
#Looks like tech replicates for ABC_055 and ABC_059 did not cluster well together. When binning tech replicates during demux, be sure that these 4 samples are kept seperate.

#justify definition of outliers with stats
lower_tri_abc_replicate_tend_unifrac.dist.matrix <- abc_replicate_tend_unifrac.dist.matrix[lower.tri(abc_replicate_tend_unifrac.dist.matrix, diag=T)] #extract lower triangle
hist(lower_tri_abc_replicate_tend_unifrac.dist.matrix) #visualize histogram. Twou outliers on the far right?
sum(lower_tri_abc_replicate_tend_unifrac.dist.matrix >= .75) #confirmed two outliers above unifrac=.75
sd_lower_tri <- sd(lower_tri_abc_replicate_tend_unifrac.dist.matrix) #calculate SD
mean_lower_tri <- mean(lower_tri_abc_replicate_tend_unifrac.dist.matrix) #calculate mean
mean_lower_tri + 3*sd_lower_tri #ctoff will be average unifrac dist between tech replicates + 3*SD. This cutoff is .7227101
