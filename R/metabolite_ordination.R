# #non metric multidimentional scaling, ABCT0Dom T0
#ordination of data
library(BiodiversityR)
library(openxlsx)
library(dendextend)
library(vegan)
library(pairwiseAdonis)
library(pheatmap)
library(viridis)
library(dplyr)
library(ggplot2)

df.area.ABCT0 <- df.area %>% filter(Timepoint_char == "T0") %>%
   filter(Treatment != "Inoculum") #  %>%
  # filter(`Sample Name` != "ABC_022") %>%
  # filter(`Sample Name` != "PC_Blank")
write_csv(df.area.ABCT0, paste0(dirOutput, "/df.area.ABCT0.csv"))

# #lets also do this with the normalized data
# df.norm.area.ABCT0 <- df.norm.area %>% filter(Timepoint_char == "T0") %>%
#   filter(Treatment != "Inoculum") %>%
#   filter(`Sample Name` != "ABC_022") %>%
#   filter(`Sample Name` != "PC_Blank")
# write_csv(df.norm.area.ABCT0, paste0(dirOutput, "/df.norm.area.ABCT0.csv"))
#
# df.norm.smpl.ABCT0 <- df.norm.smpl %>% filter(Timepoint_char == "T0") %>%
#   filter(Treatment != "Inoculum") %>%
#   filter(`Sample Name` != "ABC_022") %>%
#   filter(`Sample Name` != "PC_Blank")
# write_csv(df.norm.smpl.ABCT0, paste0(dirOutput, "/df.norm.smpl.ABCT0.csv"))

#ordination df.area

ord.mod.area<- metaMDS(df.area.ABCT0[M:ncol(df.area.ABCT0)], distance = 'bray', k = 2)
NMDS.ABCDom.T0<- bind_cols(df.area.ABCT0[1:(M-1)], as.data.frame(ord.mod.area$points))
#extract stress
stress <- round(ord.mod.area$stress, 4)
stress <- paste0("ABCDom T0 bray curtis dissimilarity k=2 \nstress = ", stress)

cost.col.fill<-c("dodgerblue1","firebrick1", "white", "white", "grey70", "grey70")
cost.col.line<-c("dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3")
fact.all.treat<-factor(NMDS.ABCDom.T0$Treatment, levels = c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", 'Bleached + Heated', "Ambient Water Control", "Heated Water Control"))
treat.labels <- c("Control", "Heated", "Bleached", 'Bleached + Heated', "Negative Control", "Negative Control + Heated")

#add a stress_status_v1 column to NMDS.ABCDom.T0
NMDS.ABCDom.T0$stress_status_v1 <- as.character(NMDS.ABCDom.T0$Treatment) #duplicate treatment
NMDS.ABCDom.T0$stress_status_v1[NMDS.ABCDom.T0$stress_status_v1 == "Non-bleached + Ambient"] <- "Control"
NMDS.ABCDom.T0$stress_status_v1[NMDS.ABCDom.T0$stress_status_v1 == "Ambient Water Control" |  NMDS.ABCDom.T0$stress_status_v1 == "Heated Water Control"] <- NA
NMDS.ABCDom.T0$stress_status_v1[NMDS.ABCDom.T0$stress_status_v1 != "Control" & is.na(NMDS.ABCDom.T0$stress_status_v1) == FALSE] <- "Stressed"

ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data=NMDS.ABCDom.T0,
             aes(x=MDS1, y=MDS2, fill = fact.all.treat, colour=fact.all.treat),
             size=5, stroke=1.5, shape = 21) +
  scale_color_manual(labels = treat.labels, values = cost.col.line, name = "Treatment")+
  scale_fill_manual(labels = treat.labels, values = cost.col.fill, name = "Treatment", guide = guide_legend(override.aes = list(shape = 21)))+
  #ggtitle(stress)+
  annotate("text", label="p > 0.001 \n stress = 0.0467", x=-.11, y=.05)+
  theme_bw()+
  coord_fixed(ratio=1.2)+
  ggforce::geom_mark_ellipse(data=NMDS.ABCDom.T0, aes(x=MDS1, y=MDS2, linetype=stress_status_v1, label=stress_status_v1), con.type="none", label.buffer=unit(4,'mm'), show.legend=F)
ggsave("NMDS_ABCDom_T0.jpg", path = dirFigs, width = 6.75, height = 5, units = "in", dpi = 320)

fig4A<-ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data=NMDS.ABCDom.T0,
             aes(x=MDS1, y=MDS2, fill = fact.all.treat, colour=fact.all.treat),
             size=5, stroke=1.5, shape = 21) +
  scale_color_manual(labels = treat.labels, values = cost.col.line, name = "Treatment")+
  scale_fill_manual(labels = treat.labels, values = cost.col.fill, name = "Treatment", guide = guide_legend(override.aes = list(shape = 21)))+
  #ggtitle(stress)+
  annotate("text", label="p > 0.001 \n stress = 0.0467", x=-.11, y=.05)+
  theme_bw()+
  theme(axis.ticks.x.top = element_blank(), axis.ticks.y.right = element_blank())+
  coord_fixed(ratio=1.2)+
  ggforce::geom_mark_ellipse(data=NMDS.ABCDom.T0, aes(x=MDS1, y=MDS2, linetype=stress_status_v1, label=stress_status_v1), con.type="none", label.buffer=unit(4,'mm'), show.legend=F)
fig4A


#visualize as clusterding dendrogram
df.area.ABCT0.bray <- vegdist(df.area.ABCT0[,M:ncol(df.area.ABCT0)], method="bray") #generate bray curtis dist matrix

clust.metabolites.t0 <- hclust(df.area.ABCT0.bray,method="ward.D") #cluster t0 samples
plot(clust.metabolites.t0)
clust.metabolites.t0.1 <- rotate(clust.metabolites.t0, order=c("15","17","16","8","7","6", "12", "11", "4", "9", "14","13", "5", "10", "1","2","3")) #rotate leaves of dendogram
plot(clust.metabolites.t0.1)
clust.metabolites.t0.1$labels <- df.area.ABCT0$Treatment #adjust labels
plot(clust.metabolites.t0.1)




#plot heirarchical clustering plot, save as jpeg
dev.off()
df.area.ABCT0$Treatment <- factor(df.area.ABCT0$Treatment, levels=levels(fact.all.treat)) #make treatment a factor
jpeg("../figures/metabolomics_dendogram.jpeg",width=2100, height=2100, res=300)
clust.metabolites.t0.2 <- as.dendrogram(clust.metabolites.t0.1, hang=-1) #convert to dendrogram
plot(clust.metabolites.t0.2, type="rectangle", dLeaf=.3, ylim=c(-1.5,2.5))
symbols(1:17, rep(-.15,17), circles=rep(1,17), add=T, fg=cost.col.line[df.area.ABCT0$Treatment][clust.metabolites.t0.1$order], bg=cost.col.fill[df.area.ABCT0$Treatment][clust.metabolites.t0.1$order], inches=.09, xpd=T, lwd=3)


# plot legend (separately) for hierarchical clustering plot and export as jpeg. Can just reuse 16S dendrogram legend.

jpeg("../figures/metabolomics_dendogram_legend.jpg", width=1500, height=2100, res=300)
plot(NULL, axes=F, bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(0, .5, legend=levels(df.area.ABCT0$Treatment), col=cost.col.line, pt.bg=cost.col.fill, pch=21, pt.cex=2, pt.lwd=2.75, cex=.75, y.intersp=1.5, bty='n')
dev.off()


#lets do stats
adon.results<-adonis(df.area.ABCT0[,M:ncol(df.area.ABCT0)] ~ fact.all.treat, method="bray",perm=999)
print(adon.results)

adon2.results<-adonis2(df.area.ABCT0[,M:ncol(df.area.ABCT0)] ~ fact.all.treat, method="bray",perm=999)
print(adon2.results)

fact.coral.treat <- fact.all.treat[fact.all.treat!="Ambient Water Control" & fact.all.treat!="Heated Water Control"] #subset for just coral treatments
adon2.coral.results <- adonis2(df.area.ABCT0[df.area.ABCT0$Origin_PlanC!="control",M:ncol(df.area.ABCT0)] ~ fact.coral.treat, method="bray",perm=999)
print(adon2.coral.results)


#make 3 groups
fact.coral.control.treat <- fact.all.treat
levels(fact.coral.control.treat)[levels(fact.coral.control.treat) == "Ambient Water Control"] <- "Control"
levels(fact.coral.control.treat)[levels(fact.coral.control.treat) == "Heated Water Control"] <- "Control"
levels(fact.coral.control.treat)[levels(fact.coral.control.treat) == "Non-bleached + Ambient"]<- "Ambient"
levels(fact.coral.control.treat)[levels(fact.coral.control.treat) != "Ambient" & levels(fact.coral.control.treat) != "Control"]<- "Stressed"
fact.coral.control.treat
adon2.coral.control.results <- adonis2(df.area.ABCT0[,M:ncol(df.area.ABCT0)] ~ fact.coral.control.treat, method="bray",perm=999)
print(adon2.coral.control.results)





#run pairwise adonis
pairwise.adonis.mb.t0 <- pairwise.adonis(df.area.ABCT0.bray, factors=df.area.ABCT0$Treatment)
write.csv(pairwise.adonis.mb.t0, file.path(dirOutput, "pairwise.adonis.mb.t0.csv"), ) #export

#read in square dist matrix and work up for visualization
pairwise.adonis.mb.t0.square.matrix <- read.csv(file.path(dirOutput, "pairwise.adonis.mb.t0.square.matrix.csv"))
rownames(pairwise.adonis.mb.t0.square.matrix) <- pairwise.adonis.mb.t0.square.matrix$X #update rownames
pairwise.adonis.mb.t0.square.matrix1 <- pairwise.adonis.mb.t0.square.matrix[,-1] #remove extra column
colnames(pairwise.adonis.mb.t0.square.matrix1) <- rownames(pairwise.adonis.mb.t0.square.matrix1) #update colnames

#visualize
dev.off()
jpeg("../figures/MB_pairwise_adonis_clust.jpeg",width=2100, height=2000, res=300)
pheatmap(pairwise.adonis.mb.t0.square.matrix1, color=viridis(n=256, alpha = 1, begin = 0, end = 1, direction = 1, option="B"))
dev.off()


# dis <- vegdist(df.area.ABCT0[,M:ncol(df.area.ABCT0)])
# mod <- betadisper(dis, df.area.ABCT0$Treatment)
# mod
# # extract the centroids and the site points in multivariate space.
# centroids<-data.frame(grps=rownames(mod$centroids),data.frame(mod$centroids))
# vectors<-data.frame(group=mod$group,data.frame(mod$vectors))
#
# # to create the lines from the centroids to each point we will put it in a format that ggplot can handle
# seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
# names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
# grp1.hull<-seg.data[seg.data$group=="Non-bleached + Ambient",1:3][chull(seg.data[seg.data$group=="Non-bleached + Ambient",2:3]),]
# grp2.hull<-seg.data[seg.data$group=="Non-bleached + Heated",1:3][chull(seg.data[seg.data$group=="Non-bleached + Heated",2:3]),]
# grp3.hull<-seg.data[seg.data$group=="Bleached + Ambient",1:3][chull(seg.data[seg.data$group=="Bleached + Ambient",2:3]),]
# grp4.hull<-seg.data[seg.data$group=="Bleached + Heated",1:3][chull(seg.data[seg.data$group=="Bleached + Heated",2:3]),]
# grp5.hull<-seg.data[seg.data$group=="Ambient Water Control",1:3][chull(seg.data[seg.data$group=="Ambient Water Control",2:3]),]
# grp6.hull<-seg.data[seg.data$group=="Heated Water Control",1:3][chull(seg.data[seg.data$group=="Heated Water Control",2:3]),]
# all.hull<-rbind(grp1.hull,grp2.hull,grp3.hull, grp4.hull, grp5.hull,grp6.hull)
#
#
# panel.d<-ggplot() +
#   geom_polygon(data=all.hull,aes(x=v.PCoA1,y=v.PCoA2),colour="black",alpha=0,linetype="dashed") +
#   geom_segment(data=seg.data,aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) +
#   geom_point(data=centroids[,1:3], aes(x=PCoA1,y=PCoA2,shape=grps),size=4,colour="red") +
#   geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2,shape=group),size=2) +
#   labs(title="All",x="",y="") +
#   theme(legend.position="none")
# panel.d

#####   not run because takes forever.
# # redundancy analyis ABCDomT0

# ABCDT0.hellinger <- disttransform(df.area.ABCT0[M:ncol(df.area.ABCT0)], method = 'hellinger')
# ord.mod.area2<-rda(ABCDT0.hellinger ~ Treatment, data = df.area.ABCT0[1:(M-1)], scaling = "species")
# summary(ord.mod.area2)
# ordi.plot2<-ordiplot(ord.mod.area2, choices = c(1,2))
# sites.long2<- sites.long(ordi.plot2, env.data = df.area.ABCT0[1:(M-1)])
# species.long2<-species.long(ordi.plot2)
# axis.long2<-axis.long(ord.mod.area2, choices = c(1,2))
#
# spec.envfit<-envfit(ordi.plot2, env = ABCDT0.hellinger, silent = TRUE)
# spec.data.envit <- data.frame(r=spec.envfit$vector$r, p=spec.envfit$vectors$pvals)
# species.long2 <- species.long(ordi.plot2, spec.data=spec.data.envit)
# species.long3 <- species.long2[species.long2$r >= 0.7, ]
# Redun.ABCDom.T0 <- sites.long2
# drivers.redund.ABCDom.T0 <- species.long2
#
# ggplot()+
#   geom_vline(xintercept = c(0), color = "grey70", linetype = 2)+
#   geom_hline(yintercept = c(0), color = "grey70", linetype = 2)+
#   xlab(axis.long2[1, "label"])+
#   ylab(axis.long2[2, "label"])+
#   scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL))+
#   scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL))+
#   geom_point(data=Redun.ABCDom.T0, aes(x=axis1, y=axis2, fill = fact.all.treat, colour=fact.all.treat),
#              size=5, stroke=1.5, shape = 21)+
#   geom_segment(data=drivers.redund.ABCDom.T0[drivers.redund.ABCDom.T0$r >= 0.8,],
#                aes(x=0, y=0, xend=axis1*4, yend=axis2*4), colour="red", size=0.7, arrow=arrow())+
#   scale_shape_manual(values = 21)+
#   scale_color_manual(labels = treat.labels, values = cost.col.line, name = "Treatment")+
#   scale_fill_manual(labels = treat.labels, values = cost.col.fill, name = "Treatment", guide = guide_legend(override.aes = list(shape = 21)))+
#   ggtitle("ABCDom T0, r>0.8")+
#   theme_bw()+
#   coord_fixed(ratio=1)
# ggsave("redundancy_ABCDom_T0_Treatment.jpg", path = dirFigs, width = 6.75, height = 5, units = "in", dpi = 320)
#
