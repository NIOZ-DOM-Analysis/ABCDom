# #non metric multidimentional scaling, ABCT0Dom T0
#ordination of data
library(BiodiversityR)
library(openxlsx)

df.area.ABCT0 <- df.area %>% filter(Experiment == "ABCDOM", Timepoint_char == "T0") %>% filter(Treatment != "Inoculum")
ord.mod.area<- metaMDS(df.area.ABCT0[M:ncol(df.area.ABCT0)], distance = 'bray', k = 2)
NMDS.ABCDom.T0<- bind_cols(df.area.ABCT0[1:(M-1)], as.data.frame(ord.mod.area$points))
#extract stress
stress <- round(ord.mod.area$stress, 4)
stress <- paste0("ABCDom T0 bray curtis dissimilarity k=2 \nstress = ", stress)

cost.col.fill<-c("blue","red", "white", "white", "grey70", "grey70")
cost.col.line<-c("dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3")
fact.all.treat<-factor(NMDS.ABCDom.T0$Treatment, levels = c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", 'Bleached + Heated', "Ambient Water Control", "Heated Water Control"))
treat.labels <- c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", 'Bleached + Heated', "Ambient Water Control", "Heated Water Control")

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
  ggtitle(stress)+
  coord_fixed(ratio=1)
ggsave("NMDS_ABCDom_T0_axis1vs2.jpg", path = dirFigs, width = 6.75, height = 5, units = "in", dpi = 320)

#lets do stats
adon.results<-adonis(df.area.ABCT0[,M:ncol(df.area.ABCT0)] ~ fact.all.treat, method="bray",perm=999)
print(adon.results)

adon2.results<-adonis2(df.area.ABCT0[,M:ncol(df.area.ABCT0)] ~ fact.all.treat, method="bray",perm=999)
print(adon2.results)

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

# # redundancy analyis ABCDomT0
ABCDT0.hellinger <- disttransform(df.area.ABCT0[M:ncol(df.area.ABCT0)], method = 'hellinger')
ord.mod.area2<-rda(ABCDT0.hellinger ~ Treatment, data = df.area.ABCT0[1:(M-1)], scaling = "species")
summary(ord.mod.area2)
ordi.plot2<-ordiplot(ord.mod.area2, choices = c(1,2))
sites.long2<- sites.long(ordi.plot2, env.data = df.area.ABCT0[1:(M-1)])
species.long2<-species.long(ordi.plot2)
axis.long2<-axis.long(ord.mod.area2, choices = c(1,2))

spec.envfit<-envfit(ordi.plot2, env = ABCDT0.hellinger, silent = TRUE)
spec.data.envit <- data.frame(r=spec.envfit$vector$r, p=spec.envfit$vectors$pvals)
species.long2 <- species.long(ordi.plot2, spec.data=spec.data.envit)
species.long3 <- species.long2[species.long2$r >= 0.7, ]
Redun.ABCDom.T0 <- sites.long2
drivers.redund.ABCDom.T0 <- species.long2

ggplot()+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2)+
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2)+
  xlab(axis.long2[1, "label"])+
  ylab(axis.long2[2, "label"])+
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL))+
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL))+
  geom_point(data=Redun.ABCDom.T0, aes(x=axis1, y=axis2, fill = fact.all.treat, colour=fact.all.treat),
             size=5, stroke=1.5, shape = 21)+
  geom_segment(data=drivers.redund.ABCDom.T0[drivers.redund.ABCDom.T0$r >= 0.8,],
               aes(x=0, y=0, xend=axis1*4, yend=axis2*4), colour="red", size=0.7, arrow=arrow())+
  scale_shape_manual(values = 21)+
  scale_color_manual(labels = treat.labels, values = cost.col.line, name = "Treatment")+
  scale_fill_manual(labels = treat.labels, values = cost.col.fill, name = "Treatment", guide = guide_legend(override.aes = list(shape = 21)))+
  ggtitle("ABCDom T0, r>0.8")+
  coord_fixed(ratio=1)
ggsave("redundancy_ABCDom_T0_Treatment.jpg", path = dirFigs, width = 6.75, height = 5, units = "in", dpi = 320)


