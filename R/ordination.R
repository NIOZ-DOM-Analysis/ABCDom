#ordination of data
library(BiodiversityR)
library(openxlsx)

# first make everything for df.area, later we can repeat everything for the normalized and logtransformed data.
# frist all, then we split the dataset in two

#non metric multidimentional scaling, all data
ord.mod.area<- metaMDS(df.area[M:ncol(df.area)], distance = 'bray', k = 2)
ordi.plot1.1<-ordiplot(ord.mod.area, choices = c(1,2))
sites.long1<-sites.long(ordi.plot1.1, env.data = df.area[1:(M-1)])

#create object/list to store the sites.lists in.
ordi.data <- list()
ordi.data$NMDS.all.data <- sites.long1

ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data=sites.long1,
             aes(x=axis1, y=axis2, colour=Treatment, shape=Experiment),
             size=5) +
  coord_fixed(ratio=1)
ggsave("NMDS_all_data.jpg", path = dirFigs, width = 8.5, height = 6.5, units = "in", dpi = 320)

#non metric multidimentional scaling, PlanC
df.area.planC <- df.area %>% filter(Experiment == "PLANC")
ord.mod.area<- metaMDS(df.area.planC[M:ncol(df.area.planC)], distance = 'bray', k = 2)
ordi.plot1.1<-ordiplot(ord.mod.area, choices = c(1,2))
sites.long1<-sites.long(ordi.plot1.1, env.data = df.area.planC[1:(M-1)])
ordi.data$NMDS.PlanC <- sites.long1

ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data=sites.long(ordi.plot1.1, env.data=df.area.planC[1:(M-1)]),
             aes(x=axis1, y=axis2, colour=Treatment, shape=Experiment),
             size=5) +
  coord_fixed(ratio=1)
ggsave("NMDS_PlanC.jpg", path = dirFigs, width = 8.5, height = 6.5, units = "in", dpi = 320)


#non metric multidimentional scaling, ABCDom
df.area.ABC <- df.area %>% filter(Experiment == "ABCDOM")
ord.mod.area<- metaMDS(df.area.ABC[M:ncol(df.area.ABC)], distance = 'bray', k = 2)
ordi.plot1.1<-ordiplot(ord.mod.area, choices = c(1,2))
sites.long1<-sites.long(ordi.plot1.1, env.data = df.area.ABC[1:(M-1)])
ordi.data$NMDS.ABCDom <- sites.long1

ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data=sites.long(ordi.plot1.1, env.data=df.area.ABC[1:(M-1)]),
             aes(x=axis1, y=axis2, colour=Treatment, shape=Timepoint_char),
             size=5) +
  coord_fixed(ratio=1)
ggsave("NMDS_ABCDom.jpg", path = dirFigs, width = 8.5, height = 6.5, units = "in", dpi = 320)

#non metric multidimentional scaling, ABCDom T0 and PLanC
df.area.T0 <- df.area %>% filter(Timepoint == 0)
ord.mod.area<- metaMDS(df.area.T0[M:ncol(df.area.T0)], distance = 'bray', k = 2)
ordi.plot1.1<-ordiplot(ord.mod.area, choices = c(1,2))
sites.long1<-sites.long(ordi.plot1.1, env.data = df.area.T0[1:(M-1)])
ordi.data$NMDS.ABCDom.PlanC.T0 <- sites.long1

ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data=sites.long(ordi.plot1.1, env.data=df.area.T0[1:(M-1)]),
             aes(x=axis1, y=axis2, colour=Treatment, shape=Experiment),
             size=5) +
  coord_fixed(ratio=1)
ggsave("NMDS_T0_PlanC_ABCDom.jpg", path = dirFigs, width = 8.5, height = 6.5, units = "in", dpi = 320)

#non metric multidimentional scaling, ABCT0Dom T0
df.area.ABCT0 <- df.area %>% filter(Experiment == "ABCDOM", Timepoint_char == "T0")
ord.mod.area<- metaMDS(df.area.ABCT0[M:ncol(df.area.ABCT0)], distance = 'bray', k = 2)
ordi.plot1.1<-ordiplot(ord.mod.area, choices = c(1,2))
sites.long1<-sites.long(ordi.plot1.1, env.data = df.area.ABCT0[1:(M-1)])
ordi.data$NMDS.ABCDom.T0 <- sites.long1

ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data=sites.long(ordi.plot1.1, env.data=df.area.ABCT0[1:(M-1)]),
             aes(x=axis1, y=axis2, colour=Treatment, shape=Timepoint_char),
             size=5) +
  coord_fixed(ratio=1)
ggsave("NMDS_ABCDom_T0.jpg", path = dirFigs, width = 8.5, height = 6.5, units = "in", dpi = 320)


############################
# redundancy analyis all data

all.data.hellinger <- disttransform(df.area[M:ncol(df.area)], method = 'hellinger')
ord.mod.area2<-rda(all.data.hellinger ~ Treatment, data = df.area[1:(M-1)], scaling = "species")
summary(ord.mod.area2)
ordi.plot2<-ordiplot(ord.mod.area2, choices = c(1,2))
sites.long2<- sites.long(ordi.plot2, env.data = df.area[1:(M-1)])
species.long2<-species.long(ordi.plot2)
axis.long2<-axis.long(ord.mod.area2, choices = c(1,2))


# ggplot() +
#   geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
#   geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
#   xlab(axis.long2[1, "label"]) +
#   ylab(axis.long2[2, "label"]) +
#   scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
#   scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
#   geom_point(data=sites.long2,
#              aes(x=axis1, y=axis2, colour=Treatment, shape=Experiment),
#              size=5) +
#   geom_point(data=species.long2,
#              aes(x=axis1, y=axis2)) +
#   coord_fixed(ratio=1)

spec.envfit<-envfit(ordi.plot2, env = all.data.hellinger, silent = TRUE)
spec.data.envit <- data.frame(r=spec.envfit$vector$r, p=spec.envfit$vectors$pvals)
species.long2 <- species.long(ordi.plot2, spec.data=spec.data.envit)
species.long3 <- species.long2[species.long2$r >= 0.8, ]
ordi.data$Redun.all.data <- sites.long2
ordi.data$drivers.redund.all.data <- species.long2

ggplot()+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2)+
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2)+
  xlab(axis.long2[1, "label"])+
  ylab(axis.long2[2, "label"])+
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL))+
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL))+
  geom_point(data=sites.long2, aes(x=axis1, y=axis2, colour=Treatment, shape=Timepoint_char, size=5))+
  geom_segment(data=species.long3, aes(x=0, y=0, xend=axis1*4, yend=axis2*4), colour="red", size=0.7, arrow=arrow())+
  #geom_text(data=species.long3, aes(x=axis1*4, y=axis2*4, label=labels),colour="red")+
  coord_fixed(ratio=1)
ggsave("redundancy_all_data_Treatment_Timepoint.jpg", path = dirFigs, width = 8.5, height = 6.5, units = "in", dpi = 320)
############################
# redundancy analyis ABCDom

ABCD.hellinger <- disttransform(df.area.ABC[M:ncol(df.area.ABC)], method = 'hellinger')
ord.mod.area2<-rda(ABCD.hellinger ~ Treatment, data = df.area.ABC[1:(M-1)], scaling = "species")
summary(ord.mod.area2)
ordi.plot2<-ordiplot(ord.mod.area2, choices = c(1,2))
sites.long2<- sites.long(ordi.plot2, env.data = df.area.ABC[1:(M-1)])
species.long2<-species.long(ordi.plot2)
axis.long2<-axis.long(ord.mod.area2, choices = c(1,2))
# ggplot() +
#   geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
#   geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
#   xlab(axis.long2[1, "label"]) +
#   ylab(axis.long2[2, "label"]) +
#   scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
#   scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
#   geom_point(data=sites.long2,
#              aes(x=axis1, y=axis2, colour=Treatment, shape=Timepoint_char),
#              size=5) +
#   geom_point(data=species.long2,
#              aes(x=axis1, y=axis2)) +
#   coord_fixed(ratio=1)

spec.envfit<-envfit(ordi.plot2, env = ABCD.hellinger, silent = TRUE)
spec.data.envit <- data.frame(r=spec.envfit$vector$r, p=spec.envfit$vectors$pvals)
species.long2 <- species.long(ordi.plot2, spec.data=spec.data.envit)
species.long3 <- species.long2[species.long2$r >= 0.9, ]
ordi.data$Redun.ABCDom <- sites.long2
ordi.data$drivers.redund.ABCDom <- species.long2

ggplot()+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2)+
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2)+
  xlab(axis.long2[1, "label"])+
  ylab(axis.long2[2, "label"])+
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL))+
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL))+
  geom_point(data=sites.long2, aes(x=axis1, y=axis2, colour=Treatment, shape=Timepoint_char, size=5))+
  geom_segment(data=species.long3, aes(x=0, y=0, xend=axis1*4, yend=axis2*4), colour="red", size=0.7, arrow=arrow())+
  #geom_text(data=species.long3, aes(x=axis1*4, y=axis2*4, label=labels),colour="red")+
  coord_fixed(ratio=1)
ggsave("redundancy_ABCDom_Treatment_Timepoint.jpg", path = dirFigs, width = 8.5, height = 6.5, units = "in", dpi = 320)

############################
# redundancy analyis PlanC

PlanC.hellinger <- disttransform(df.area.planC[M:ncol(df.area.planC)], method = 'hellinger')
ord.mod.area2<-rda(PlanC.hellinger ~ Treatment, data = df.area.planC[1:(M-1)], scaling = "species")
summary(ord.mod.area2)
ordi.plot2<-ordiplot(ord.mod.area2, choices = c(1,2))
sites.long2<- sites.long(ordi.plot2, env.data = df.area.planC[1:(M-1)])
species.long2<-species.long(ordi.plot2)
axis.long2<-axis.long(ord.mod.area2, choices = c(1,2))
# ggplot() +
#   geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
#   geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
#   xlab(axis.long2[1, "label"]) +
#   ylab(axis.long2[2, "label"]) +
#   scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
#   scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
#   geom_point(data=sites.long2,
#              aes(x=axis1, y=axis2, colour=Treatment, shape=Timepoint_char),
#              size=5) +
#   geom_point(data=species.long2,
#              aes(x=axis1, y=axis2)) +
#   coord_fixed(ratio=1)

spec.envfit<-envfit(ordi.plot2, env = PlanC.hellinger, silent = TRUE)
spec.data.envit <- data.frame(r=spec.envfit$vector$r, p=spec.envfit$vectors$pvals)
species.long2 <- species.long(ordi.plot2, spec.data=spec.data.envit)
species.long3 <- species.long2[species.long2$r >= 0.8, ]
ordi.data$Redun.PlanC <- sites.long2
ordi.data$drivers.redund.PlanC <- species.long2

ggplot()+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2)+
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2)+
  xlab(axis.long2[1, "label"])+
  ylab(axis.long2[2, "label"])+
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL))+
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL))+
  geom_point(data=sites.long2, aes(x=axis1, y=axis2, colour=Treatment, shape=Timepoint_char, size=5))+
  geom_segment(data=species.long3, aes(x=0, y=0, xend=axis1*4, yend=axis2*4), colour="red", size=0.7, arrow=arrow())+
  #geom_text(data=species.long3, aes(x=axis1*4, y=axis2*4, label=labels),colour="red")+
  coord_fixed(ratio=1)
ggsave("redundancy_PlanC_Treatment_Timepoint.jpg", path = dirFigs, width = 8.5, height = 6.5, units = "in", dpi = 320)

############################
# redundancy analyis T0

T0.hellinger <- disttransform(df.area.T0[M:ncol(df.area.T0)], method = 'hellinger')
ord.mod.area2<-rda(T0.hellinger ~ Treatment, data = df.area.T0[1:(M-1)], scaling = "species")
summary(ord.mod.area2)
ordi.plot2<-ordiplot(ord.mod.area2, choices = c(1,2))
sites.long2<- sites.long(ordi.plot2, env.data = df.area.T0[1:(M-1)])
species.long2<-species.long(ordi.plot2)
axis.long2<-axis.long(ord.mod.area2, choices = c(1,2))
# ggplot() +
#   geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
#   geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
#   xlab(axis.long2[1, "label"]) +
#   ylab(axis.long2[2, "label"]) +
#   scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
#   scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
#   geom_point(data=sites.long2,
#              aes(x=axis1, y=axis2, colour=Treatment, shape=Timepoint_char),
#              size=5) +
#   geom_point(data=species.long2,
#              aes(x=axis1, y=axis2)) +
#   coord_fixed(ratio=1)

spec.envfit<-envfit(ordi.plot2, env = T0.hellinger, silent = TRUE)
spec.data.envit <- data.frame(r=spec.envfit$vector$r, p=spec.envfit$vectors$pvals)
species.long2 <- species.long(ordi.plot2, spec.data=spec.data.envit)
species.long3 <- species.long2[species.long2$r >= 0.8, ]
ordi.data$Redun.ABCDom.PlanC.T0 <- sites.long2
ordi.data$drivers.redund.ABCDom.PlanC.T0 <- species.long2

ggplot()+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2)+
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2)+
  xlab(axis.long2[1, "label"])+
  ylab(axis.long2[2, "label"])+
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL))+
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL))+
  geom_point(data=sites.long2, aes(x=axis1, y=axis2, colour=Treatment, shape=Experiment, size=5))+
  geom_segment(data=species.long3, aes(x=0, y=0, xend=axis1*4, yend=axis2*4), colour="red", size=0.7, arrow=arrow())+
  #geom_text(data=species.long3, aes(x=axis1*4, y=axis2*4, label=labels),colour="red")+
  coord_fixed(ratio=1)
ggsave("redundancy_ABCDom_PlanC_T0_Treatment_Experiment.jpg", path = dirFigs, width = 8.5, height = 6.5, units = "in", dpi = 320)

########
# redundancy analyis ABCDomT0
ABCDT0.hellinger <- disttransform(df.area.ABCT0[M:ncol(df.area.ABCT0)], method = 'hellinger')
ord.mod.area2<-rda(ABCDT0.hellinger ~ Treatment, data = df.area.ABCT0[1:(M-1)], scaling = "species")
summary(ord.mod.area2)
ordi.plot2<-ordiplot(ord.mod.area2, choices = c(1,2))
sites.long2<- sites.long(ordi.plot2, env.data = df.area.ABCT0[1:(M-1)])
species.long2<-species.long(ordi.plot2)
axis.long2<-axis.long(ord.mod.area2, choices = c(1,2))
# ggplot() +
#   geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
#   geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
#   xlab(axis.long2[1, "label"]) +
#   ylab(axis.long2[2, "label"]) +
#   scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
#   scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
#   geom_point(data=sites.long2,
#              aes(x=axis1, y=axis2, colour=Treatment, shape=Timepoint_char),
#              size=5) +
#   geom_point(data=species.long2,
#              aes(x=axis1, y=axis2)) +
#   coord_fixed(ratio=1)

spec.envfit<-envfit(ordi.plot2, env = ABCDT0.hellinger, silent = TRUE)
spec.data.envit <- data.frame(r=spec.envfit$vector$r, p=spec.envfit$vectors$pvals)
species.long2 <- species.long(ordi.plot2, spec.data=spec.data.envit)
species.long3 <- species.long2[species.long2$r >= 0.7, ]
ordi.data$Redun.ABCDom.T0 <- sites.long2
ordi.data$drivers.redund.ABCDom.T0 <- species.long2

ggplot()+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2)+
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2)+
  xlab(axis.long2[1, "label"])+
  ylab(axis.long2[2, "label"])+
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL))+
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL))+
  geom_point(data=sites.long2, aes(x=axis1, y=axis2, colour=Treatment, size=5))+
  geom_segment(data=species.long3, aes(x=0, y=0, xend=axis1*4, yend=axis2*4), colour="red", size=0.7, arrow=arrow())+
  #geom_text(data=species.long3, aes(x=axis1*4, y=axis2*4, label=labels),colour="red")+
  coord_fixed(ratio=1)
ggsave("redundancy_ABCDom_T0_Treatment.jpg", path = dirFigs, width = 8.5, height = 6.5, units = "in", dpi = 320)


write.xlsx(ordi.data, paste0(dirOutput,"/ordination_data.xlsx"), overwrite = TRUE)

rm(ordi.plot1.1, ordi.plot2, species.long2, species.long3, spec.data.envit, spec.envfit, axis.long2, sites.long1, sites.long2, ord.mod.area, ord.mod.area2, ABCD.hellinger, ABCDT0.hellinger, T0.hellinger, PlanC.hellinger, all.data.hellinger, df.area.ABC, df.area.ABCT0, df.area.planC, df.area.T0)
