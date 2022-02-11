'DOC_analysis.R'

#library(ggplot2) are already in the package list that you should activate at the beginning
#library(cowplot)

# only uncomment if not in the work environement.
# DOC_dat <- read.csv(file.path(dirOutput, "DOC_dat.csv"))
# sym_dat_aquaria_final <- read.csv(file.path(dirOutput, "sym_dat_aquaria_final.csv"))

# Mege sym_dat_aquaria_final and DOC_dat.
DOC_dat1 <- merge(DOC_dat, sym_dat_aquaria_final, by.x="PLANC_aquaria", by.y="Group.1", all.x=T, all.y=F)

# Normalize DOC to surface area.
DOC_dat1$DOC_SA_Normalized = DOC_dat1$uMC/DOC_dat1$SA_no_outliers #normalize raw DOC
DOC_dat1$Control_Corrected_DOC_SA_Normalized = DOC_dat1$Control_Corrected_DOC/DOC_dat1$SA_no_outliers #normalize control corrected DOC


# Reorder treatment factors.
DOC_dat1$Treatment.x = factor(DOC_dat1$Treatment.x, levels=c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", "Bleached + Heated", "Ambient Water Control", "Heated Water Control"))

# Visualize raw DOC data with water controls.
#first subset for just t0 samples.
DOC_dat_t0 <- subset(DOC_dat1, Timepoint==0)
DOC_dat_t0$Treatment = DOC_dat_t0$Treatment.x

#create manual
cost.col.fill<-c("dodgerblue3","firebrick3", "white", "white", "grey70", "grey70", "grey70")
cost.col.line<-c("dodgerblue1", "firebrick1", "dodgerblue1", "firebrick1", "dodgerblue1", "firebrick1", "lightgrey")
fact.all.treat<-factor(DOC_dat_t0$Treatment, levels = c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", "Bleached + Heated", "Ambient Water Control", "Heated Water Control", "not applicable"))
cost.shape <- c(21,24)

#visualize
ggplot(DOC_dat_t0,aes(x=Treatment,y=uMC,color=Treatment,fill=Treatment))+
  geom_boxplot(size = 1.2)+
  geom_point(size = 3)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  # theme(legend.key.height=unit(0.5,"in"))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylab("DOC (uM)")+
  xlab("")
#use this figure for a potential supplement
ggsave('DOC_per treatment.jpeg', path = dirFigs, width = 6.75, height = 5, units = "in", dpi = 320)

# old plot
# ggplot(DOC_dat_t0,aes(x=Treatment,y=uMC,color=Treatment,fill=Treatment))+
#   geom_boxplot(size=2)+
#   geom_point(size=4)+
#   scale_color_manual(values=cost.col.line)+
#   scale_fill_manual(values=cost.col.fill )+
#   theme(axis.text=element_text(size=12),axis.title=element_text(size=12), axis.text.x=element_text(size=12), legend.text=element_text(size=12), legend.key.height=unit(0.8,"in"))+
#   scale_x_discrete(guide = guide_axis(n.dodge = 2))+
#   labs(y="DOC (uM)")
# #use this figure for a potential supplement
# ggsave('DOC_per treatment.jpeg', path = dirFigs, width = 8, height = 5, units = "in", dpi = 320)


#Check the distribution of raw DOC values.
hist(DOC_dat_t0$uMC) #not very normal
hist(log10(DOC_dat_t0$uMC)) #still not normal

# Since raw DOC values don't seem to be normal, test the effect of Treatment using a K-W test.
kruskal.test(DOC_dat_t0$uMC, DOC_dat_t0$Treatment) #marginally significant

#double check with an ANOVA
mod.DOC.t0 = aov(uMC ~ Treatment, data=DOC_dat_t0)
summary(mod.DOC.t0) #not significant


# Next, plot raw control corrected DOC values.
ggplot(DOC_dat_t0,aes(x=Treatment,y=Control_Corrected_DOC,color=Treatment,fill=Treatment))+
  geom_boxplot(size = 1.2)+
  geom_point(size = 3)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  # theme(legend.key.height=unit(0.5,"in"))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylab("Control Corrected DOC (uM)")+
  xlab("")
ggsave('DOC_control corrected_per treatment.jpeg', path = dirFigs, width = 6.75, height = 5, dpi = 300)

#Next plot coral specifc DOC values, first SA normalized DOC and then SA normalized control corrected DOC.
ggplot(DOC_dat_t0[DOC_dat_t0$Origin_PlanC.x != "control",],
       aes(x=Treatment,y=DOC_SA_Normalized,color=Treatment,fill=Treatment))+
  geom_boxplot(size = 1.2)+
  geom_point(size = 3)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  # theme(legend.key.height=unit(0.5,"in"))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylab ("Surface Area Normalized DOC (uM cm-2)")+
  xlab ("")
ggsave('DOC_Surface area normalized_per treatment.jpeg', path = dirFigs, width = 6.75, height = 5, dpi = 300)

ggplot(DOC_dat_t0[DOC_dat_t0$Origin_PlanC.x != "control",],aes(x=Treatment,y=Control_Corrected_DOC_SA_Normalized,color=Treatment,fill=Treatment))+
  geom_boxplot(size = 1.2)+
  geom_point(size = 3)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  # theme(legend.key.height=unit(0.5,"in"))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylab ("Surface Area Normalized DOC (uM cm-2)")+
  xlab("")
ggsave('DOC_Surface area normalized_control corrected_per treatment.jpeg', path = dirFigs, width = 6.75, height = 5, dpi = 300)
#Use as panel A for Figure 2

#Next, run stats on the control corrected SA normalized DOC data.

#first check the distribution of the data.
hist(DOC_dat_t0$Control_Corrected_DOC_SA_Normalized) #kind of normal looking

#run a K-W test.
kruskal.test(DOC_dat_t0$Control_Corrected_DOC_SA_Normalized, DOC_dat_t0$Treatment) #not significant.

#run an ANOVA
mod.DOC.t0.corrected.SA.normalized = aov(Control_Corrected_DOC_SA_Normalized ~ Treatment, data=DOC_dat_t0)
summary(mod.DOC.t0.corrected.SA.normalized) #still not significant.

