'DOC_analysis.R'

library(ggplot2) #are already in the package list that you should activate at the beginning
library(cowplot)

# only uncomment if not in the work environement.
# DOC_dat <- read.csv(file.path(dirOutput, "DOC_dat.csv"))
# sym_dat_aquaria_final <- read.csv(file.path(dirOutput, "sym_dat_aquaria_final.csv"))

# Mege sym_dat_aquaria_final and DOC_dat.
DOC_dat1 <- merge(DOC_dat, sym_dat_aquaria_final, by.x="PLANC_aquaria", by.y="Group.1", all.x=T, all.y=F)

# Normalize DOC to surface area.
DOC_dat1$DOC_SA_Normalized = DOC_dat1$uMC/DOC_dat1$SA_no_outliers #normalize raw DOC
DOC_dat1$Control_Corrected_DOC_SA_Normalized = DOC_dat1$Control_Corrected_DOC/DOC_dat1$SA_no_outliers #normalize control corrected DOC
DOC_dat1$Control_Corrected_DOC_flux_SA_Normalized <- DOC_dat1$Control_Corrected_DOC_flux/DOC_dat1$SA_no_outliers #normalize control corrected DOC fluxes

# Reorder treatment factors.
DOC_dat1$Treatment.x = factor(DOC_dat1$Treatment.x, levels=c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", "Bleached + Heated", "Ambient Water Control", "Heated Water Control"))

# Visualize raw DOC data with water controls.
#first subset for just t0 samples.
DOC_dat_t0 <- subset(DOC_dat1, Timepoint==0)
DOC_dat_t0$Treatment = DOC_dat_t0$Treatment.x

#create manual
cost.col.fill<-c("dodgerblue1","firebrick1", "white", "white", "grey70", "grey70", "grey70")
cost.col.line<-c("dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "lightgrey")
fact.all.treat<-factor(DOC_dat_t0$Treatment, levels = c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", "Bleached + Heated", "Ambient Water Control", "Heated Water Control"))
cost.shape <- c(21,24)

#visualize
ggplot(DOC_dat_t0,aes(x=Treatment,y=uMC,color=Treatment,fill=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  # geom_point(size = 3)+
  scale_color_manual(values=cost.col.line)+
  # scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  scale_fill_manual(values=cost.col.fill)+
  theme_classic()+
  # theme(legend.key.height=unit(0.5,"in"))+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  theme(legend.position = "none")+
  ylab("DOC (uM)")+
  xlab("")
#use this figure for a potential supplement
ggsave('DOC_per treatment.jpeg', path = dirFigs, width = 9, height = 5.5, units = "in", dpi = 320)

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
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  # geom_point(size = 3)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill)+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  theme_classic()+
  theme(legend.position = "none")+
  # theme(legend.key.height=unit(0.5,"in"))+
  ylab("Control Corrected DOC (uM)")+
  xlab("")
ggsave('DOC_control corrected_per treatment.jpeg', path = dirFigs, width = 9, height = 5.5, dpi = 300)

#Next plot coral specifc DOC values, first SA normalized DOC and then SA normalized control corrected DOC.
ggplot(DOC_dat_t0[DOC_dat_t0$Origin_PlanC.x != "control",],
       aes(x=Treatment,y=DOC_SA_Normalized,color=Treatment,fill=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  # geom_point(size = 3)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill)+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  theme_classic()+
  theme(legend.position = "none")+
  # theme(legend.key.height=unit(0.5,"in"))+
  ylab ("Surface Area Normalized DOC (uM cm-2)")+
  xlab ("")
ggsave('DOC_Surface area normalized_per treatment.jpeg', path = dirFigs, width = 9, height = 5.5, dpi = 300)

ggplot(DOC_dat_t0[DOC_dat_t0$Origin_PlanC.x != "control",],aes(x=Treatment,y=Control_Corrected_DOC_SA_Normalized,color=Treatment,fill=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  # geom_point(size = 3)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill)+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  theme_classic()+
  theme(legend.position = "none")+
  ylab ("Surface Area Normalized DOC (uM cm-2)")+
  xlab("")
ggsave('DOC_Surface area normalized_control corrected_per treatment.jpeg', path = dirFigs, width = 9, height = 5.5, dpi = 300)

#first check the distribution of the data.
hist(DOC_dat_t0$Control_Corrected_DOC_SA_Normalized) #kind of normal looking

#run a K-W test.
kruskal.test(DOC_dat_t0$Control_Corrected_DOC_SA_Normalized, DOC_dat_t0$Treatment) #not significant.

#run an ANOVA
mod.DOC.t0.corrected.SA.normalized = aov(Control_Corrected_DOC_SA_Normalized ~ Treatment, data=DOC_dat_t0)
summary(mod.DOC.t0.corrected.SA.normalized) #still not significant.

#Next, plot the DOC control corrected SA normalized fluxes for coral samples
DOC_dat_t0$Control_Corrected_DOC_flux_SA_Normalized_v1 <- DOC_dat_t0$Control_Corrected_DOC_flux_SA_Normalized*100 #convert from cm-2 to dm-2

fig2A<-ggplot(DOC_dat_t0[DOC_dat_t0$Origin_PlanC.x != "control",],aes(x=Treatment,y=Control_Corrected_DOC_flux_SA_Normalized_v1,color=Treatment,fill=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  # geom_point(size = 3)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill)+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  theme_classic()+
  theme(legend.position = "none")+
  ylab ("Surface Area Normalized DOC flux (uM dm-2 h-1)")+
  xlab("")
ggsave('DOC_flux Surface area normalized_control corrected_per treatment.jpeg', path = dirFigs, width = 7.5, height = 5.5, dpi = 300)
#Use as panel A for Figure 2

#now run stats on DOC flux
#first check the distribution of the data.
hist(DOC_dat_t0$Control_Corrected_DOC_flux_SA_Normalized) #kind of normal looking

#run a K-W test.
kruskal.test(DOC_dat_t0$Control_Corrected_DOC_flux_SA_Normalized, DOC_dat_t0$Treatment) #not significant.

#run an ANOVA
mod.DOC.t0.corrected.SA.normalized.flux = aov(Control_Corrected_DOC_flux_SA_Normalized ~ Treatment, data=DOC_dat_t0)
summary(mod.DOC.t0.corrected.SA.normalized.flux) #still not significant.

#calculate DOC drawdown
DOC_drawdown <- subset(DOC_dat1, Timepoint==0) #subset for t0
DOC_drawdown1 <- DOC_drawdown[,c(1:20)] #subset columns
colnames(DOC_drawdown1)[7:8] <- c("DOC_t0", "stdev_t0")
DOC_drawdown2 <-merge(DOC_drawdown1, subset(DOC_dat1, Timepoint==36)[,c(7:8,20)], by.x="Bottle_NR.x", by.y="Bottle_NR.x") #merge t36 DOC data
colnames(DOC_drawdown2)[21:22] <- c("DOC_t36", "stdev_t36")
DOC_drawdown2$DOC_drawdown <- DOC_drawdown2$DOC_t0 - DOC_drawdown2$DOC_t36 #calculate drawdown
DOC_drawdown2$DOC_percent_drawdown <- 100*DOC_drawdown2$DOC_drawdown/DOC_drawdown2$DOC_t0

#visualize DOC drawdown
ggplot(DOC_drawdown2[DOC_drawdown2$Treatment.x!="Bleached + Heated",], aes(x=Treatment.x, y=DOC_drawdown, color=Treatment.x, fill=Treatment.x))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  # geom_point(size = 3)+
  scale_color_manual(values=cost.col.line[-4])+
  scale_fill_manual(values=cost.col.fill[-4])+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Negative Control", "Negative Control + Heated"))+
  theme_classic()+
  theme(legend.position = "none")+
  ylab ("DOC Drawdown (uM)")+
  xlab("")+
  labs(x="", color="Treatment", fill="Treatment")

#visualize percent DOC drawdown
ggplot(DOC_drawdown2[DOC_drawdown2$Treatment.x!="Bleached + Heated",], aes(x=Treatment.x, y=DOC_percent_drawdown, color=Treatment.x, fill=Treatment.x))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  # geom_point(size = 3)+
  scale_color_manual(values=cost.col.line[-4])+
  scale_fill_manual(values=cost.col.fill[-4])+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Negative Control", "Negative Control + Heated"))+
  theme_classic()+
  theme(legend.position = "none")+
  ylab ("DOC Percent Drawdown")+
  xlab("")
ggsave('DOC_percent_drawdown_per_treatment.jpeg', path = dirFigs, width = 11, height = 8, dpi = 300)

fig2B<-ggplot(DOC_drawdown2[], aes(x=Treatment.x, y=DOC_percent_drawdown, color=Treatment.x, fill=Treatment.x))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  # geom_point(size = 3)+
  scale_color_manual(values=c("dodgerblue3", "firebrick3", "dodgerblue3", "white", "dodgerblue3", "firebrick3", "lightgrey"))+
  scale_fill_manual(values=cost.col.fill[])+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  theme_classic()+
  theme(legend.position = "none")+
  ylab ("DOC Percent Drawdown")+
  xlab("")
fig2B


#run stats on percent drawdown
#first visualize distribution
hist(DOC_drawdown2$DOC_percent_drawdown) #normal ish

#run model
mod.doc.percent.drawdown <- aov(DOC_percent_drawdown ~ Treatment.x, data=DOC_drawdown2)
summary(mod.doc.percent.drawdown) #not significant

mod.doc.percent.drawdown1 <- aov(DOC_percent_drawdown ~ Treatment.x, data=DOC_drawdown2[DOC_drawdown2$Treatment.x!="Bleached + Heated",])
summary(mod.doc.percent.drawdown1) #not significant

kruskal.test(DOC_drawdown2$DOC_percent_drawdown, DOC_drawdown2$Treatment.x)

kruskal.test(DOC_drawdown2$DOC_percent_drawdown[DOC_drawdown2$Treatment.x!="Bleached + Heated"], DOC_drawdown2$Treatment.x[DOC_drawdown2$Treatment.x!="Bleached + Heated"])
