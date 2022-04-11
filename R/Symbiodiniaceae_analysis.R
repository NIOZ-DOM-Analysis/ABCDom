'Symbiodiniaceae_analysis.R

Goal: Analyze and visualize the differences in nubbin symbiodiniaceae between treatments and timepoints.
'

#libraries used, should be loaded already
# library(ggplot2)
# library(cowplot)
# library(FSA)

sym_dat <- read.csv(file.path(dirOutput, "sym_dat.csv"), ) #load the processed symbiodiniaceae data
sym_dat_aquaria_final <- read.csv(file.path(dirOutput, "sym_dat_aquaria_final.csv")) #load in the processed symbiodiniaceae data by aquaria
# metadata <- read.csv(file.path(dirRAW, "Metadata.csv"), ) #load in metadata



# Visualize the t0 corals sym densities (SA normalized) excluding stringent outliers (outliers) and relaxed outliers (outliers1).
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & !is.na(sym_dat$Collection_Bleaching_Level),],(aes(x=Collection_Bleaching_Level,y=log10.sym.SA,color=Collection_Bleaching_Level)))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot()+
  theme_classic()+
  ggtitle("Stringent Outliers")
ggsave('Symbiont density_log10_per collection bleaching level_stringent outliers.jpeg', path = dirFigs, dpi = 300)

ggplot(sym_dat[sym_dat$Outlier1!="Y" & sym_dat$Timepoint=="T0" & !is.na(sym_dat$Collection_Bleaching_Level),],(aes(x=Collection_Bleaching_Level,y=log10.sym.SA,color=Collection_Bleaching_Level)))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot()+
  theme_classic()+
  ggtitle("Relaxed Outliers")
ggsave('Symbiont density_log10_per collection bleaching level_relaxed outliers.jpeg', path = dirFigs, dpi = 300)


# Sym densities at T0 behave as you would expect.

# Visualize with BLEACHED and PARTIALLY BLEACHED lumped as BLEACHED treatment.
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & !is.na(sym_dat$Collection_Bleaching_Level),],(aes(x=Collection_Bleaching_Level1,y=log10.sym.SA,color=Collection_Bleaching_Level1)))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size=1.2)+
  theme_classic()+
  ggtitle("Stringent Outliers")
ggsave('Symbiont density_log10_Bleached_Healthy_Stringent outliers.jpeg', path = dirFigs, dpi = 300)


ggplot(sym_dat[sym_dat$Outlier1!="Y" & sym_dat$Timepoint=="T0" & !is.na(sym_dat$Collection_Bleaching_Level),],(aes(x=Collection_Bleaching_Level1,y=log10.sym.SA,color=Collection_Bleaching_Level1)))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size=1.2)+
  theme_classic()+
  ggtitle("Relaxed Outliers")
ggsave('Symbiont density_log10_Bleached_Healthy_Relaxed outliers.jpeg', path = dirFigs, dpi = 300)


ggplot(sym_dat[sym_dat$Outlier1!="Y" & sym_dat$Timepoint=="T0" & !is.na(sym_dat$Collection_Bleaching_Level),],(aes(x=Collection_Bleaching_Level1,y=sym.SA,color=Collection_Bleaching_Level1)))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size=1.2)+
  theme_classic()+
  ggtitle("Relaxed Outliers and not log transformed")
ggsave('Symbiont density_Bleached_Healthy_Relaxed outliers.jpeg', path = dirFigs, dpi = 300)


ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & !is.na(sym_dat$Collection_Bleaching_Level),],(aes(x=Species,y=sym.SA,color=Collection_Bleaching_Level1)))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size=1.2)+
  theme_classic()+
  ggtitle("Stringent Outliers and not log transformed")

ggsave('Symbiont density_Bleached_Healthy_Stringent outliers.jpeg', path = dirFigs, dpi = 300)

sym_dat$Collection_Bleaching_Level1 <- factor(sym_dat$Collection_Bleaching_Level1, levels = c("HE", "BL"))
#visualize nubbin sym densities at t0, split out by species.
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & is.na(sym_dat$Species)!=TRUE,],(aes(x=Collection_Bleaching_Level1,y=log10(sym.SA),color=Collection_Bleaching_Level1, fill=Collection_Bleaching_Level1)))+
  facet_wrap(.~Species, scales="fixed")+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size=1.2)+
  scale_color_manual(labels = c("Non-Bleached", "Bleached"), values=c("dodgerblue1", "dodgerblue1"))+
  scale_fill_manual(labels = c("Non-Bleached", "Bleached"), values=c("dodgerblue3", "white"))+
  scale_x_discrete(labels=NULL)+
  theme_classic()+
  theme(text=element_text(size=15),legend.key.height=unit(1.75,"cm"))+
  labs(y="Log10 Symbiodiniaceae cells per cm^2",x="Bleaching Status at Collection",color="Bleaching Status at Collection",fill="Bleaching Status at Collection")
##save and use as panel A for Fig 1. Add posthoc values manually
ggsave('Symbiont cells per cm2_Bleaching status at collection_v1.jpeg', path = dirFigs, dpi = 300, width=13, height=8)

#make the same plot but with free scales
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & is.na(sym_dat$Species)!=TRUE,],(aes(x=Collection_Bleaching_Level1,y=log10(sym.SA),color=Collection_Bleaching_Level1, fill=Collection_Bleaching_Level1)))+
  facet_wrap(.~Species, scales="free")+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size=1.2)+
  scale_color_manual(values=c("dodgerblue1", "dodgerblue1"))+
  scale_fill_manual(values=c("white","dodgerblue3"))+
  theme_classic()+
  theme(text=element_text(size=15),legend.key.height=unit(1.75,"cm"))+
  labs(y="Log10 Symbiodiniaceae cells per cm^2",x="Bleaching Status at Collection",color="Bleaching Status at Collection",fill="Bleaching Status at Collection")
##save and use as panel A for Fig 1. Add posthoc values manually
ggsave('Symbiont cells per cm2_Bleaching status at collection_v2.jpeg', path = dirFigs, dpi = 300, width=13, height=8)


# Next, visualize the t0 aquaria aggregated sym densities.

ggplot(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T0",],aes(x=Bleaching_Status_at_Collection,y=sym.SA.normalized.no.outliers,color=Bleaching_Status_at_Collection,fill=Bleaching_Status_at_Collection))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size=2)+
  # geom_point(size=3)+
  scale_color_manual(values=c("dodgerblue1", "dodgerblue1"))+
  scale_fill_manual(values=c("white","dodgerblue3"))+
  theme_classic()+
  theme(text=element_text(size=15),legend.key.height=unit(1.75,"cm"))+
  labs(y="Symbiodiniaceae cells per cm^2",x="Bleaching Status at Collection",color="Bleaching Status at Collection",fill="Bleaching Status at Collection")
ggsave('Symbiont cells per cm2_Bleaching status at collection.jpeg', path = dirFigs, dpi = 300)


#try with log10 transformed data
ggplot(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T0",],aes(x=Bleaching_Status_at_Collection,y=log10(sym.SA.normalized.no.outliers1),color=Bleaching_Status_at_Collection,fill=Bleaching_Status_at_Collection))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size=2)+
  # geom_point(size=3)+
  scale_color_manual(values=c("dodgerblue1","dodgerblue1"))+
  scale_fill_manual(values=c("white", "dodgerblue3"))+
  theme_classic()+
  theme(text=element_text(size=15),legend.key.height=unit(1.75,"cm"))+
  labs(y="log10 Symbiodiniaceae cells per cm^2",x="Bleaching Status at Collection",color="Bleaching Status at Collection",fill="Bleaching Status at Collection")
ggsave('Symbiont cells per cm2_logtrans_Bleaching status at collection.jpeg', path = dirFigs, dpi = 300)


# Now, visualize the distribution and run the statistics on the t0 aquaria data.

hist(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T0",]$sym.SA.normalized.no.outliers) #not normal
hist(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T0",]$sym.SA.normalized.no.outliers1) #try with the other outlier criteria. still not normal
hist(log10(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T0",]$sym.SA.normalized.no.outliers)) #not normal
hist(log10(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T0",]$sym.SA.normalized.no.outliers1)) #not normal

# Since the data are not normal, test the effect of Bleaching_Status_at_Collection using a k-w test.

kruskal.test(sym.SA.normalized.no.outliers ~ Bleaching_Status_at_Collection, sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T0",]) #significant


# Check with ANOVA,

mod.t0.sym=aov(sym.SA.normalized.no.outliers ~ Bleaching_Status_at_Collection, data=sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T0",])
summary(mod.t0.sym) #significant.

#try with log 10 data
mod.t0.sym.log=aov(log10(sym.SA.normalized.no.outliers) ~ Bleaching_Status_at_Collection, data=sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T0",])
summary(mod.t0.sym.log) #significant


#Next, plot the sym densities (SA normalized) for the T7 aquaria.

#reorder levels
sym_dat_aquaria_final$Treatment_v1<-factor(sym_dat_aquaria_final$Treatment_v1, levels = c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", "Bleached + Heated"))


ggplot(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",],aes(x=Treatment_v1,y=sym.SA.normalized.no.outliers,color=Treatment_v1,fill=Treatment_v1))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size=1.2)+
  # geom_point(size=3)+
  scale_color_manual(values=c("dodgerblue1","firebrick1","dodgerblue1","firebrick1"))+
  scale_fill_manual(values=c("dodgerblue3","firebrick3","white","white"))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), name = "")+
  theme_classic()+
  theme(text=element_text(size=15),legend.key.height=unit(1.75,"cm"))+
  labs(y="Symbiodiniaceae cells per cm^2",x="Treatment",color="Treatment",fill="Treatment")
#save and use as panel B for Figure 1.
ggsave('Symbiont cells per cm2_per treatment.jpeg', path = dirFigs, dpi = 300)



# Now, visualize the distribution and run the statistics on the t7 aquaria data.

hist(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",]$sym.SA.normalized.no.outliers) #normal-ish
hist(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",]$sym.SA.normalized.no.outliers1) #a little less normal looking
hist(log10(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",]$sym.SA.normalized.no.outliers)) #normal ish
hist(log10(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",]$sym.SA.normalized.no.outliers1)) #more normal.

# Since the data look normal enough, test the effect of Treatment_v1 using an ANOVA.

mod.t7.sym=aov(sym.SA.normalized.no.outliers ~ Treatment_v1, data=sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",])
summary(mod.t7.sym) #significant

#run tukey posthoc
TukeyHSD(mod.t7.sym)


