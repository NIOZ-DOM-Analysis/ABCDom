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
#ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & !is.na(sym_dat$Collection_Bleaching_Level),],(aes(x=Collection_Bleaching_Level,y=log10.sym.SA,color=Collection_Bleaching_Level)))+
#  stat_boxplot(geom = 'errorbar', size = 2)+
#  geom_boxplot()+
#  theme_classic()+
#  ggtitle("Stringent Outliers")
#ggsave('Symbiont density_log10_per collection bleaching level_stringent outliers.jpeg', path = dirFigs, dpi = 300)

#ggplot(sym_dat[sym_dat$Outlier1!="Y" & sym_dat$Timepoint=="T0" & !is.na(sym_dat$Collection_Bleaching_Level),],(aes(x=Collection_Bleaching_Level,y=log10.sym.SA,color=Collection_Bleaching_Level)))+
#  stat_boxplot(geom = 'errorbar', size = 2)+
#  geom_boxplot()+
#  theme_classic()+
#  ggtitle("Relaxed Outliers")
#ggsave('Symbiont density_log10_per collection bleaching level_relaxed outliers.jpeg', path = dirFigs, dpi = 300)


# Sym densities at T0 behave as you would expect.

# Visualize with BLEACHED and PARTIALLY BLEACHED lumped as BLEACHED treatment.
#ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & !is.na(sym_dat$Collection_Bleaching_Level),],(aes(x=Collection_Bleaching_Level1,y=log10.sym.SA,color=Collection_Bleaching_Level1)))+
#  stat_boxplot(geom = 'errorbar', size = 2)+
#  geom_boxplot(size=1.2)+
#  theme_classic()+
#  ggtitle("Stringent Outliers")
#ggsave('Symbiont density_log10_Bleached_Healthy_Stringent outliers.jpeg', path = dirFigs, dpi = 300)


#ggplot(sym_dat[sym_dat$Outlier1!="Y" & sym_dat$Timepoint=="T0" & !is.na(sym_dat$Collection_Bleaching_Level),],(aes(x=Collection_Bleaching_Level1,y=log10.sym.SA,color=Collection_Bleaching_Level1)))+
#  stat_boxplot(geom = 'errorbar', size = 2)+
#  geom_boxplot(size=1.2)+
#  theme_classic()+
#  ggtitle("Relaxed Outliers")
#ggsave('Symbiont density_log10_Bleached_Healthy_Relaxed outliers.jpeg', path = dirFigs, dpi = 300)


#ggplot(sym_dat[sym_dat$Outlier1!="Y" & sym_dat$Timepoint=="T0" & !is.na(sym_dat$Collection_Bleaching_Level),],(aes(x=Collection_Bleaching_Level1,y=sym.SA,color=Collection_Bleaching_Level1)))+
#  stat_boxplot(geom = 'errorbar', size = 2)+
#  geom_boxplot(size=1.2)+
#  theme_classic()+
#  ggtitle("Relaxed Outliers and not log transformed")
#ggsave('Symbiont density_Bleached_Healthy_Relaxed outliers.jpeg', path = dirFigs, dpi = 300)


#ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & !is.na(sym_dat$Collection_Bleaching_Level),],(aes(x=Species,y=sym.SA,color=Collection_Bleaching_Level1)))+
#  stat_boxplot(geom = 'errorbar', size = 2)+
#  geom_boxplot(size=1.2)+
#  theme_classic()+
#  ggtitle("Stringent Outliers and not log transformed")

#ggsave('Symbiont density_Bleached_Healthy_Stringent outliers.jpeg', path = dirFigs, dpi = 300)

sym_dat$Collection_Bleaching_Level1 <- factor(sym_dat$Collection_Bleaching_Level1, levels = c("HE", "BL"))
#visualize nubbin sym densities at t0, split out by species

#first, generate a combined species_bleaching level for TukeyHSD purposes.
sym_dat$species_bleaching <- paste(sym_dat$Species, sym_dat$Collection_Bleaching_Level1, sep="_")

#calculate maximum for each species_bleaching group. CANT FIGURE OUT HOW TO DO THIS WITH FACETS.
#sym_dat_max <- aggregate(sym_dat[,40:41], by=list(sym_dat$species_bleaching), FUN=max, na.rm=T) #be sure to exclude NAs
#sym_dat_max1 <- sym_dat_max[-3,] #remove NA row
#sym_dat_max1$posthoc <- c("A","C","A/B","C","B/C","C")
#sym_dat_max1$Collection_Bleaching_Level1 <- c("BL", "HE", "BL", "HE", "BL", "HE")

#first with raw densities
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & is.na(sym_dat$Species)!=TRUE,],(aes(x=Collection_Bleaching_Level1,y=sym.SA,color=Collection_Bleaching_Level1, fill=Collection_Bleaching_Level1)))+
  facet_wrap(.~Species, scales="fixed")+
  stat_boxplot(geom = 'errorbar', size = 1.2)+
  geom_boxplot(size=1.2)+
  scale_color_manual(labels = c("Non-bleached", "Bleached"), values=c("dodgerblue3", "dodgerblue3"))+
  scale_fill_manual(labels = c("Non-bleached", "Bleached"), values=c("dodgerblue1", "white"))+
  scale_x_discrete(labels=NULL)+
  theme_classic()+
  # theme(text=element_text(size=24),legend.key.height=unit(2,"cm"))+
  labs(y="Symbiodiniaceae cells per cm^2",x="Bleaching Status at Collection",color="Bleaching Status at Collection",fill="Bleaching Status at Collection")

#Then with log10 densities. This one looks better than raw data and is probably more normal, so save this for figure 1
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & is.na(sym_dat$Species)!=TRUE,],(aes(x=Collection_Bleaching_Level1,y=log10(sym.SA),color=Collection_Bleaching_Level1, fill=Collection_Bleaching_Level1)))+
  facet_wrap(.~Species, scales="fixed")+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  geom_boxplot(size=2)+
  scale_color_manual(labels = c("Unbleached", "Bleached"), values=c("dodgerblue3", "dodgerblue3"))+
  scale_fill_manual(labels = c("Unbleached", "Bleached"), values=c("dodgerblue1", "white"))+
  scale_x_discrete(labels=NULL)+
  theme_classic()+
  theme(text=element_text(size=34),legend.key.height=unit(2,"cm"), plot.margin=unit(c(1,1,1,1.5), "cm"), strip.text.x=element_text(size=18))+
  ylab(expression(paste("Log10 Symbiodiniaceae\ncells per cm^2")))+
  xlab("")+
  coord_cartesian(ylim = c(3, 6))+
  labs(color="Bleaching Status\nat Collection",fill="Bleaching Status\nat Collection")
  ##save and use as panel A for Fig 1. Add posthoc values manually
ggsave('Symbiont cells per cm2_Bleaching status at collection_v3.jpeg', path = dirFigs, dpi = 300, width=15, height=9)

#make the same plot but with free scales
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & is.na(sym_dat$Species)!=TRUE,],(aes(x=Collection_Bleaching_Level1,y=log10(sym.SA),color=Collection_Bleaching_Level1, fill=Collection_Bleaching_Level1)))+
  facet_wrap(.~Species, scales="free")+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  geom_boxplot(size=2)+
  scale_color_manual(values=c("dodgerblue3", "dodgerblue3"))+
  scale_fill_manual(values=c("dodgerblue1","white"))+
  theme_classic()+
  theme(text=element_text(size=24),legend.key.height=unit(2,"cm"))+
  labs(y="Log10 Symbiodiniaceae cells per cm^2",x="Bleaching Status at Collection",color="Bleaching Status at Collection",fill="Bleaching Status at Collection")
ggsave('Symbiont cells per cm2_Bleaching status at collection_v2.jpeg', path = dirFigs, dpi = 300, width=15, height=9)

# Now, visualize the distribution and run the statistics on the t0 sym_dat
hist(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & is.na(sym_dat$Species)!=TRUE,]$sym.SA) #not normal
hist(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & is.na(sym_dat$Species)!=TRUE,]$log10.sym.SA) #normal

# Since the log10 data are normal, run linear model

#run model
mod.t0.sym=aov(log10.sym.SA ~ Collection_Bleaching_Level1*Species, data=sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & is.na(sym_dat$Species)!=TRUE,])
summary(mod.t0.sym) #bleaching level is significant, species is significant, and the interaction is significant.

#run model with species_bleaching
mod.t0.sym.species.bleaching=aov(log10.sym.SA ~ species_bleaching, data=sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T0" & is.na(sym_dat$Species)!=TRUE,])

#run tukeyHSD for visualization
TukeyHSD(mod.t0.sym.species.bleaching, "species_bleaching")

#Next, plot the sym densities (SA normalized) for the T7 aquaria.

#reorder levels
sym_dat_aquaria_final$Treatment_v1<-factor(sym_dat_aquaria_final$Treatment_v1, levels = c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", "Bleached + Heated"))

#visualize
ggplot(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",],aes(x=Treatment_v1,y=sym.SA.normalized.no.outliers,color=Treatment_v1,fill=Treatment_v1))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  geom_boxplot(size=2)+
  # geom_point(size=3)+
  scale_color_manual( values=c("dodgerblue3","firebrick3","dodgerblue3","firebrick3"))+
  scale_fill_manual(values=c("dodgerblue1","firebrick1","white","white"))+
  scale_x_discrete( guide = guide_axis(n.dodge = 2), name = "")+
  theme_classic()+
  theme(text=element_text(size=24),legend.key.height=unit(2,"cm"))+
  labs(y="Symbiodiniaceae cells per cm^2",x="Treatment",color="Treatment",fill="Treatment")

#visualize with log10 densities. Use for panel 1D
ggplot(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",],aes(x=Treatment_v1,y=log10(sym.SA.normalized.no.outliers),color=Treatment_v1,fill=Treatment_v1))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  geom_boxplot(size=2)+
  # geom_point(size=3)+
  scale_color_manual(values=cost.col.line, labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_fill_manual(values=cost.col.fill, labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Bleached + Heated"), name = "")+
  theme_classic()+
  theme(text=element_text(size=34),legend.key.height=unit(2,"cm"),axis.text.x=element_blank(), plot.margin=unit(c(1,1,1,1.5), "cm"))+
  ylab(expression(paste("Log10 Symbiodiniaceae\ncells per cm^2")))+
  labs(color="Treatment",fill="Treatment")
ggsave('Symbiont cells per cm2_per treatment v2.jpeg', path = dirFigs, dpi = 300, width=15, height=9)

#visualize again with only PLANC aquaria that were included in abcDOM
ggplot(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7" & is.na(sym_dat_aquaria_final$Sample.Name)==FALSE,],aes(x=Treatment_v1,y=sym.SA.normalized.no.outliers,color=Treatment_v1,fill=Treatment_v1))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  geom_boxplot(size=2)+
  # geom_point(size=3)+
  scale_color_manual(values=c("dodgerblue3","firebrick3","dodgerblue3","firebrick3"),  labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_fill_manual( values=c("dodgerblue1","firebrick1","white","white"),  labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Bleached + Heated"), name = "")+
  theme_classic()+
  theme(text=element_text(size=24),legend.key.height=unit(2,"cm"))+
  labs(y="Symbiodiniaceae cells per cm^2",x="Treatment",color="Treatment",fill="Treatment")

#now with log10 sym densities
ggplot(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7" & is.na(sym_dat_aquaria_final$Sample.Name)==FALSE,],aes(x=Treatment_v1,y=log10(sym.SA.normalized.no.outliers),color=Treatment_v1,fill=Treatment_v1))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  geom_boxplot(size=2)+
  # geom_point(size=3)+
  scale_color_manual( values=c("dodgerblue3","firebrick3","dodgerblue3","firebrick3"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_fill_manual( values=c("dodgerblue1","firebrick1","white","white"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_x_discrete( name = "", guide = guide_axis(n.dodge = 2),labels = c("Control", "Heated", "Bleached", "Bleached + Heated") )+
  theme_classic()+
  theme(text=element_text(size=24),legend.key.height=unit(2,"cm"))+
  labs(y="log10 mean aquaria Symbiodiniaceae cells per cm^2",x="Treatment",color="Treatment",fill="Treatment")
#save and
ggsave('Symbiont cells per cm2_per treatment v1.jpeg', path = dirFigs, dpi = 300, width=15, height=9)

#next, plot a line graph of sym densities for synthesis figure
sym_dat_aquaria_final_t7 <- sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",]
sym_dat_aquaria_final_t7_mean <- aggregate(sym_dat_aquaria_final_t7[,7], by=list(sym_dat_aquaria_final_t7$Treatment_v1), FUN=mean) #calculate mean
sym_dat_aquaria_final_t7_mean_v1 <- cbind(sym_dat_aquaria_final_t7_mean, aggregate(sym_dat_aquaria_final_t7[,7], by=list(sym_dat_aquaria_final_t7$Treatment_v1), FUN=stderror)[,2])
colnames(sym_dat_aquaria_final_t7_mean_v1) <- c("Treatment", "mean_sym.SA.normalized.no.outliers", "se")
sym_dat_aquaria_final_t7_mean_v1$Treatment <- factor(sym_dat_aquaria_final_t7_mean_v1$Treatment,levels=c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Heated", "Bleached + Ambient"))

ggplot(sym_dat_aquaria_final_t7_mean_v1, aes(color=Treatment, fill=Treatment, x=Treatment, y=log10(mean_sym.SA.normalized.no.outliers)))+
  geom_point(size=10, shape=21)+
  geom_pointrange(aes(x=Treatment, y=log10(mean_sym.SA.normalized.no.outliers), ymin=log10(mean_sym.SA.normalized.no.outliers-se), ymax=log10(mean_sym.SA.normalized.no.outliers+se)))+
  scale_color_manual(values=c("dodgerblue3","firebrick3","firebrick3","dodgerblue3"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_fill_manual(values=c("dodgerblue1","firebrick1","NA","NA","NA","NA"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  ylab("log10 Symbiodiniaceae cells per cm^2")+
  theme_classic()+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position="none")
ggsave('Symbiont cells per cm2_per treatment dotplot.png', path = dirFigs, dpi = 600, width=12, height=3)








# Now, visualize the distribution and run the statistics on the t7 aquaria data with ONLY abcDOM samples.

hist(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7" & is.na(sym_dat_aquaria_final$Sample.Name)==FALSE,]$sym.SA.normalized.no.outliers)#normal-ish
hist(log10(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7" & is.na(sym_dat_aquaria_final$Sample.Name)==FALSE,]$sym.SA.normalized.no.outliers))#normal-ish

# Since the data look normal enough, test the effect of Treatment_v1 using an ANOVA on both log10 and raw dat.
#run on log10 data
mod.t7.sym=aov(log10(sym.SA.normalized.no.outliers) ~ Treatment_v1, data=sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7" & is.na(sym_dat_aquaria_final$Sample.Name)==FALSE,])
summary(mod.t7.sym) #not significant

#run on raw data
mod.t7.sym.raw=aov(sym.SA.normalized.no.outliers ~ Treatment_v1, data=sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7" & is.na(sym_dat_aquaria_final$Sample.Name)==FALSE,])
summary(mod.t7.sym.raw) #marginally significant

#run two wat anova on log10 data
mod.t7.sym.2way=aov(log10(sym.SA.normalized.no.outliers) ~ Stress_status*Bleaching_Status, data=sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7" & is.na(sym_dat_aquaria_final$Sample.Name)==FALSE,])
summary(mod.t7.sym.2way) #bleaching status significant, stress status (heated) not significant.

#run two wat anova on raw data
mod.t7.sym.raw.2way=aov(sym.SA.normalized.no.outliers ~ Stress_status*Bleaching_Status, data=sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7" & is.na(sym_dat_aquaria_final$Sample.Name)==FALSE,])
summary(mod.t7.sym.raw.2way) #bleaching status significant, stress status (heated) not significant.

#Now repeat for t7 aquaria data with ALL samples.

hist(log10(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",]$sym.SA.normalized.no.outliers))#normal-ish
hist(sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",]$sym.SA.normalized.no.outliers)#normal-ish

# Since the data look normal enough, test the effect of Treatment_v1 using a one way and two way ANOVA on raw and log10 transformed values
#first run on log10 values
mod.t7.sym.log.1way=aov(log10(sym.SA.normalized.no.outliers) ~ Treatment_v1, data=sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",])
summary(mod.t7.sym.log.1way) #not significant
#run tukeyhsd
TukeyHSD(mod.t7.sym.log.1way, "Treatment_v1")

mod.t7.sym.log.2way=aov(log10(sym.SA.normalized.no.outliers) ~ Stress_status*Bleaching_Status, data=sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",])
summary(mod.t7.sym.log.2way) #bleaching status is significant.

#next run on raw values
mod.t7.sym.raw.1way=aov(sym.SA.normalized.no.outliers ~ Treatment_v1, data=sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",])
summary(mod.t7.sym.raw.1way) #significant

mod.t7.sym.raw.2way=aov(sym.SA.normalized.no.outliers ~ Stress_status*Bleaching_Status, data=sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",])
summary(mod.t7.sym.raw.2way) #bleaching status is significant.

kruskal.test(sym.SA.normalized.no.outliers ~ Treatment_v1, data=sym_dat_aquaria_final[sym_dat_aquaria_final$Timepoint_char=="T7",])
#marginally significant.

#Next, repeat visualization and analysis of the t7 data but by nubbin (merged and seperated by species).

#First, create a treatment column in the sym_dat df.
sym_dat$Treatment <- paste(sym_dat$Collection_Bleaching_Level1, sym_dat$Heat, sep="_") #merge bleaching and heat treatments
sym_dat$Treatment[sym_dat$Treatment=="HE_C"] <- "Non-bleached + Ambient"
sym_dat$Treatment[sym_dat$Treatment=="HE_H"] <- "Non-bleached + Heated"
sym_dat$Treatment[sym_dat$Treatment=="BL_C"] <- "Bleached + Ambient"
sym_dat$Treatment[sym_dat$Treatment=="BL_H"] <- "Bleached + Heated"
sym_dat$Treatment <- factor(sym_dat$Treatment, levels=c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", "Bleached + Heated")) #update levels

#Next, create a species_treatment column
sym_dat$species_treatment <- paste(sym_dat$Species, sym_dat$Treatment, sep="_")

#next, add a column if the nubbins were used in abcDOM
sym_dat$ABCDOM <- sym_dat$Aquaria_Timepoint_Heat
sym_dat$ABCDOM[sym_dat$ABCDOM %in% metadata1$PLANC_aquaria] = "T"
sym_dat$ABCDOM[sym_dat$ABCDOM!="T"] = "F"

#Visualize sym.SA by treatment for nubbins at T7
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T7" & is.na(sym_dat$Species)!=TRUE,],(aes(x=Treatment,y=sym.SA,color=Treatment, fill=Treatment)))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  geom_boxplot(size=2)+
  scale_color_manual( values=c("dodgerblue3","firebrick3","dodgerblue3","firebrick3"),labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_fill_manual( values=c("dodgerblue1","firebrick1","white","white"),labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_x_discrete( name = "", guide = guide_axis(n.dodge = 2),labels = c("Control", "Heated", "Bleached", "Bleached + Heated") )+
  theme_classic()+
  theme(text=element_text(size=24),legend.key.height=unit(2,"cm"))+
  labs(y="Symbiodiniaceae cells per cm^2",x="Treatment",color="Treatment",fill="Treatment")

#Visualize log10sym.SA by treatment for nubbins at T7
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T7" & is.na(sym_dat$Species)!=TRUE,],(aes(x=Treatment,y=log10(sym.SA),color=Treatment, fill=Treatment)))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  geom_boxplot(size=2)+
  scale_color_manual( values=c("dodgerblue3","firebrick3","dodgerblue3","firebrick3"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_fill_manual( values=c("dodgerblue1","firebrick1","white","white"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_x_discrete( name = "", guide = guide_axis(n.dodge = 2), labels = c("Control", "Heated", "Bleached", "Bleached + Heated") )+
  theme_classic()+
  theme(text=element_text(size=24),legend.key.height=unit(2,"cm"))+
  labs(y="log10 Symbiodiniaceae cells per cm^2",x="Treatment",color="Treatment",fill="Treatment")

#Visualize sym.SA by treatment for nubbins at T7, faceted by species
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T7" & is.na(sym_dat$Species)!=TRUE,],(aes(x=Treatment,y=sym.SA,color=Treatment, fill=Treatment)))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  facet_wrap(.~Species)+
  geom_boxplot(size=2)+
  scale_color_manual( values=c("dodgerblue3","firebrick3","dodgerblue3","firebrick3"),labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_fill_manual( values=c("dodgerblue1","firebrick1","white","white"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_x_discrete( name = "", guide = guide_axis(angle = 90), labels = c("Control", "Heated", "Bleached", "Bleached + Heated") )+
  theme_classic()+
  theme(text=element_text(size=24),legend.key.height=unit(2,"cm"))+
  labs(y="log10 mean aquaria Symbiodiniaceae cells per cm^2",x="",color="Treatment",fill="Treatment")

#Visualize log10sym.SA by treatment for nubbins at T7, faceted by species
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T7" & is.na(sym_dat$Species)!=TRUE,],(aes(x=Treatment,y=log10(sym.SA),color=Treatment, fill=Treatment)))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  facet_wrap(.~Species)+
  geom_boxplot(size=2)+
  scale_color_manual( values=c("dodgerblue3","firebrick3","dodgerblue3","firebrick3"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_fill_manual( values=c("dodgerblue1","firebrick1","white","white"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_x_discrete( name = "", guide = guide_axis(angle = 90), labels = c("Control", "Heated", "Bleached", "Bleached + Heated") )+
  theme_classic()+
  theme(text=element_text(size=24),legend.key.height=unit(2,"cm"))+
  labs(y="log10 mean aquaria Symbiodiniaceae cells per cm^2",x="",color="Treatment",fill="Treatment")

#Visualize sym.SA by treatment for only abcDOM nubbins at T7
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T7" & is.na(sym_dat$Species)!=TRUE & sym_dat$ABCDOM=="T",],(aes(x=Treatment,y=sym.SA,color=Treatment, fill=Treatment)))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  geom_boxplot(size=2)+
  scale_color_manual( values=c("dodgerblue3","firebrick3","dodgerblue3","firebrick3"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_fill_manual( values=c("dodgerblue1","firebrick1","white","white"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_x_discrete( name = "", guide = guide_axis(n.dodge = 2), labels = c("Control", "Heated", "Bleached", "Bleached + Heated") )+
  theme_classic()+
  theme(text=element_text(size=24),legend.key.height=unit(2,"cm"))+
  labs(y="Symbiodiniaceae cells per cm^2",x="Treatment",color="Treatment",fill="Treatment")

#Visualize log10sym.SA by treatment for only abcDOM nubbins at T7
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T7" & is.na(sym_dat$Species)!=TRUE & sym_dat$ABCDOM=="T",],(aes(x=Treatment,y=log10(sym.SA),color=Treatment, fill=Treatment)))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  geom_boxplot(size=2)+
  scale_color_manual( values=c("dodgerblue3","firebrick3","dodgerblue3","firebrick3"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_fill_manual( values=c("dodgerblue1","firebrick1","white","white"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_x_discrete( name = "", guide = guide_axis(n.dodge = 2), labels = c("Control", "Heated", "Bleached", "Bleached + Heated") )+
  theme_classic()+
  theme(text=element_text(size=24),legend.key.height=unit(2,"cm"))+
  labs(y="Symbiodiniaceae cells per cm^2",x="Treatment",color="Treatment",fill="Treatment")

#Visualize sym.SA by treatment for only abcDOM nubbins at T7, faceted by species
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T7" & is.na(sym_dat$Species)!=TRUE & sym_dat$ABCDOM=="T",],(aes(x=Treatment,y=sym.SA,color=Treatment, fill=Treatment)))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  facet_wrap(.~Species)+
  geom_boxplot(size=2)+
  scale_color_manual( values=c("dodgerblue3","firebrick3","dodgerblue3","firebrick3"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_fill_manual( values=c("dodgerblue1","firebrick1","white","white"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_x_discrete( name = "", guide = guide_axis(angle = 90), labels = c("Control", "Heated", "Bleached", "Bleached + Heated") )+
  theme_classic()+
  theme(text=element_text(size=24),legend.key.height=unit(2,"cm"))+
  labs(y="Symbiodiniaceae cells per cm^2",x="Treatment",color="Treatment",fill="Treatment")

#Visualize log10sym.SA by treatment for only abcDOM nubbins at T7, faceted by species
ggplot(sym_dat[sym_dat$Outlier!="Y" & sym_dat$Timepoint=="T7" & is.na(sym_dat$Species)!=TRUE & sym_dat$ABCDOM=="T",],(aes(x=Treatment,y=log10(sym.SA),color=Treatment, fill=Treatment)))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  facet_wrap(.~Species)+
  geom_boxplot(size=2)+
  scale_color_manual( values=c("dodgerblue3","firebrick3","dodgerblue3","firebrick3"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_fill_manual( values=c("dodgerblue1","firebrick1","white","white"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_x_discrete( name = "", guide = guide_axis(angle = 90), labels = c("Control", "Heated", "Bleached", "Bleached + Heated") )+
  theme_classic()+
  theme(text=element_text(size=24),legend.key.height=unit(2,"cm"))+
  labs(y="Symbiodiniaceae cells per cm^2",x="Treatment",color="Treatment",fill="Treatment")

#After this, it still looks better to present things as an aquaria wide average.
