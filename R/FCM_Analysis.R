'FCM_Analysis.R '

#load libraries
library(ggplot2)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)

# Check distribution of concentration data.
hist(FCM_dat$Concentration) #not normal, strong right skew
hist(sqrt(FCM_dat$Concentration)) #more normal

#calculate the mean and SE for bleaching susceptability treatments and timepoint
# but first create the SE function
stderror=function(x){
  se=sd(x)/sqrt(length(x))
  return(se)
}

FCM_dat_mean <- as.data.frame(aggregate(FCM_dat$Concentration,by=list(paste(FCM_dat$Treatment.y,FCM_dat$Timepoint..h.,sep="_")),FUN=mean)) #calcualte mean

FCM_dat_mean$SE <- as.data.frame(aggregate(FCM_dat$Concentration,by=list(paste(FCM_dat$Treatment.y,FCM_dat$Timepoint..h.,sep="_")),FUN=stderror))[[2]] #add se data

colnames(FCM_dat_mean) <- c("Sample_Timepoint", "Mean_Concentration", "SE") #adjust colnames

FCM_dat_mean$Timepoint <- unlist(strsplit(FCM_dat_mean$Sample_Timepoint,split="_"))[seq(from=2,to=110,by=2)] #extract timepoint values and add to new column
FCM_dat_mean$Timepoint = as.numeric(FCM_dat_mean$Timepoint)

FCM_dat_mean$Treatment <- unlist(strsplit(FCM_dat_mean$Sample_Timepoint,split="_"))[seq(from=1,to=110,by=2)] #extract timepoint values and add to new column

FCM_dat_mean$Treatment <- factor(FCM_dat_mean$Treatment, levels=c("Non-bleached + Ambient","Non-bleached + Heated","Bleached + Ambient","Bleached + Heated", "Ambient Water Control","Heated Water Control")) #set treatment as a factor

FCM_dat_mean <- FCM_dat_mean[-37,] #remove NA row

# Visualize.

#calculate maximum concentration for each treatment at t 24
FCM_dat_mean_24 <- subset(FCM_dat_mean, Timepoint==24) #subset for t24
FCM_dat_mean_24$posthoc <- c("C/D","B","A/B","D","B/C","A")

ggplot(FCM_dat_mean,aes(x=Timepoint,y=Mean_Concentration,fill=Treatment,shape=Treatment,color=Treatment,group=Treatment))+
  geom_point(size=10)+
  geom_line(size=3)+
  scale_shape_manual(values=c(21,21,21,21,21,21), labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  geom_pointrange(aes(ymin=Mean_Concentration-SE,ymax=Mean_Concentration+SE))+
  scale_color_manual(values=cost.col.line, labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  scale_fill_manual(values=cost.col.fill, labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  theme_classic()+
  theme(text=element_text(size=30),legend.key.height=unit(1.75,"cm"),complete=FALSE)+
  ylab(label="Concentration (cells per uL)")+
  xlab(label="Time (hours)")+
  geom_text(data=FCM_dat_mean_24, aes(x=Timepoint, y=Mean_Concentration+120, hjust=1.25, label=posthoc), color="black", size=8)
ggsave('Microbialgrowthcurve_per_treatment.jpeg', path = dirFigs, dpi = 300, height=10, width=14)

fig2C<-ggplot(FCM_dat_mean,aes(x=Timepoint,y=Mean_Concentration,fill=Treatment,shape=Treatment,color=Treatment,group=Treatment))+
  geom_point(size=8)+
  geom_line(size=3)+
  scale_shape_manual(values=c(21,21,21,21,21,21), labels = c("Control [C]", "Heated [A]", "Bleached [B/C]", "Bleached + Heated [B]", "Negative Control [D]", "Negative Control + Heated [E]"))+
  geom_pointrange(aes(ymin=Mean_Concentration-SE,ymax=Mean_Concentration+SE))+
  scale_color_manual(values=cost.col.line, labels = c("Control [C]", "Heated [A]", "Bleached [B/C]", "Bleached + Heated [B]", "Negative Control [D]", "Negative Control + Heated [E]"))+
  scale_fill_manual(values=cost.col.fill, labels = c("Control [C]", "Heated [A]", "Bleached [B/C]", "Bleached + Heated [B]", "Negative Control [D]", "Negative Control + Heated [E]"))+
  theme_classic()+
  theme(legend.position="right")+
  ylab(label="Concentration (cells/µL)")+
  xlab(label="Time (hours)")
  #geom_text(data=FCM_dat_mean_24, aes(x=Timepoint, y=Mean_Concentration+120, hjust=1.25, label=posthoc), color="black", size=5)
fig2C


#Visualize without water controls.

ggplot(subset(FCM_dat_mean, Treatment != "Ambient Water Control" & Treatment != "Heated Water Control"),aes(x=Timepoint,y=Mean_Concentration,fill=Treatment,shape=Treatment,color=Treatment,group=Treatment))+
  geom_point(size=10)+
  geom_line(size=3)+
  scale_shape_manual(values=c(21,21,21,21,21,21), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  geom_pointrange(aes(ymin=Mean_Concentration-SE,ymax=Mean_Concentration+SE))+
  scale_color_manual(values=c("dodgerblue3","firebrick3","dodgerblue3","firebrick3","black","gray"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  scale_fill_manual(values=c("dodgerblue1","firebrick1","NA","NA","NA","NA"), labels = c("Control", "Heated", "Bleached", "Bleached + Heated"))+
  theme_classic()+
  theme(text=element_text(size=30),legend.key.height=unit(1.75,"cm"),complete=FALSE)+
  ylab(label="Concentration (cells per uL)")+
  xlab(label="Time (hours)")

#Run statistics on the 24h timepoint with all treatments.
#subset FCM dat for just t24.
FCM_dat_24 <- subset(FCM_dat, Timepoint..h. == 24)

#check distribution
hist(FCM_dat_24$Concentration) #not super normal
hist(sqrt(FCM_dat_24$Concentration)) #more normal

#add the sqrt transformed concentratio data as a new column
FCM_dat_24$sqrt_concentration <- sqrt(FCM_dat_24$Concentration)

#run an ANOVA on the sqrt transformed concentration data
mod.24.sqrt <- aov(sqrt_concentration ~ Treatment.y, data=FCM_dat_24)
summary(mod.24.sqrt) #treatment is significant

#run a tukey HSD on mod.24.sqrt
TukeyHSD(mod.24.sqrt, "Treatment.y")

#repeat without water controls.
mod.24.sqrt.nowater <- aov(Concentration ~ Treatment.y, data=subset(FCM_dat_24, Origin_PlanC != "control"))
summary(mod.24.sqrt.nowater) #marginally significant

#Run statistics again with raw concentration values.
#run an ANOVA on the sqrt transformed concentration data
mod.24 <- aov(Concentration ~ Treatment.y, data=FCM_dat_24)
summary(mod.24) #treatment is significant

#run a tukey HSD on mod.24
TukeyHSD(mod.24, "Treatment.y")
#For now we will stick with these stats since in the post hoc there are significant differences within the coral treatment.

#run a test on the 24 data without neg controls
mod.24.nowater <- aov(Concentration ~ Treatment.y, data=subset(FCM_dat_24, Origin_PlanC != "control"))
summary(mod.24.nowater) #treatment is marginally significant. Don't bother using tests without the water controls.

#visualize T0 vs Tmax(24h) cell concentrations for Fig 1
FCM_dat_24_0 <- subset(FCM_dat, Timepoint..h. == 0 | Timepoint..h. == 24) #subset for the appropriate timepoints

FCM_dat_24_0$Treatment <- FCM_dat_24_0$Treatment.y #duplicate Treatment.y

FCM_dat_24_0$Treatment <- factor(FCM_dat_24_0$Treatment, levels(fact.all.treat)) #set as correct factor levels

ggplot(FCM_dat_24_0, aes(x=as.factor(Timepoint..h.), y=Concentration, color=Treatment, fill=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  geom_boxplot(size = 2)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  ylab(label="Concentration (cells per uL)")+
  xlab(label="Time of dark incubation (hours)")+
  theme_classic()+
  theme(text=element_text(size=34),legend.key.height=unit(2,"cm"),complete=FALSE)
ggsave('Microbialgrowth_after_24h_per_treatment.jpeg', path = dirFigs, dpi = 300, width=15, height = 13)

#run stats on 24_0 data, just testing for timepoint.
mod.24.0 <- aov(Concentration ~ Timepoint..h., data=FCM_dat_24_0)
summary(mod.24.0)

#calculate maximum SGR for each treatment.
FCM_dat_growth_max <- aggregate(FCM_dat_growth[,25], by=list(FCM_dat_growth$Treatment), FUN=max, na.rm=T) #be sure to exclude NAs
FCM_dat_growth_max$posthoc <- c("A/B","A","A","A","A/B","B")
colnames(FCM_dat_growth_max)[1:2] <- c("Treatment", "max")

#visualize specific growth rate.
ggplot(FCM_dat_growth[!is.na(FCM_dat_growth$Specific_Growth_Rate),],aes(x=Treatment,y=Specific_Growth_Rate,color=Treatment,fill=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  geom_boxplot(size = 2)+
  # geom_point(size = 3)+
  scale_color_manual(values=cost.col.line, labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  scale_fill_manual(values=cost.col.fill, labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  theme_classic()+
  # theme(legend.key.height=unit(2,"cm"))+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  ylab ("Specific Growth Rate (log10 cells per hour)")+
  xlab("")+
  theme(text=element_text(size=24), legend.position = "none")+
  geom_text(data=FCM_dat_growth_max, aes(x=Treatment, y=max+.005, label=posthoc), color="black", size=6)
  # annotate("rect", xmin=-Inf, xmax=-Inf, ymin=, ymax =)
ggsave('Specific_growth_rate.jpeg', path = dirFigs, dpi = 300, height=9, width=15)

fig2D<-ggplot(FCM_dat_growth[!is.na(FCM_dat_growth$Specific_Growth_Rate),],aes(x=Treatment,y=Specific_Growth_Rate,color=Treatment,fill=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  geom_boxplot(size = 2)+
  # geom_point(size = 3)+
  scale_color_manual(values=cost.col.line, labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  scale_fill_manual(values=cost.col.fill, labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  theme_classic()+
  # theme(legend.key.height=unit(2,"cm"))+
  scale_x_discrete(labels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  ylab ("Specific Growth Rate (log10 cells per hour)")+
  xlab("")+
  theme(legend.position = "none")+
  geom_text(data=FCM_dat_growth_max, aes(x=Treatment, y=max+.005, label=posthoc), color="black", size=5)
fig2D

#now make figure 2 complete
fig2<-plot_grid(fig2A, fig2B, fig2C, fig2D, labels="AUTO", nrow =2, rel_widths = c(1, 1.4), align = c("hv"), axis = "lt")
fig2
ggsave("figure2.jpg", fig2, path = dirFigs, dpi = 300, height=12, width=15)
fig2_v1<-plot_grid(fig2A, fig2C, labels="AUTO", nrow =2, rel_widths = c(1, 1.4), align = c("hv"), axis = "lt")
fig2_v1
ggsave("figure2_v1.jpg", fig2_v1, path = dirFigs, dpi = 300, height=10, width=7)
fig2_v2<-plot_grid(fig2A_v1, fig2C, labels="AUTO", nrow =2, rel_widths = c(1, 1.4), align = c("hv"), axis = "lt")
fig2_v2
ggsave("figure2_v2.jpg", fig2_v2, path = dirFigs, dpi = 300, height=10, width=7)
fig2_v3<-plot_grid(fig2A_v1, fig2C, labels="AUTO", nrow =1, rel_widths = c(1, 2), align = c("hv"), axis = "lt")
fig2_v3
ggsave("figure2_v3.jpg", fig2_v3, path = dirFigs, dpi = 300, height=6, width=13)

#run stats on specific growth rate.
hist(FCM_dat_growth$Specific_Growth_Rate) #check distribution. Looks normal emough.

mod.growth.rate <- aov(Specific_Growth_Rate ~ Treatment, data=FCM_dat_growth) #run model
summary(mod.growth.rate) #examine model
TukeyHSD(mod.growth.rate, "Treatment") #run posthoc

#run mixed mod on fcm data
FCM_dat$sqrt_concentration <- sqrt(FCM_dat$Concentration)
FCM_dat_v1 <- FCM_dat[is.na(FCM_dat$Timepoint..h.)==FALSE,] #remove controls

#for full dataset
mod.fcm.lm <- lmer(sqrt_concentration ~ Treatment.y*Timepoint..h. + (1|Bottle_NR), data=FCM_dat)
summary(mod.fcm.lm, ddf="Kenward-Roger") #ranef explains 2.5x compared to residiuals (fixed)
anova(mod.fcm.lm, ddf="Kenward-Roger") #Time and Treatment*Time are significant, no treatment on its own.
mod.fcm.lm.ph <- lsmeans(mod.fcm.lm, "Treatment.y", adjust="tukey") #Post Hoc analysis for the Gut_SubSection factor. This function averages across Species so that it is ONLY testing the differences between Gut_SubSection.
contrast(mod.fcm.lm.ph, alpha=.05, method="pairwise", adjust="tukey")
mod.fcm.lm.ph1 <- lsmeans(mod.fcm.lm, ~Treatment.y*Timepoint..h., adjust="tukey") #Post Hoc analysis for the Gut_SubSection factor. This function averages across Species so that it is ONLY testing the differences between Gut_SubSection.
contrast(mod.fcm.lm.ph1, alpha=.05, method="pairwise", adjust="tukey")

#for data just up to 24 hours (peak growth)
mod.fcm.lm.v1 <- lmer(sqrt_concentration ~ Treatment.y*Timepoint..h. + (1|Bottle_NR), data=subset(FCM_dat, Timepoint..h.!=32 & Timepoint..h.!=36))
summary(mod.fcm.lm.v1, ddf="Kenward-Roger") #ranef explains similar variation compared to residiuals (fixed)
anova(mod.fcm.lm.v1, ddf="Kenward-Roger") #Time and Treatment*Time are significant, no treatment on its own.

#try just using lm
mod.fcm.lm.v3 <- aov(sqrt_concentration ~ Treatment.y*Timepoint..h., data=FCM_dat)
anova(mod.fcm.lm.v3) #Time and Treatment*Time and Ttreatment are significant.
TukeyHSD(mod.fcm.lm.v3, "Treatment.y")

#for just coral treatments
mod.fcm.lm.justcoral <- lmer(sqrt_concentration ~ Treatment.y*Timepoint..h. + (1|Bottle_NR), data=subset(FCM_dat_v1, Origin_PlanC!="control"))
summary(mod.fcm.lm.justcoral, ddf="Kenward-Roger") #ranef explains 2.5x compared to residiuals (fixed)
anova(mod.fcm.lm.justcoral, ddf="Kenward-Roger") #Time and Treatment*Time are significant, no treatment on its own.

#try just using lm with timepoint as a factor
mod.fcm.lm.v3.factor <- aov(sqrt_concentration ~ Treatment.y*as.factor(Timepoint..h.), data=FCM_dat)
anova(mod.fcm.lm.v3.factor) #Time and Treatment*Time and Ttreatment are significant. Use this model for final stats!!
TukeyHSD(mod.fcm.lm.v3.factor, "Treatment.y")


