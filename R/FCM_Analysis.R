'FCM_Analysis.R '

#load libraries
library(ggplot2)

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
FCM_dat_mean_24$posthoc <- c("B/C","A/B/C","A/B","C","B/C","A")

ggplot(FCM_dat_mean,aes(x=Timepoint,y=Mean_Concentration,fill=Treatment,shape=Treatment,color=Treatment,group=Treatment))+
  geom_point(size=10)+
  geom_line(size=3)+
  scale_shape_manual(values=c(21,21,21,21,21,21))+
  geom_pointrange(aes(ymin=Mean_Concentration-SE,ymax=Mean_Concentration+SE))+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill)+
  theme_classic()+
  theme(text=element_text(size=30),legend.key.height=unit(1.75,"cm"),complete=FALSE)+
  ylab(label="Concentration (cells per uL)")+
  xlab(label="Time (hours)")+
  geom_text(data=FCM_dat_mean_24, aes(x=Timepoint, y=Mean_Concentration+120, hjust=1.25, label=posthoc), color="black", size=3.5)
ggsave('Microbialgrowthcurve_per_treatment.jpeg', path = dirFigs, dpi = 300, height=10, width=14)

#Visualize without water controls.

ggplot(subset(FCM_dat_mean, Treatment != "Ambient Water Control" & Treatment != "Heated Water Control"),aes(x=Timepoint,y=Mean_Concentration,fill=Treatment,shape=Treatment,color=Treatment,group=Treatment))+
  geom_point(size=10)+
  geom_line(size=3)+
  scale_shape_manual(values=c(21,21,21,21,21,21))+
  geom_pointrange(aes(ymin=Mean_Concentration-SE,ymax=Mean_Concentration+SE))+
  scale_color_manual(values=c("#00BFC4","#F8766D","#00BFC4","#F8766D","black","gray"))+
  scale_fill_manual(values=c("dodgerblue3","firebrick3","NA","NA","NA","NA"))+
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
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  theme_classic()+
  # theme(legend.key.height=unit(2,"cm"))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylab ("Specific Growth Rate (log10 cells per hour)")+
  xlab("")+
  theme(text=element_text(size=24))+
  geom_text(data=FCM_dat_growth_max, aes(x=Treatment, y=max+.005, label=posthoc), color="black", size=6)+
  annotate("rect", xmin=-Inf, xmax=-Inf, ymin=)
ggsave('Specific_growth_rate.jpeg', path = dirFigs, dpi = 300, height=9, width=15)

#run stats on specific growth rate.
hist(FCM_dat_growth$Specific_Growth_Rate) #check distribution. Looks normal emough.

mod.growth.rate <- aov(Specific_Growth_Rate ~ Treatment, data=FCM_dat_growth) #run model
summary(mod.growth.rate) #examine model
TukeyHSD(mod.growth.rate, "Treatment") #run posthoc
