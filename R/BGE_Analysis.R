'BGE Analysis.R'

#load libraries
#library(ggplot2)

#Work up BGE by subsetting t0 and tfinal FCM data
FCM_dat <- FCM_dat[is.na(FCM_dat)!=TRUE,] #remove NAs
BGE <- subset(FCM_dat, Timepoint..h.==0) #subset FCM dat for just t0, save as BGE df
FCM_dat_t36_concentration <- FCM_dat$Concentration[FCM_dat$Timepoint..h.==36] #subset for t36 concentrations
FCM_dat_t36_concentration1 <- FCM_dat_t36_concentration[is.na(FCM_dat_t36_concentration)!=TRUE] #remove NAs
BGE$Concentration_t36 <- FCM_dat_t36_concentration1 #add in t36 concentrations to BGE df.
colnames(BGE)[3] <- "Concentration_t0" #update colnames

#Work up BGE by adding DOC data for t0 and tfinal
BGE1 <- merge(BGE, DOC_dat1[DOC_dat1$Timepoint==0,c(5,18)], by.x="Bottle", by.y="Bottle_NR.x", all.x=T, all.y=F) #merge BGE data with t0 DOC data, by bottle number. Merging allows to include NAs for missing DOC data.
BGE2 <- merge(BGE1, DOC_dat1[DOC_dat1$Timepoint==36,c(5,18)], by.x="Bottle", by.y="Bottle_NR.x", all.x=T, all.y=F) #repeat for t36
colnames(BGE2)[23:24] <- c("DOC_t0", "DOC_t36")

#add column corresponding to bottle volume.
BGE2$Bottle_volume_mL <- c(1200, 1180, 1200, 1200, 1250, 1370, 1010, 1350, 1220, 1200, 900, 1200, 1250, 1250, 1400, 1400, 1400)

#calculate changes in cell c from t0 to t36.
BGE2$FCM_concentration_change <- BGE2$Concentration_t36 - BGE2$Concentration_t0 #calculate concentratin change
BGE2$FCM_cell_change <- BGE2$FCM_concentration_change*BGE2$Bottle_volume_mL*1000 #calculate cell count change
BGE2$FCM_C_change_fg <- BGE2$FCM_cell_change*20 #calculate change in cell C using a conversion ration of 20fgC per cell
BGE2$FCM_C_change_uG <- BGE2$FCM_C_change_fg*10^-9 #convert to cell C in uG.

#calculate changes in DOC (g) from t0 to t36.
BGE2$DOC_change <- BGE2$DOC_t0 - BGE2$DOC_t36 #calculate change in DOC concentration
BGE2$DOC_change_umolC <- (BGE2$DOC_change*BGE2$Bottle_volume_mL)/1000 #calculate DOC change in uMol C.
BGE2$DOC_change_uG <- BGE2$DOC_change_umolC*12.0107 #convert to DOC change in uG.

#calculate BGE
BGE2$BGE <- BGE2$FCM_C_change_uG/BGE2$DOC_change_uG

#update factor levels
BGE2$Treatment <- factor(BGE2$Treatment.y, levels=levels(fact.all.treat))

#visualize
ggplot(BGE2[BGE2$Treatment!="Bleached + Heated",], aes(x=Treatment, y=BGE, color=Treatment, fill=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(size = 1.2)+
  # geom_point(size = 3)+
  scale_color_manual(values=cost.col.line[-4])+
  scale_fill_manual(values=cost.col.fill[-4], guide = guide_legend(override.aes = list(size = 1)))+
  # theme(legend.key.height=unit(0.5,"in"))+
  #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylab ("Bacterial Growth Efficiancy")+
  xlab("")
ggsave('BGE_per_treatment.jpeg', path = dirFigs, width = 11, height = 8, dpi = 300)

#run stats on BGE data
hist(BGE2$BGE) #check distribution. Looks normal ish.

mod.BGE <- aov(BGE ~ Treatment, data=BGE2[BGE2$Treatment!="Bleached + Heated",]) #run the model, removing the n=1 bleached+heated
summary(mod.BGE) #Treatment has a marginally significant effect.
