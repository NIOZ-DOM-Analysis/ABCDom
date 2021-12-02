'FCM_Analysis.R '

#make a function to calculate standard error
stderror=function(x) {
  se=sd(x)/sqrt(length(x))
  return(se)
}


# Check distribution of concentration data.
hist(FCM_dat$Concentration) #not normal, strong right skew
hist(sqrt(FCM_dat$Concentration)) #more normal

#calculate the mean and SE for bleaching susceptability treatments and timepoint
FCM_dat_mean <- as.data.frame(aggregate(FCM_dat$Concentration,by=list(paste(FCM_dat$Treatment.y,FCM_dat$Timepoint..h.,sep="_")),FUN=mean)) #calcualte mean

FCM_dat_mean$SE <- as.data.frame(aggregate(FCM_dat$Concentration,by=list(paste(FCM_dat$Treatment.y,FCM_dat$Timepoint..h.,sep="_")),FUN=stderror))[[2]] #add se data

colnames(FCM_dat_mean) <- c("Sample_Timepoint", "Mean_Concentration", "SE") #adjust colnames

FCM_dat_mean$Timepoint <- unlist(strsplit(FCM_dat_mean$Sample_Timepoint,split="_"))[seq(from=2,to=110,by=2)] #extract timepoint values and add to new column
FCM_dat_mean$Timepoint = as.numeric(FCM_dat_mean$Timepoint)

FCM_dat_mean$Treatment <- unlist(strsplit(FCM_dat_mean$Sample_Timepoint,split="_"))[seq(from=1,to=110,by=2)] #extract timepoint values and add to new column

FCM_dat_mean$Treatment <- factor(FCM_dat_mean$Treatment, levels=c("Ambient Water Control","Heated Water Control","Non-bleached + Ambient","Non-bleached + Heated","Bleached + Ambient","Bleached + Heated")) #set treatment as a factor

FCM_dat_mean <- FCM_dat_mean[-37,] #remove NA row


# Visualize.

ggplot(FCM_dat_mean,aes(x=Timepoint,y=Mean_Concentration,fill=Treatment,shape=Treatment,color=Treatment,group=Treatment))+
  geom_point(size=7)+
  geom_line(size=2.5)+
  scale_shape_manual(values=c(21,21,21,21,21,21))+
  geom_pointrange(aes(ymin=Mean_Concentration-SE,ymax=Mean_Concentration+SE))+
  scale_color_manual(values=c("#F8766D","#00BFC4","#F8766D","#00BFC4","black","gray"))+
  scale_fill_manual(values=c("NA","NA","firebrick3","dodgerblue3","NA","NA"))+
  theme_classic()+
  theme(text=element_text(size=15),legend.key.height=unit(1.75,"cm"),complete=FALSE)+
  ylab(label="Concentration (cells per uL)")+
  xlab(label="Time (hours)")
ggsave('Microbialgrowthcurve_per_treatment.jpeg', path = dirFigs, dpi = 300)



