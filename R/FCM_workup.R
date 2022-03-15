'FCM_workup.R'

PlateII.Map <- read.csv(file.path(dirRAW, "FCM", "abcDOM.SYBR.NOV.15.Plate.II.Map.csv"))
PlateI.Map <- read.csv(file.path(dirRAW, "FCM", "abcDOM.SYBR.NOV.15.Plate.I.Map.csv"))
PlateII <- read.csv(file.path(dirRAW, "FCM", "abcDOM.SYBR.NOV.15.Plate.II.csv"))
PlateI <- read.csv(file.path(dirRAW, "FCM", "abcDOM.SYBR.NOV.15.Plate.I.csv"))

#subset for only SYBR Polygon
PlateI <- PlateI[PlateI$Gate.Name=="SYBR Polygon",c(4,21)]
PlateII <- PlateII[PlateII$Gate.Name=="SYBR Polygon",c(4,21)]

#merge plate map and data
FCM.I <- merge(PlateI,PlateI.Map,by.x="Sample.Name",by.y="Well")
FCM.II <- merge(PlateII,PlateII.Map,by.x="Sample.Name",by.y="Well")
FCM.dat <- rbind(FCM.I,FCM.II) #combine the merged dfs

FCM.dat$Concentration <- as.numeric(as.character(sub(",","",FCM.dat$Concentration))) #commas were causing a probem when converting to numeric, have to convert to character first and sub out commas and then convert to numeric.

#merge FCM.dat with metadata
metadata.fcm <- subset(metadata, Sample_Type == "FCM")
FCM_dat <- merge(FCM.dat, metadata.fcm, by.x="SampleID", by.y="Sample Name", all.x=T, all.y=F)

#write.csv(FCM_dat, file.path(dirOutput, "FCM_dat.csv"), )

#calculate specific growth rate between t0 and t24.
#wor up the data
FCM_dat_growth <- FCM_dat[FCM_dat$Timepoint..h.==0,] #subset FCM_dat for just t0 samples, save in new df
colnames(FCM_dat_growth)[3] <- "T0_Concentration" #rename concentration column
FCM_dat_growth$T24_Concentration <- FCM_dat[FCM_dat$Timepoint..h.==24,3] #extract concentration values at t24 and add to new column

#calculate specific growth rate
FCM_dat_growth$log10_T0_Concentration <- log10(FCM_dat_growth$T0_Concentration) #convert to log10 concentrations
FCM_dat_growth$log10_T24_Concentration <- log10(FCM_dat_growth$T24_Concentration) #convert to log10 concentrations
FCM_dat_growth$Specific_Growth_Rate <- ((FCM_dat_growth$log10_T24_Concentration - FCM_dat_growth$log10_T0_Concentration)*2.303)/24 #calculate converted growth rate.
FCM_dat_growth <- FCM_dat_growth[is.na(FCM_dat_growth$SampleID)!=TRUE,] #remove NAs
FCM_dat_growth$Treatment <- factor(FCM_dat_growth$Treatment.y, levels=levels(fact.all.treat)) #reset factor levels

rm(PlateI, PlateII,PlateI.Map,PlateII.Map, FCM.dat,FCM.I,FCM.II)
