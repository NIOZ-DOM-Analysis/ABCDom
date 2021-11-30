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
FCM_dat <- merge(FCM.dat, metadata.fcm, by.x="SampleID", by.y="Sample.Name", all.x=T, all.y=F)

write.csv(FCM_dat, file.path(dirOutput, "FCM_dat.csv"), )


