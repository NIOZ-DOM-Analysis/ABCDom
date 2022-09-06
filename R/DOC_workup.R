'DOC_workup.R'

DOC_raw <- read.csv(file.path(dirRAW, "DOC", "DOC_with_SampleName.csv"))
if(!exists("metadata") ){
  metadata <- read.csv(file.path(dirRAW, "Metadata.csv"))}

dat <- DOC_raw[c(-1,-11,-14,-17,-19,-23,-34),]

# subset for only DOC metadata
# metadata.doc <- metadata[metadata$Sample_Type=="DOC",]

# merge with DOC dat
DOC_dat <- merge(dat, metadata, by.x="Sample.Name", by.y="Sample Name", all.x=T, all.y=F)


# Calculate the background water control corrected DOC for the t0 samples.

# calcualte background water control DOC values
ambient_background = mean(DOC_dat$uMC[DOC_dat$Treatment=="Ambient Water Control" & DOC_dat$Timepoint==0]) #calcualte ambient background
heated_background = mean(DOC_dat$uMC[DOC_dat$Treatment=="Heated Water Control" & DOC_dat$Timepoint==0]) #calcualte ambient background

# Add a new column corresponding to DOC background values for every sample.
DOC_dat$Background_DOC = DOC_dat$uMC #first copy DOC values
DOC_dat$Background_DOC[DOC_dat$Stress_status == "Thermal Stress"] = heated_background #in the background DOC column, for any sample which was thermally stressed add corresponding heated background value.
DOC_dat$Background_DOC[DOC_dat$Stress_status == "Ambient"] = ambient_background #in the background DOC column, for any sample which was ambient add corresponding ambient background value.

# Calculate control corrected DOC by subtracting the background value from the DOC value.
DOC_dat$Control_Corrected_DOC = DOC_dat$uMC - DOC_dat$Background_DOC

#Calculate DOC flux for control corrected DOC values by dividing by 3.5
DOC_dat$Control_Corrected_DOC_flux <- DOC_dat$Control_Corrected_DOC/3.5

# export to output folder
#write.csv(DOC_dat, file.path(dirOutput, "DOC_dat.csv"), )
rm(dat)

