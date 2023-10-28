library(readr)
library(textclean)
library(expss)
library(vegan)
library(rstatix)
library(multcomp)
library(tidyverse)

mutate <- dplyr::mutate
select <- dplyr::select

# this script automatically reads and writes to folders
# the working directory is the R folder, otherwise change the code accordingly

#Which MZmine version did you use?
Mzmine<-3

#Define if you your GNPS network was "Classic" or  "FBMN" or "IIN" (ion identity networking)
networking.type <- "IIN"
#Did you run MolNetEnhancer
MolNetEnh <- "YES"
#Did you run Dereplicator separately?
Dereplicator <- "NO"
#Did you run Dereplicator+?
Dereplicator_plus <- "NO"

background_noise <- 5e4
#of how many samples does the smallest group consist?
W<-2 #this will be the minimum for transient_feature removal

#when logtransforming the data should your 0's in the dataset be replaced by something? (1000 = 3, 1 = 0, 15=1.17)
zeroes<-15


'--------------------------------------------------------------------------'

#set seed for reproducibility
set.seed(25)

# import datafiles
orbitrapsequence <- read_csv(paste0(dirRAW, "/Metabolomics/Orbitrapsequence.csv"))

  if("File Name" %in% colnames(orbitrapsequence) &
       "Sample Name" %in% colnames(orbitrapsequence) &
       "Injection_Type" %in% colnames(orbitrapsequence)){
      print("orbitrapsequence file correct")
    }else{
      stop("Wrong column names in orbitrapsequense.csv")}

rawpeakareas <- read_csv(paste0(dirRAW, "/Metabolomics/Moorea2021_IIN_quant.csv"))

# metadata <- read_csv(paste0(dirRAW, "/Metadata.csv"))
#
#     if( "Sample Name" %in% colnames(metadata) ){
#         print("metadatasheet accepted")
#       }else{
#         stop("No column with \"Sample Name\" as columnname. check metadata.csv ")}
#

# load in files from Cytoscape based on the type of networking
  dir_library_hits <- paste0(dir_analogs_on, '/DB_result')
  dir_analogs_hits <- paste0(dir_analogs_on, '/DB_analogresult')
  dir_node_info <- paste0(dir_analogs_on, '/clusterinfo_summary')

# Molnetenhancer results
  dir_library_hits <- paste0(dirCyto, '/MolNetEnhancer/DB_result')
  dir_ClassyFire_hits <-paste0(dirCyto, '/MolNetEnhancer/output_network')
  ClassyFire_hits <- read_tsv(paste0(dir_ClassyFire_hits, "/ClassyFireResults_Network.txt"))
  dir_node_info <- paste0(dirCyto, '/MolNetEnhancer/clusterinfo_summary')

# read in the files from GNPS with irregular names
setwd(dir_node_info)
temp <- list.files(pattern = "*.tsv")
node_info <- read_tsv(temp)
setwd(wd.project)

setwd(dir_library_hits)
temp <- list.files(pattern = "*.tsv")
library_hits <- read_tsv(temp)
setwd(wd.project)

setwd(dir_analogs_hits)
temp <- list.files(pattern = "*.tsv")
analogs_hits <- read_tsv(temp, col_types = "c")
setwd(wd.project)


####PREP FILES####
##### MZmine OUTPUT####

# remove the .mzXML.Peak.area string from the filenames
temp <- as.character(colnames(rawpeakareas)) #move first row to temp
tmp <-  as.data.frame(matrix(ncol = ncol(rawpeakareas), nrow = 1))#make dataframe to store colnames in
tmp[1, ] <- sub(".mzXML Peak area", "", temp)  #change temp by removing .mzXML.Peak area in insert in rawpeakareas
tmp[1, ] <- sub(".mzML Peak area", "", tmp)  #change temp by removing .mzML.Peak area in insert in rawpeakareas

colnames(rawpeakareas) <- tmp[1, ]

# if there is an extra empty column imported, that we don't need now.
# there are standard some columns imported that we don't want.
# I'm going to delete them by name, to be sure that if in the future column names change we catch it and not accidentally remove the first samples.

rawpeakareas <-
  rawpeakareas %>% select(
    - "row ion mobility",
    - "row ion mobility unit",
    - "row CCS",
    - "correlation group ID" ,
    - "annotation network number",
    - "best ion",
    - "auto MS2 verify",
    - "identified by n=",
    - "partners",
    - "neutral M mass"
  )


# if there are other columns with NA's that we also want to remove now,
# therefore we loop through the first 5 lines to find all columns with NA
# often the last column
temp <- c()
for (i in 1:(ncol(rawpeakareas))) {
  if (all(is.na(rawpeakareas[1:5, i])))
  {
    y <- i
    temp <- c(temp, y)
  }
}

rawpeakareas <- rawpeakareas %>%
  subset(select = -c(temp)) %>%
  rename("feature_nr" = "row ID")



##### PREP Cytoscape output #####

# rename the column of feature number/scan nr/clusterindex
node_info <- node_info %>%
  rename("feature_nr" = "cluster index") %>%
  rename("network" = "componentindex") %>%
  select(!contains("GNPSGROUP:")) %>%
  select(!contains("ATTRIBUTE")) %>%
  select(!num_range("G", 1:6))

library_hits <- library_hits %>%
  rename("feature_nr" = "#Scan#")

analogs_hits <- analogs_hits %>%
  rename("feature_nr" = "#Scan#")%>%
  rename("Analog_LibraryID" = "Compound_Name")

ClassyFire_hits <- ClassyFire_hits %>%
  rename("feature_nr" = "cluster index")


# remove the brackets in analog hits and library hits
analogs_hits$Analog_LibraryID <-
  textclean::replace_non_ascii(analogs_hits$Analog_LibraryID)
library_hits$Compound_Name <-
  textclean::replace_non_ascii(library_hits$Compound_Name)
node_info$LibraryID <-
  textclean::replace_non_ascii(node_info$LibraryID)



###### Make feature info file #####
feature_info1 <- as.data.frame(rawpeakareas[, 1:3])

# join library hits, analog hits and super computer predictions
feature_info2 <-
  full_join(node_info,
            full_join(
              library_hits,
              analogs_hits,
              by = "feature_nr",
              suffix = c("_Library", "_Analog")
            ),
            by = "feature_nr")

feature_info1$feature_nr <- as.numeric(feature_info1$feature_nr)
feature_info2$feature_nr <- as.numeric(feature_info2$feature_nr)
feature_info <-
  full_join(feature_info1, feature_info2, by = "feature_nr")
rm(feature_info1)
rm(feature_info2)

if (exists("ClassyFire_hits")) {
  feature_info <-
    left_join(feature_info, ClassyFire_hits, by = "feature_nr")
}

feature_info$feature_nr <- as.character(feature_info$feature_nr)


###### Make feature ID #####
# we first make a dataframe with all the columns that we want to use.

df.featureID <- feature_info %>%
  select(1:3, network, LibraryID, Analog_LibraryID, CF_kingdom, CF_superclass, CF_class, CF_subclass)

featureID <- df.featureID %>%
  mutate(across(`row m/z`:`row retention time`, as.numeric)) %>%
  mutate(`row m/z` = round(`row m/z`, 4)) %>%
  mutate(`row retention time` = round(`row retention time`, 4)) %>%
  unite(featureID, everything())


# in between we will save all the info we have on the features
featureID_info <- bind_cols(featureID, df.featureID)
df.featureID <- featureID_info
featureID_info <- left_join(featureID_info, feature_info, by = "feature_nr")
# setwd(wd.project)
# setwd(dirOutput)
# write.csv(featureID_info, "feature_info_final.csv", row.names = FALSE) #write all feature info available
# setwd(wd.project)

# make a new df with combined names
df <- cbind(featureID, rawpeakareas[, 4:ncol(rawpeakareas)])
df <- t(df)
df <- as.data.frame(df)

# we name the columns and the rows by naming columns and transposing
colnames(df) <- as.character(unlist(df[1, ]))
df <- df[-1, ]

# combine rawpeakareas with orbitrap sequence
shared.name <- as.data.frame(rownames(df))
colnames(shared.name) <- 'File Name'
df <- bind_cols(shared.name, df)
df1 <- right_join(orbitrapsequence, df, by = "File Name") #join matching info from orbitrap sequence

# setwd(wd.project)
# setwd(dirOutput)
# write.csv(df1, "rawpeakareas.csv", row.names = TRUE) #write third version of data, not cleaned
# setwd(wd.project)

rm(rawpeakareas, library_hits, analogs_hits, node_info, featureID_info, shared.name, df, feature_info, featureID, tmp, tmp2)

# blank correction:
# we are going to flag and remove features if max(blanks) >= mean(peak area)*0.5, thus over all samples from experiment (not the 800 samples)
#devide the dataset in blanks and samples
InjectionType<-as.factor(df1$Injection_Type)
blanks<-subset(df1, df1$Injection_Type!= "Sample")
samples<-subset(df1, df1$Injection_Type == "Sample")
blanks<-as.data.frame(t(blanks))
samples<-as.data.frame(t(samples))


#this part sucks because of all the steps before there are lists loaded etc and the numbers are read as characters or factors
#in order to fix it I have to save the file and re-load the files, in this way they are imported as numerical.
#there is probably an easier way with dplyr but I havent coded that yet since this works
blanks.name<-as.data.frame(unlist(blanks[1,]))
blanks.name<-blanks.name$`unlist(blanks[1, ])`

samples.name<-as.data.frame(unlist(samples[1,]))
samples.name<-samples.name$`unlist(samples[1, ])`

colnames(blanks)<- blanks.name
colnames(samples)<-samples.name

blanks<-dplyr::slice(blanks, 4:nrow(blanks))
samples<-dplyr::slice(samples, 4:nrow(samples))
setwd(wd.project)
setwd(dirWrite)
write.csv(blanks,"blanks.csv",row.names = TRUE)
write.csv(samples,"samples.csv",row.names = TRUE)

rm(blanks)
rm(samples)

blanks<-read.csv("blanks.csv", sep=",", header = TRUE, row.names = 1, check.names = TRUE, stringsAsFactors=FALSE)
samples<-read.csv("samples.csv", sep=",", header = TRUE, row.names = 1, check.names = TRUE, stringsAsFactors=FALSE)

setwd(wd.project)

## flagging background features and subtraction from samples ------------------------------------------------
# The idea here is to flag and remove features where max(blanks) >= mean(samples)*0.5
# so flagged are removed, non flagged you keep.
# things are flagged with 1, not flagged = 0, this is the opposite than flagging_per_sample
blanks$max_blanks <- apply(blanks, 1, max)
samples$mean_samples <- apply (samples, 1, mean)

max_blanks<-dplyr::select(blanks, select = "max_blanks") %>% dplyr::rename("max_blanks" = "select")
mean_samples<-dplyr::select(samples, select = "mean_samples") %>% dplyr::rename("mean_samples" = "select")
flagging<-cbind(max_blanks, mean_samples)

rm(max_blanks)
rm(mean_samples)

#flag if blanks is more than half of the mean in the sample
flagging$flag<-if_else((flagging$mean_samples*0.5) < flagging$max_blanks, 1, 0, missing = NULL)
flagging<-rownames_to_column(flagging, var = "feature")

#split in flagged and not flagged
flagged<-dplyr::filter(flagging, flagging$flag == 1)
notflagged <-dplyr::filter(flagging, flagging$flag == 0)
flagged<-tibble::column_to_rownames(flagged, 'feature')

#keep the non flagged features
notflagged<-(as.character(notflagged$feature))
df1.filtered<- df1 %>% dplyr::select("File Name", one_of(notflagged))

df1.filtered<-tibble::column_to_rownames(df1.filtered, 'File Name')
df1.filtered<-as.data.frame(t(df1.filtered))

sample_files<-dplyr::filter(orbitrapsequence, Injection_Type == "Sample")
sample_files<-sample_files$`File Name`

#select and keep the columns that are samples
df1.filtered_noblanks<-df1.filtered %>% select(all_of(sample_files))


df.trans<-df1.filtered_noblanks

#create an empty matrix to fill in for every feature in how many samples (count) the area under the peak is higher than the set background noise.
df.count<-as.data.frame(matrix(ncol=1, nrow = nrow(df.trans)))
for (i in 1:dim(df.trans)[1]){
  df.count[i,1] <-sum(df.trans[i,]>background_noise)}

#Give th counts row names so we can add that column to the original table
rownames(df.count)<-rownames(df.trans)
df.filtered<-bind_cols(df.trans, df.count)
df.filtered<-rownames_to_column(df.filtered, 'feature')

#we now split the file into the transient features (occuring less than W times) in df.trans (W is the size of the smallest group)
#and the ones that are above the backgroud noise level more than W times (filtered)
df.trans<-dplyr::filter(df.filtered, df.filtered$V1 < W)
df.trans<-dplyr::select(df.trans, -V1)
df.filtered<-dplyr::filter(df.filtered, df.filtered$V1 >= W)
df.filtered<-dplyr::select(df.filtered, -V1)

df.trans<-column_to_rownames(df.trans, 'feature')
df.filtered<-column_to_rownames(df.filtered, 'feature')

write.csv(df.trans,"transient.feat_filtered_out.csv",row.names = TRUE)
write.csv(df.filtered,"rawpeaks_no-background_no-transientfeat.csv",row.names = TRUE)

#remove some objects we don't need anymore
rm(df.trans, df.count)
rm(flagging, blanks, samples)

#combine filtered data with metadata to create df.area but also save one without metadata but with samples in rows.
df.area<-as.data.frame(t(df.filtered))
df.area <- df.area %>% mutate_if(is.character, as.numeric)
df.area<-rownames_to_column(df.area, var = "File Name")
df.area.no.meta<-df.area

#create metadata and add.
if (!exists("full_metadata")){
  full_metadata<-dplyr::right_join(orbitrapsequence, metadata, by = "Sample Name")}
df.area<-dplyr::right_join(full_metadata, df.area, by = "File Name")
df.area<-if_na(df.area, "not applicable")

setwd(wd.project)


### normalize by TIC
df1.trans.feat<-df.filtered #create backup of your data, although this should be the same as "rawpeaks_no-background_no-transientfeat.csv"
df.filtered<-as.data.frame(t(df.filtered))

df.filtered<-as.data.frame(sapply(df.filtered,as.numeric))
#Sum all raw area's under the peak to get a Total Ion Current
df.filtered$TIC<-apply(df.filtered, 1 ,sum)
rownames(df.filtered)<-colnames(df1.trans.feat)

#Normalize by TIC (area under the peak/TIC)
#define function
FUNTIC<-function (x) {
  x / df.filtered$TIC
}

#create matrix to put data in
TICNORM.smpl<-data.frame(matrix(ncol=ncol(df.filtered), nrow = nrow(df.filtered)))
#apply function
TICNORM.smpl<-sapply(df.filtered, FUNTIC)
#change to data frame and add rownames etc.
TICNORM.smpl<-as.data.frame(TICNORM.smpl)
rownames(TICNORM.smpl)<-rownames(df.filtered)
TICNORM.smpl<-TICNORM.smpl[,1:(ncol(TICNORM.smpl)-1)]

#create function to transform data by asin of sqare root of the normalized data.
FUNSQRT<-function(x){
  asin(sqrt(x))
}

#create matrix to put data in
ASIN_sqrt_smpl<-data.frame(matrix(ncol=ncol(TICNORM.smpl), nrow = nrow(TICNORM.smpl)))
#apply function
ASIN_sqrt_smpl<-sapply(TICNORM.smpl, FUNSQRT)
#change to data frame and add rownames etc.
ASIN_sqrt_smpl<-as.data.frame(ASIN_sqrt_smpl)
rownames(ASIN_sqrt_smpl)<-rownames(TICNORM.smpl)

df.norm.smpl<-ASIN_sqrt_smpl
df.norm.smpl<-rownames_to_column(df.norm.smpl, "File Name")
df.norm.smpl<-dplyr::right_join(full_metadata, df.norm.smpl, by = "File Name")
df.norm.smpl<-if_na(df.norm.smpl, "not applicable")

#number of metadata columns
M<-ncol(full_metadata)+1


#if you normalized by sample you have a backup of your filtered data as df.1.trans.feat (and your df.filtered will be changed!))
df2.filtered<-df1.trans.feat

#Sum all raw area's under the peak to get a Total Ion Current
df2.filtered <- df2.filtered %>% mutate_if(is.character, as.numeric)
df2.filtered$TIC<-apply(df2.filtered, 1 ,sum)

#Normalize by TIC (area under the peak/TIC)
#define function
FUNTIC<-function (x) {
  x / df2.filtered$TIC
}

#create matrix to put data in
TICNORM.ft<-data.frame(matrix(ncol=ncol(df2.filtered), nrow = nrow(df2.filtered)))
#apply function
TICNORM.ft<-sapply(df2.filtered, FUNTIC)
#change to data frame and add rownames etc.
TICNORM.ft<-as.data.frame(TICNORM.ft)
rownames(TICNORM.ft)<-rownames(df2.filtered)
TICNORM.ft<-TICNORM.ft[,1:(ncol(TICNORM.ft)-1)]

#create function to transform data by asin of sqare root of the normalized data.
FUNSQRT<-function(x){
  asin(sqrt(x))
}

#create matrix to put data in
ASIN_sqrt.ft<-data.frame(matrix(ncol=ncol(TICNORM.ft), nrow = nrow(TICNORM.ft)))
#apply function
ASIN_sqrt.ft<-sapply(TICNORM.ft, FUNSQRT)
#change to data frame and add rownames etc.
ASIN_sqrt.ft<-as.data.frame(ASIN_sqrt.ft)
rownames(ASIN_sqrt.ft)<-rownames(TICNORM.ft)

df.norm.ft<-as.data.frame(t(ASIN_sqrt.ft))
df.norm.ft<-rownames_to_column(df.norm.ft, "File Name")

#join with full metadata
df.norm.ft<-dplyr::right_join(full_metadata, df.norm.ft, by = "File Name")
df.norm.ft<-if_na(df.norm.ft, "not applicable")

rm(df2.filtered)

# #create backup of your data, although this should be the same as "rawpeaks_no-background_no-transientfeat.csv"
#
#   df.filtered<-df1.trans.feat
#   df.filtered<-as.data.frame(t(df.filtered))
#   df.filtered<-as.data.frame(sapply(df.filtered,as.numeric))
#   rownames(df.filtered)<-colnames(df1.trans.feat)
#
#


#Normalize by log 10 transformation
#define function
FUNLOG<-function (x) {
  log10(x)
}

#create matrix to put data in
LOGNORM<-data.frame(matrix(ncol=ncol(df.filtered), nrow = nrow(df.filtered)))
#apply function
LOGNORM<-sapply(df.filtered, FUNLOG)
#change to data frame and add rownames etc.
LOGNORM<-as.data.frame(LOGNORM)
rownames(LOGNORM)<-rownames(df.filtered)

#replace the 0s with the number you want. - i couldve done this before log transformation. but now i also have a raw version)
tmp<-LOGNORM
set.seed(25)
log10zeroes<-runif(sum(LOGNORM == -Inf))
log10zeroes<-log10zeroes + zeroes
log10zeroes<-log10(log10zeroes)
tmp[which(LOGNORM == -Inf, arr.ind = TRUE)]<- log10zeroes

df.norm.area<-tmp
df.norm.area<-rownames_to_column(df.norm.area, "File Name")
rm(log10zeroes)

#now create metadata and add to the normalized and transformed data
df.norm.area<-dplyr::right_join(full_metadata, df.norm.area, by = "File Name")
df.norm.area<-if_na(df.norm.area, "not applicable")

# create z-score dataset
tmp <- df.norm.area[M:ncol(df.norm.area)]
tmp <- tmp %>% mutate(across(where(is.numeric), ~ as.numeric(scale(.))))
df.zscore <- bind_cols(df.norm.area[1:(M-1)], tmp)


setwd(wd.project)
#we can remove these for now
rm(ASIN_sqrt_smpl, ASIN_sqrt.ft, df.filtered, df1, df1.filtered, df1.filtered_noblanks, df1.trans.feat, flagged, metadata, tmp, ClassyFire_hits)



# #non metric multidimentional scaling, ABCT0Dom T0
#ordination of data

df.area.ABCT0 <- df.area %>% filter(Timepoint_char == "T0") %>%
  filter(Treatment != "Inoculum")

#ordination df.area
ord.mod.area<- metaMDS(df.area.ABCT0[M:ncol(df.area.ABCT0)], distance = 'bray', k = 2)
NMDS.ABCDom.T0<- bind_cols(df.area.ABCT0[1:(M-1)], as.data.frame(ord.mod.area$points))
#extract stress
stress <- round(ord.mod.area$stress, 4)
stress <- paste0("ABCDom T0 bray curtis dissimilarity k=2 \nstress = ", stress)

cost.col.fill<-c("dodgerblue1","firebrick1", "white", "white", "grey70", "grey70")
cost.col.line<-c("dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3")
fact.all.treat<-factor(NMDS.ABCDom.T0$Treatment, levels = c("Control", "Heated", "Bleached", 'Bleached + Heated', "Negative Control", "Negative Control + Heated"))
treat.labels <- c("Control", "Heated", "Bleached", 'Bleached + Heated', "Negative Control", "Negative Control + Heated")

NMDS.ABCDom.T0 <- NMDS.ABCDom.T0 %>%
  mutate(stress_status_v1 = case_when(Treatment == "Control" ~ "Control",
                                      Treatment == "Negative Control" ~ NA,
                                      Treatment == "Negative Control + Heated" ~ NA,
                                      .default = "Stressed"))

fig4A<-ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data=NMDS.ABCDom.T0,
             aes(x=MDS1, y=MDS2, fill = fact.all.treat, colour=fact.all.treat),
             size=5, stroke=1.5, shape = 21) +
  scale_color_manual(labels = treat.labels, values = cost.col.line, name = "Treatment")+
  scale_fill_manual(labels = treat.labels, values = cost.col.fill, name = "Treatment", guide = guide_legend(override.aes = list(shape = 21)))+
  annotate("text", label="p > 0.001 \n stress = 0.0467", x=-.11, y=.05)+
  theme_bw()+
  theme(axis.ticks.x.top = element_blank(), axis.ticks.y.right = element_blank())+
  coord_fixed(ratio=1.2)+
  ggforce::geom_mark_ellipse(data=NMDS.ABCDom.T0, aes(x=MDS1, y=MDS2, linetype=stress_status_v1, label=stress_status_v1), con.type="none", label.buffer=unit(4,'mm'), show.legend=F)
fig4A
ggsave("4A_NMDS_ABCDom_T0.jpg", path = dirFigs, width = 6.75, height = 5, units = "in", dpi = 320)

write_csv(NMDS.ABCDom.T0, file = file.path(dirOutput, "fig4A_data.csv"))

#PERMANOVA // ADONIS
# all treatments incl. water controls
adonis2(df.area.ABCT0[,M:ncol(df.area.ABCT0)] ~ fact.all.treat, method="bray",perm=999) %>%
  print()
# in paper

tmp <- df.area.ABCT0  %>%
  mutate(stress_status_v1 = case_when(Treatment == "Negative Control" ~ "Water",
                                      Treatment == "Negative Control + Heated" ~ "Water",
                                      .default = "Coral"))
factor(tmp$stress_status_v1)

# water controls vs. corals
adonis2(tmp[(M):(ncol(tmp)-1)] ~ factor(tmp$stress_status_v1), method="bray",perm=999)%>%
  print()


# watercontrols  vs. control(ambient) vs stressed

tmp <- df.area.ABCT0  %>%
  mutate(stress_status_v1 = case_when(
    Treatment == "Control" ~ "Ambient",
    Treatment == "Negative Control" ~ "Control",
    Treatment == "Negative Control + Heated" ~ "Control",
    .default = "Stressed"))
factor(tmp$stress_status_v1)

# watercontrols  vs. control(ambient) vs stressed
adonis2(tmp[(M):(ncol(tmp)-1)] ~ factor(tmp$stress_status_v1), method="bray",perm=999)%>%
  print()


# only corals
tmp <- df.area.ABCT0  %>%
  mutate(stress_status_v1 = case_when(Treatment == "Control" ~ "Ambient",
                                      Treatment == "Negative Control" ~ "Control",
                                      Treatment == "Negative Control + Heated" ~ "Control",
                                      .default = "Stressed")) %>%
  filter(Bleaching_Status != "not applicable")

# control vs. heated vs. bleached vs. bleached + heated
set.seed(86)
tmp2<-adonis2(tmp[(M):(ncol(tmp)-1)] ~ factor(tmp$Treatment), method="bray",perm=999)%>%
  print()
# in paper

# control vs. stressed
adonis2(tmp[(M):(ncol(tmp)-1)] ~ factor(tmp$stress_status_v1), method="bray",perm=999)%>%
  print()



