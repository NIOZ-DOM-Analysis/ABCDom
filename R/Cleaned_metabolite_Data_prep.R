'Cleaned_metabolite_Data_prep.R'
setwd(wd.project)

#import file
peakareas<-read_csv(paste0(dirRAW, '/peakareas_blankremoved.csv'))
orbitrapsequence <- read_csv(paste0(dirRAW,"/orbitrapsequence.csv"))
metadata <- read_csv(paste0(dirRAW,"/metadata.csv"))
# SIRIUS <-read_table(paste0(dirCyto, '/MzMine output Blank removed/Moorea2021_IIN_SIRIUS_BlanksRemoved.mgf.txt'))

# remove the .mzXML.Peak.area string from the filenames
temp <- as.character(colnames(peakareas)) #move first row to temp
tmp <-  as.data.frame(matrix(ncol = ncol(peakareas), nrow = 1))#make dataframe to store colnames in
tmp[1, ] <- sub(".mzXML Peak area", "", temp)  #change temp by removing .mzXML.Peak area in insert in rawpeakareas
tmp[1, ] <- sub(".mzML Peak area", "", tmp)  #change temp by removing .mzML.Peak area in insert in rawpeakareas

colnames(peakareas) <- tmp[1, ]

# if there is an extra empty column imported, that we don't need now.
# there are standard 12 columns imported that we don't want.
# I'm going to delete them by name, to be sure that if in the future column names change we catch it and not accedentally remove the first samples.

peakareas <-
  peakareas %>% select(
    # -"row identity (main ID)",
    # - "row identity (all IDs)",
    # - "row identity (main ID + details)",
    # - "row comment",
    # - "row number of detected peaks",
    - "correlation group ID" ,
    - "annotation network number",
    - "best ion",
    - "auto MS2 verify",
    - "identified by n=",
    - "partners",
    - "neutral M mass"
  )


# if there are other columns with NA's that we also want to remove now, therefore we loop through the first 5 lines to find all columns with NA
temp <- c()

for (i in 1:(ncol(peakareas))) {
  if (all(is.na(peakareas[1:5, i])))
  {
    y <- i
    temp <- c(temp, y)
  }
}
peakareas <- subset(peakareas, select = -c(temp))

# rename the naming of feature number
peakareas <- dplyr::rename(peakareas, "feature_nr" = "row ID")

analysis_info<-list()
analysis_info$total_nr_features <- nrow(peakareas)
analysis_info$total_nr_runs <- ncol(peakareas)-3

#load cytoscape data
temp <- list.files(path = paste0(dirCyto, '/MzMine output Blank removed/clusterinfo_summary'), pattern = "*.tsv")
node_info <- read_tsv(paste0(dirCyto, '/MzMine output Blank removed/clusterinfo_summary/', temp))
temp <- list.files(path = paste0(dirCyto, '/MzMine output Blank removed/DB_result'), pattern = "*.tsv")
library_hits <- read_tsv(paste0(dirCyto, '/MzMine output Blank removed/DB_result/', temp))
temp <- list.files(path = paste0(dirCyto, '/MzMine output Blank removed/DB_analogresult'), pattern = "*.tsv")
analogs_hits <- read_tsv(paste0(dirCyto, '/MzMine output Blank removed/DB_analogresult/', temp))

# rename the column of feature number/scan nr/clusterindex
tmp <-
  which(colnames(node_info) == "cluster index") #lookup which colnumber is clusterindex to rename
colnames(node_info)[tmp] <-
  "feature_nr" #rename cluster.index to feature_nr

tmp <-
  which(colnames(node_info) == "componentindex") #lookup which colnumber is componentindex to rename
colnames(node_info)[tmp] <-
  "network" #rename componentindex to network

node_info <- dplyr::select(node_info,!contains("GNPSGROUP:"))
node_info <- dplyr::select(node_info,!contains("ATTRIBUTE"))
node_info <- dplyr::select(node_info,!num_range("G", 1:6))

tmp <- which(colnames(library_hits) == "#Scan#")
colnames(library_hits)[tmp] <- "feature_nr"

tmp <- which(colnames(analogs_hits) == "#Scan#")
colnames(analogs_hits)[tmp] <- "feature_nr"

tmp <- which(colnames(analogs_hits) == "Compound_Name")
colnames(analogs_hits)[tmp] <- "Analog_LibraryID"

# remove the brackets in analog hits and library hits
library_hits$Compound_Name <-
  textclean::replace_non_ascii(library_hits$Compound_Name)
analogs_hits$Analog_LibraryID <-
  textclean::replace_non_ascii(analogs_hits$Analog_LibraryID)
node_info$LibraryID <-
  textclean::replace_non_ascii(node_info$LibraryID)

analysis_info$nr_analog_hits<-nrow(analogs_hits)
analysis_info$nr_library_hits<-nrow(library_hits)

###### Make feature info file #####
feature_info1 <- as.data.frame(peakareas[, 1:3])

# join library hits
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
  join(feature_info1, feature_info2, by = "feature_nr")
rm(feature_info1)
rm(feature_info2)

feature_info$feature_nr <- as.character(feature_info$feature_nr)

df.featureID <- feature_info[, 1:3]
df.featureID <-
  cbind(
    df.featureID,
    feature_info$network,
    feature_info$LibraryID,
    feature_info$Analog_LibraryID
  )

# remove all feature_info$ parts of the column names
temp <- as.character(colnames(df.featureID)) #move first row to temp
tmp <-
  as.data.frame(matrix(ncol = ncol(df.featureID), nrow = 1))#make dataframe to store colnames in
temp <- gsub("[^0-9A-Za-z///' ]", "'" , temp)
temp <- sub("feature'info'", "", temp)
tmp[1, ] <- gsub("'", "_", temp)
colnames(df.featureID) <- tmp[1, ]

# concatenate the feature names from the different files (row.ID, row.m.z, row.retention.time, Library ID and Analog ID and sirius_formula when present, Classifire classes)
featureID <- c(df.featureID, sep = "_") #get the columns to combine
featureID$`row m/z` <- as.numeric(featureID$`row m/z`, digits = 4)
featureID$`row m/z` <-
  round(featureID$`row m/z`, digits = 4) #round digits of m/z
featureID$`row retention time` <-
  as.numeric(featureID$`row retention time`)
featureID$`row retention time` <-
  round(featureID$`row retention time`, digits = 4) #round digits of retention time
featureID <- do.call(paste, featureID) #concatenate
featureID <- as.data.frame(featureID) #make data frame

# in between we will save all the info we have on the features
featureID_info <- cbind(featureID, df.featureID)
df.featureID <- featureID_info
featureID_info <-
  join(featureID_info, feature_info, by = "feature_nr", type = "left")
write.csv(featureID_info, paste0(dirOutput, "/feature_info_final.csv"), row.names = FALSE) #write all feature info availible


# make a new df with combined names
df <- cbind(featureID, peakareas[, 4:ncol(peakareas)])
df <- t(df)
df <- as.data.frame(df)
# we name the columns and the rows by naming columns and transposing
colnames(df) <- as.character(unlist(df[1, ]))
df <- df[-1, ]



# combine rawpeakareas with orbitrap sequence
shared.name <- as.data.frame(rownames(df))
colnames(shared.name) <- 'File Name'
df <- cbind(shared.name, df)
df1 <-
  right_join(orbitrapsequence, df, by = "File Name") #join matching info from orbitrap sequence
write.csv(df1, paste0(dirOutput, "/peakareas.csv"), row.names = TRUE) #write third version of data, not cleaned

analysis_info$nr_selected_runs <- sum(!is.na(df1$Injection_Type))
analysis_info$nr_not_selected_runs <- sum(is.na(df1$Injection_Type))

#select only the samples that have a sample name
df1<- df1 %>% filter(!is.na(Injection_Type))

#also kick out all the columns that are 0
#first make everything numeric
sapply(df1, class)
df1[4:ncol(df1)] <- sapply(df1[4:ncol(df1)], as.numeric)
sapply(df1, class)

temp<-apply(df1[4:ncol(df1)], 2 ,sum)
temp<-as.data.frame(temp)
temp <- temp %>% filter(temp == 0)
temp <- rownames(temp)

analysis_info$empty_features_removed <- length(temp)

#filter out the columns that had a sum of 0
df1 <- df1 %>% select(!all_of(temp))


#####
# add another transient feature removal?
#create an empty matrix to fill in for every feature in how many samples (count) the area under the peak is higher than the set background noise.
df.trans <- as.data.frame(t(df1))
colnames(df.trans)<-df.trans[1,]
df.trans <- df.trans[-1,]
df.trans <- df.trans[-1,]
df.trans <- df.trans[-1,]

background_noise <- 5E4
W <- 3

df.count<-as.data.frame(matrix(ncol=1, nrow = nrow(df.trans)))
for (i in 1:dim(df.trans)[1]){
  df.count[i,1] <-sum(df.trans[i,]>background_noise)}

rownames(df.count)<-rownames(df.trans)
df.trans<-cbind(df.trans, df.count)
df.trans<-rownames_to_column(df.trans, 'feature')
df.filtered <- df.trans

df.trans<-dplyr::filter(df.trans, df.trans$V1 < W)
df.trans<-dplyr::select(df.trans, -V1)
df.filtered<-dplyr::filter(df.filtered, df.filtered$V1 >= W)
df.filtered<-dplyr::select(df.filtered, -V1)

df.trans<-column_to_rownames(df.trans, 'feature')
df.filtered<-column_to_rownames(df.filtered, 'feature')


write.csv(df.trans,"transient.feat_filtered_out.csv",row.names = TRUE)
write.csv(df.filtered,"rawpeaks_no-background_no-transientfeat.csv",row.names = TRUE)

analysis_info$nr_transient_features_removed<-nrow(df.trans)
analysis_info$nr_features_for_analysis<-nrow(df.filtered)

rm(df.trans, df.count)

####

df1<-as.data.frame(t(df.filtered))
df1<-rownames_to_column(df1, "File Name")
df1 <-right_join(orbitrapsequence, df1, by = "File Name") #join matching info from orbitrap sequence

#####

df1[4:ncol(df1)]<- lapply(df1[4:ncol(df1)], function(x) as.numeric(as.character(x)))
#now lets calculate the TIC
df1$TIC<-apply(df1[4:ncol(df1)], 1 ,sum)

#Normalize by TIC (area under the peak/TIC)
#define function
FUNTIC<-function (x) {
  x / df1$TIC
}


#create matrix to put data in
TICNORM.smpl<-data.frame(matrix(ncol=ncol(df1[4:(ncol(df1)-1)]), nrow = nrow(df1)))
#apply function
TICNORM.smpl<-sapply(df1[4:(ncol(df1)-1)], FUNTIC)
#change to data frame and add rownames etc.
TICNORM.smpl<-as.data.frame(TICNORM.smpl)
TICNORM.smpl<-cbind(df1[1:3],TICNORM.smpl)

write.csv(TICNORM.smpl,"TICNORM.smpl.csv",row.names = TRUE)

#create function to transform data by asin of sqare root of the normalized data.
FUNSQRT<-function(x){
  asin(sqrt(x))
}

#create matrix to put data in
ASIN_sqrt_smpl<-data.frame(matrix(ncol=ncol(df1[4:(ncol(df1)-1)]), nrow = nrow(df1)))
#apply function
ASIN_sqrt_smpl<-sapply(TICNORM.smpl[4:(ncol(TICNORM.smpl))], FUNSQRT)
#change to data frame
ASIN_sqrt_smpl<-as.data.frame(ASIN_sqrt_smpl)
ASIN_sqrt_smpl<-cbind(df1[1:3],ASIN_sqrt_smpl)

df.norm.smpl<-ASIN_sqrt_smpl

write.csv(df.norm.smpl,"df.norm.smpl_no_metadata.csv",row.names = FALSE)

#should your 0's in the dataset be replaced by something? (1000 = 3, 1 = 0)
zeroes<-1000

#or transform XIC by log 10
# here we use df1 but first delete the TIC column by saving it somewhere else
df.TIC<-df1 %>% select(1:3, TIC)
df1 <- df1 %>% select(-TIC)

#Normalize by log 10 transformation
#define function
FUNLOG<-function (x) {
  log10(x)
}

#create matrix to put data in
LOGNORM<-data.frame(matrix(ncol=ncol(df1[, 4:ncol(df1)]), nrow = nrow(df1)))
#apply function
LOGNORM<-sapply(df1[, 4:ncol(df1)], FUNLOG)
#change to data frame and add rownames etc.
LOGNORM<-as.data.frame(LOGNORM)
LOGNORM<-cbind(df1[,1:3], LOGNORM)

write.csv(LOGNORM, paste0(dirOutput, "/LOGNORM_raw.csv"),row.names = TRUE)


#replace the 0s with the number you want. - i couldve done this before log transformation. but now i also have a raw version)
tmp<-LOGNORM
set.seed(25)
log10zeroes<-runif(sum(LOGNORM == -Inf))
log10zeroes<-log10zeroes + zeroes
log10zeroes<-log10(log10zeroes)
tmp[which(LOGNORM == -Inf, arr.ind = TRUE)]<- log10zeroes

df.norm.area<-tmp
write.csv(df.norm.area, paste0(dirOutput, "/df.norm.area_no_metadata.csv"),row.names = FALSE)

analysis_info$nr_replaced_zeroes_lognorm<-length(log10zeroes)
rm(log10zeroes)

full_metadata<-dplyr::left_join(orbitrapsequence, metadata, by = "Sample Name")
df.norm.area<-df.norm.area %>% select (-'Sample Name', -Injection_Type)
df.norm.area<-dplyr::right_join(full_metadata, df.norm.area, by = "File Name")
# df.norm.area<-if_na(df.norm.area, "not applicable")

write.csv(df.norm.area,paste0(dirOutput, "/df.norm.area.csv"),row.names = FALSE)

#number of metadata columns
M<-ncol(full_metadata)+1


#now lets also add the metadata to the other datasets
df.norm.smpl<-df.norm.smpl %>% select (-'Sample Name', -Injection_Type)
df.norm.smpl<-dplyr::right_join(full_metadata, df.norm.smpl, by = "File Name")

df.TIC<-df.TIC %>% select (-'Sample Name', -Injection_Type)
df.TIC<-dplyr::right_join(full_metadata, df.TIC, by = "File Name")

TICNORM.smpl<-TICNORM.smpl %>% select (-'Sample Name', -Injection_Type)
TICNORM.smpl<-dplyr::right_join(full_metadata, TICNORM.smpl, by = "File Name")

df1<-df1 %>% select (-'Sample Name', -Injection_Type)
df.area<-dplyr::right_join(full_metadata, df1, by = "File Name")

rm(df, featureID, feature_info, df.featureID)
rm(ASIN_sqrt_smpl, df1)
# setwd(wd.project)
# getwd()
