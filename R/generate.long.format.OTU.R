#This function will take a relative abundance/raw abundance OTU feature table, a metadata table, and a taxonomy table, merge them and convert them to "longformat" for easy plotting in ggplot2.
#Note: This function currently references OTUs and taxonomy, but with some reworking it can be applied to metabolites or other feature tables. Need to work on this further.

generate.long.format.OTUs = function(abund,metadata,taxonomy) {
  #Inputs:  #abund = relative or raw abundance dataframe, with the samples as columns and OTUs as rows. Column names must correspond exactly to sample names in metadata and row names must correspond exactly to OTU names in taxonomy.
            
            #metadata = metadata dataframe, with samples as rows. Contents of dataframe may vary, but it must contain a column called "Sample_ID" that contains the sample names which will match exactly to the sample names in the abund dataframe.
  
            #taxonomy = taxonomy dataframe, with OTUs as rows. Contents of dataframe may vary, but it must contain a column called "OTUNumber" that contains the OTU names which will match exactly to the OTU names in the abund dataframe.
  
  df1=as.data.frame(unlist(abund)) #unlist abund, put unlisted abund vector as a column in new dataframe, df1. The length of this column = ncol(abund) x nrow(abund)
  
  colnames(df1)="abund" #change column 1 name in df1 to "abund".
  
  df1$OTU=rep(rownames(abund),times=ncol(abund)) #add a new column, OTU, which contains the original OTU names ( rownames(abund) ), repeated by the number of columns ( times=ncol(abund) ). Because unlist will list values sequentially BY COLUMN, to create the "OTU" column the original OTU names in the original order are simply repeated X times, where X is the number of columns. This new column will have a length = nrow(abund) x ncol(abund), and now each value in the "abund" column has an associated value in the "OTU" column. 
  
  df1$Sample=rep(colnames(abund),each=nrow(abund))
  
  df2=merge(df1,taxonomy,by.x="OTU",by.y="OTUNumber")
  
  df3=merge(df2,metadata,by.x="Sample",by.y="Sample_ID")
  
  abund.longformat <<- df3
  
}