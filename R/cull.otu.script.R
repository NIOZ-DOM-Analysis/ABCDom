cull.otu=function(relabund.df, min.num, min.abund, min.single.abund) {
  #Inputs:  relabund.df = dataframe containing ONLY relative abundance data, no metadata or other info. Samples in rows and OTUs in columns.
            
            #min.num = the minimum number of samples an OTU needs to be present in to not be culled.
            
            #min.abund = the minimum relative abundance an OTU needs to have in (the min.num) samples to not be culled.
 
            #min.single.abund = the minimum relative abundance an OTU needs to have in a SINGLE sample to not be culled.
   
  sub=c() #create a empty vector
  
  cull=function(x) { #make a function that says for any input, generate a logical vector of TRUEs and FALSEs that will be used for subsetting, selecting OTUs that
    sub=ifelse(length(x[x>=min.abund])>=min.num #have a relabund>"min.abund" in "min.num" samples 
               | length(x[x>=min.single.abund])>0,TRUE,FALSE) #or have a relabund>"min.single.abund" in at least one sample
    return(sub)
  }
  
  cull.vec=apply(relabund.df,2,FUN=cull) #apply cull function to relabund.df, save output as a vector.
  relabund.df.cull=relabund.df[,cull.vec] #Use cull.vec to subset the columns of relabund.df for OTUs that passed the cull threshold.
  
  relabund.df.cull<<-relabund.df.cull
}