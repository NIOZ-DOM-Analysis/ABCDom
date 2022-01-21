log2.fold.change = function(df,col.num,control.name,col.range) {
  
  #Inputs:  df is a data frame (relative abundance or abundance, with OTUs/features as columns and samples as rows). The df must contain one column corresponding to sample treatment metadata in which controls and treatments are specified. All OTU/feature columns must come before any metadata columns.
  
            #col.num is an integer that corresponds to the column number which contains the sample treatment metadata info
  
            #control.name is a character string corresponding to the name given to the control treatments.
  
            #col.range is a integer that corresponds to the final column number in df that contains OTU/feature information.
  
  log2.fold.df=df #create a duplicate of df called "log.2.fold.df". Output log2 fold-change values will be stored in this df.
  
  controlvec=df[,col.num]==control.name #subset df for column = col.num, which corresponds to treatment column. Then generate a logical string for only elements that == control.name, a character string corresponding to the name of the control treatment. This now generates a logical string that selects for only the positions corresponding to the control treatments.
  
  for (j in 1:w) { #for each column in df that corresponds to an OTU/feature, ranging from 1 to col.range,
    
    column=df[,j] #subset df for column j, save as new vector "column"
    
    column.log2.fold=column #duplicate column, save as output column "column.log.2.fold". Log2 fold-change values will be stored in this vector.
    
    for (i in 1:length(column)) { #for each element in column
      
      column.log2.fold[i] = log2(column[i]/mean(column[controlvec])) # calculate the log2fold ratio of the element to the mean control value (calculated by subsetting column by controlvec and calculating the mean). Add this value into the output vector "column.log2.fold" at the corresponding position.
      
    }
    
    log2.fold.df[,j] = column.log2.fold #Add the output column column.log2.fold back into the corresponding jth column position in the out put df "log2.fold.df"
    
  }
  
  log2.fold.df1 <<- log2.fold.df #store output as new df.
  
}