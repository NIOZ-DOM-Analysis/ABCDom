#This function takes the weighted unifrac table that the CMAIKI 16S pipeline outputs, which contains pairwise comparisons and an associated unifrac values, and a metadata table and converts it to a square distance matrix.

generate.square.dist.matrix=function(dist.table,metadata,colnum) {
  #Inputs: dist.table = a dataframe corresponding to the distance matrix output from the CMAIKI pipeline. It should have 3 columns: "Tree#",	"Groups",	"WScore"
  #The CMAIKI distance matrix output is usually called otu_repr_"OTU cluster level".tre1.wsummary (weighted unifrac) or otu_repr_"OTU cluster level".uwsummary (unweighted unifrac).
  #Note that "OTU_cluster_level" is replaced by the % value you specified OTUs to be clustered at.
  #To get these CMAIKI distance matrix outputs into R, you must first open them in excel and save as a CSV. Then they can be read into R.
  
  #metadata = a dataframe corresponding to the metadata of your samples. Samples must be rows. Only include the samples for which you want to include in your 16S analysis.
  
  #colnum = the column number that corresponds to the sample names column. Note that sample names in metadata and dist.table must match exactly.
  
  #Note: vegan needs to be installed prior to running this function.
  
  #Split the "Groups" column of the dist.table into 2 columns so that each sample name in a pairwise comparison has a seperate column.
  dist.table.groups=strsplit(as.character(dist.table$Groups),split="-") #split the groups column by "-", moving the 2 samples that are being compared into seperate vectors within a list.
  
  for (i in 1:nrow(dist.table)) { #for each row in the dist.table (corresponding to a pairwise comparison)
    
    dist.table$Group1[i]=dist.table.groups[[i]][1] #make a new column called Group1, extract the character string from the unifrac.split.groups list that corresponds to the first sample in the pairwise comparison and add it to the Group1 column
    
    dist.table$Group2[i]=dist.table.groups[[i]][2] #do the same thing but with the second sample in the pairwsie comparison (Group 2)
  }
  
  #make an empty df with the number of rows and columns equal to the number of samples you want in the dist matrix. Pairwise distances from the dist.table df will be extracted and moved here.
  dist.matrix=as.data.frame(matrix(nrow=length(metadata[,1]),ncol=length(metadata[,1]))) #Use the input metadata df to count the number of samples (nrows and ncols) to include in the dist.matrix
  
  #name the rows and columns according to samples you want in the matrix. X corresponds to the number of the column containing the desired sample names.
  rownames(dist.matrix)=metadata[,colnum]
  colnames(dist.matrix)=metadata[,colnum]
  
  #for each row (pairwise comparison) in the dist.table df, extract the W score and then copy it into the correct row and column of the empty dist.matrix.
  for (i in 1:nrow(dist.table)) { #for each row in the weighted.unifrac df
    
    dist.matrix[rownames(dist.matrix)==dist.table[,4][i], #subset the dist.matrix to find the row that corresponds to the target Group1 samplename in the dist.table df. 
                
                colnames(dist.matrix)==dist.table[,5][i]]=dist.table[,3][i] #subset the dist.matrix df to find the column that corresponds to the target Group2 samplename in the dist.table.
    #the dist.matrix has now been subset to the cell that has the row and column values corresponding to the Group1 and Group2 values for the pairwise comparison in the dist.table df. For this cell, add the corresponding W score from the comparison.
  }
  
  dist.matrix[is.na(dist.matrix)]=0 #replace NAs with 0s
  
  dist.matrix1=dist.matrix #copy distance matrix df
  
  #Fill diagnols 
  for (i in 1:ncol(dist.matrix)) { #for each column in dist.matrix
    
    dist.matrix1[,i]=dist.matrix[,i]+as.vector(t(dist.matrix[i,])) #make the corresponding column in dist.matrix1 equal to the values of the column in dist.matrix AND the values of the corresponding transposed row in dist.matrix.
  } #This ensures that if values on both sides the diagnol of the dist.matrix will be captured in the final dist.matrix. Simply transposing all rows to columns (or vice versa) would erase any values across the diagnol. Adding columns and transposed rows avoids this issue. Columns and rows can be added because no value is duplicated across the diagnol (AKA no pairwise comparison value is present twice at this point).
  
  dist.matrix1=dist.matrix1[as.logical(apply(dist.matrix1,1,FUN=sum)),as.logical(apply(dist.matrix1,2,FUN=sum))] #Sum the vaues of each column and row. Convert to logical so all sums with non-0 values = TRUE and all sums with 0 values = FALSE. Subset the dist.matrix1 so that any rows or columns whose sum=0 are removed. This removes samples which are in the metadata but failed to sequence and/or pass through the pipeline.
  
  #convert to matrix
  square.dist.matrix=as.dist(dist.matrix1)
  
  dist.matrix1 <<- dist.matrix1
  square.dist.matrix <<- square.dist.matrix
  
}