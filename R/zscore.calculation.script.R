zscore.calculation=function(x) {
  #Inputs: x = dummy variable that function will be run on. Must be a vector of numbers.
  
  zscoredat=x #duplicate vector x, save as a new df zscoredat.
  
  for (i in 1:length(x)) { #for each value in vector x
    zscoredat[i]=(x[i]-mean(x))/sd(x) #subset the input data to the corresponding position, and make it equal to the zscored value of the original cell value at that position.
    #zscoring is done by subtracting the mean value for a vector from the value of the vector at position i, and then dividing by the SD for the vector.
  }
  zscoredat1 <<- zscoredat #save output as new df, zscoredat1
}
