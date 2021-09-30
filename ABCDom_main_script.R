'ABCDom main script'


#check if your wd is the R folder
getwd()
#change the wd.project to the current R folder/current folder
wd.project<-getwd()

# activate packages
source(paste0(dirR,'/packages.R'))

#first analyze DOC values
source(paste0(dirR, '/DOC_analysis.R')) # have to make this
