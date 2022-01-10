'ABCDom main script'


#check if your wd is the R folder
getwd()
#change the wd.project to the current R folder/current folder
wd.project<-getwd()

# activate packages
source(paste0(dirR,'/packages.R'))


#first workup DOC data and symbiondiniacae
source(paste0(dirR, '/DOC_workup.R'))
source(paste0(dirR, '/Symbiodiniaceae_density_workup.R'))
source(paste0(dirR, '/FCM_workup.R'))

<<<<<<< Updated upstream
#then analyze data
source(paste0(dirR, '/DOC_analysis.R'))
source(paste0(dirR, '/Symbiodiniaceae_analysis.R'))
source(paste0(dirR, '/FCM_analysis.R'))

#metabolite data cleanup
#will open new script that needs user input.
file.edit(paste0(Data.cleanup.folder, '/DataCleanup.R'))

#Symbiodiniaceae scripts
source(paste0(dirR, '/Symbiodiniaceae_density_workup.R'))

source(paste0(dirR, '/Symbiodiniaceae_analysis.R'))
>>>>>>> Stashed changes
