'DOC_Analysis.R

goal:  merge with metadata and plot the DOC data'

# import DOC data from the DOC folder in the RAW folder
# change the mu sign to u (uM) otherwise you will get an error and add the matching Sample Names like in the metadata

DOC_data<-read_csv(file.path(dirRAW, "DOC", "DOC_with_SampleName.csv"), )

# add metadata
DOC_data<-right_join(metadata, DOC_data, by = "Sample Name")


# plot something?
