'16S_ASV_fasta_workup.R'

#read in data
asv_fasta <- read.csv(file.path(dirRAW, "16S", "all_multipletonsFilter_100.fasta.csv"))

#take every other cell value (ASV name), starting with row 1, and paste in a new df with colnames="ASV".
asv_fasta1 <- as.data.frame(asv_fasta[seq(from=1, to=nrow(asv_fasta), by=2),])
colnames(asv_fasta1) <- "ASV"

#take every other cell value (ASV fasta seq), starting with row 2, and paste as a new column in asv_fasta1 df, with each sequence corresponding to the correct asv.
asv_fasta1$fasta <- asv_fasta[seq(from=2, to=nrow(asv_fasta), by=2),]

#extract just ASV number and put in a new column
asv_fasta1$ASV_num <- as.data.frame(t(as.data.frame(strsplit(as.character(unlist(asv_fasta1$ASV)), split="_"))))[,2]

#export
write.csv(asv_fasta1, file.path(dirOutput, "asv_fasta1.csv"))
