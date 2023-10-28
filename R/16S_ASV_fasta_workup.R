"16S_workup and analysis"


# library(phyloseq)
# library(DESeq2)

library(tidyverse)

#read in data
asv_fasta <- read_csv(file.path(dirRAW, "16S", "all_multipletonsFilter_100.fasta.csv"))

#take every other cell value (ASV name), starting with row 1, and paste in a new df with colnames="ASV".
asv_fasta1 <- asv_fasta %>%
  mutate(number = rep(1:(nrow(asv_fasta)/2), each=2)) %>%
  mutate(names = case_when(grepl("esv", Row1) ~ "ASV",
                           .default = "fasta") ) %>%
  pivot_wider(names_from = names, values_from = Row1) %>%
  select(-number) %>%
  separate(ASV, c("esv", "ASV_num"), sep = "_", remove = FALSE) %>%
  select(-esv)

