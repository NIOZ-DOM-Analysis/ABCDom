'Procrustes script'


# we first have to give the metabolomes and 16S data the same form/shape/samples etc.
# define which dataframes we are using
pro.16 <- unifrac.dist.matrix.tend
pro.Meta <- df.area.ABCT0 %>% dplyr::filter(Experiment == "ABCDOM")

rownames(pro.16)[1]<-"ABC_055"
rownames(pro.16)[4]<-"ABC_059"
meta.pro.16<- metadata %>% dplyr::filter(`Sample Name` %in% rownames(pro.16)) %>% select(`Sample Name`, Bottle_NR)

#give pro.16 the right name
pro.16 <- pro.16 %>% rownames_to_column("Sample Name")
pro.16 <- right_join(meta.pro.16, pro.16, by = "Sample Name")
colnames(pro.16)[3:ncol(pro.16)]<-meta.pro.16$Bottle_NR
pro.16 <- pro.16 %>% column_to_rownames("Bottle_NR")%>% select(-"Sample Name")

#now fix the metabolome dataset
pro.Meta<- pro.Meta %>% dplyr::filter(Bottle_NR %in% colnames(pro.16))
#remove the one double sample
pro.Meta<- pro.Meta %>% dplyr::filter(`File Name` != "Mo_2019_324_B")

pro.Meta.dist<-as.matrix(vegdist(pro.Meta[M:ncol(pro.Meta)]))
pro.Meta.dist<-as.data.frame(pro.Meta.dist)
rownames(pro.Meta.dist)<-pro.Meta$Bottle_NR
colnames(pro.Meta.dist)<-pro.Meta$Bottle_NR

pro<-procrustes(pro.16, pro.Meta.dist, symmetric = FALSE)
pro
pro.m2<-procrustes(pro.16, pro.Meta.dist, symmetric = TRUE)
pro.m2

plot(pro, kind = 1, type = "text")
plot(pro, kind = 2)

plot(pro.m2, kind = 1, type = "text")
plot(pro.m2, kind = 2)

protest(pro.16, pro.Meta.dist, scores = "sites", permutations = 999)
#there is significance!



#now with relabundance and bray dissim
pro.16.bray <- relabund %>% select(rownames(unifrac.dist.matrix.tend))
pro.16.bray<-as.data.frame(t(pro.16.bray))
rownames(pro.16.bray)[1]<-"ABC_055"
rownames(pro.16.bray)[4]<-"ABC_059"
pro.16.bray<-as.matrix(vegdist(pro.16.bray))
rownames(pro.16.bray)<-pro.Meta$Treatment
colnames(pro.16.bray)<-pro.Meta$Bottle_NR

pro<-procrustes(pro.16.bray, pro.Meta.dist, symmetric = FALSE)
pro
pro.m2<-procrustes(pro.16.bray, pro.Meta.dist, symmetric = TRUE)
pro.m2

plot(pro, kind = 1)
plot(pro, kind = 2)

plot(pro.m2, kind = 1)
plot(pro.m2, kind = 2)

protest(pro.16, pro.Meta.dist, scores = "sites", permutations = 999)
#there is significance!
