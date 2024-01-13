#procrustus

#first analyze 16S and metabolome data

'Procrustes script'
library(cowplot)
library(tibble)

# we first have to give the metabolomes and 16S data the same form/shape/samples etc.
# define which dataframes we are using
pro.16 <- unifrac.dist.matrix.tend
pro.Meta <- df.area %>%
  filter(Timepoint_char == "T0")

rownames(pro.16)[1]<-"ABC_055"
rownames(pro.16)[4]<-"ABC_059"
meta.pro.16<- metadata %>%
  dplyr::filter(`Sample Name` %in% rownames(pro.16)) %>%
  select(`Sample Name`, Bottle_NR)

#give pro.16 the right name
pro.16 <- pro.16 %>% rownames_to_column("Sample Name")
pro.16 <- right_join(meta.pro.16, pro.16, by = "Sample Name")
colnames(pro.16)[3:ncol(pro.16)]<-meta.pro.16$Bottle_NR
pro.16 <- pro.16 %>% column_to_rownames("Bottle_NR")%>% select(-"Sample Name")

#now fix the metabolome dataset
pro.Meta<- pro.Meta %>% dplyr::filter(Bottle_NR %in% colnames(pro.16))

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

#lets first make a nmds and then do procrustus on it.
nmds.pro16<-monoMDS(pro.16)
nmds.proMeta<-monoMDS(pro.Meta.dist)

nmds.procrus<-procrustes(nmds.pro16, nmds.proMeta)
nmds.procrus
summary(nmds.procrus)
plot(nmds.procrus)
plot(nmds.procrus, kind=2)
residuals(nmds.procrus)
protest(nmds.pro16, nmds.proMeta, scores = "sites", permutations = 999)


#now with relabundance and bray dissim
pro.16.bray <- relabund %>% select(rownames(unifrac.dist.matrix.tend))
pro.16.bray<-as.data.frame(t(pro.16.bray))
rownames(pro.16.bray)[1]<-"ABC_055"
rownames(pro.16.bray)[4]<-"ABC_059"
pro.16.bray<-as.matrix(vegdist(pro.16.bray, method = "bray"))
rownames(pro.16.bray)<-pro.Meta$Treatment
colnames(pro.16.bray)<-pro.Meta$Bottle_NR



#make the procrustus plot in GGplot!
library(ggplot2)
library(grid)
# x = 16S y= metabolome
ctest <- data.frame(rda1=nmds.procrus$X[,1], rda2=nmds.procrus$X[,2],
                    Treatment=pro.Meta$Treatment, Dataset = "16S",
                    xstart=nmds.procrus$Yrot[,1], ystart=nmds.procrus$Yrot[,2])
ctest2 <- data.frame(rda1=nmds.procrus$Yrot[,1], rda2=nmds.procrus$Yrot[,2],
                     Treatment=pro.Meta$Treatment, Dataset = "Metabolome")
ctest<-bind_rows(ctest, ctest2)
rot <- nmds.procrus$rotation


cost.col.fill<-c("dodgerblue1","firebrick1", "white", "white", "grey70", "grey70", "grey70")
cost.col.line<-c("dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "lightgrey")
fact.all.treat<-factor(ctest$Treatment, levels = c("Control", "Heated", "Bleached","Bleached + Heated","Negative Control", "Negative Control + Heated"))




ggplot(ctest) +
  geom_hline(aes(yintercept = 0), linetype = 2)+
  geom_vline(aes(xintercept = 0), linetype = 2)+
  geom_abline(aes(intercept = 0, slope = rot[1,2]/rot[1,1]), linetype = 1)+
  geom_abline(aes(intercept = 0, slope = rot[2,2]/rot[2,1]), linetype = 1)+
  geom_point(aes(x=rda1, y=rda2, fill=fact.all.treat, colour=fact.all.treat, shape = Dataset), size=3, stroke=1.5) +
  # geom_point(aes(x=xrda1, y=xrda2, fill=fact.all.treat, colour=fact.all.treat), size=3, stroke=1.5, shape = 24) +
  geom_segment(aes(x=xstart,y=ystart,xend=rda1,yend=rda2,colour=fact.all.treat),arrow=arrow(length=unit(0.2,"cm")), show.legend = FALSE)+
  scale_color_manual(values=cost.col.line, name = "Treatment", labels=c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  scale_fill_manual(values=cost.col.fill, name = "Treatment", guide = guide_legend(override.aes = list(shape = 21)), labels=c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  scale_shape_manual(values = c(21, 22), name = "Dataset", labels = c("16S", "Metabolome"))+
  scale_x_continuous(name = "Dimension 1")+
  scale_y_continuous(name = "Dimension 2")+
  theme_bw()
ggsave("16S_metabolome_nmds_procrustes.jpeg", path = dirFigs,  width = 9, height = 5.5, units = "in", dpi = 320)

fact.all.treat.proc<-factor(ctest$Treatment, levels = c("Control", "Heated", "Bleached","Bleached + Heated","Negative Control", "Negative Control + Heated"))

fig4B<-ggplot(ctest) +
  geom_hline(aes(yintercept = 0), linetype = 2)+
  geom_vline(aes(xintercept = 0), linetype = 2)+
  geom_abline(aes(intercept = 0, slope = rot[1,2]/rot[1,1]), linetype = 1)+
  geom_abline(aes(intercept = 0, slope = rot[2,2]/rot[2,1]), linetype = 1)+
  geom_point(aes(x=rda1, y=rda2, fill=fact.all.treat.proc, colour=fact.all.treat.proc, shape = Dataset), size=3, stroke=1.5) +
  # geom_point(aes(x=xrda1, y=xrda2, fill=fact.all.treat, colour=fact.all.treat), size=3, stroke=1.5, shape = 24) +
  geom_segment(aes(x=xstart,y=ystart,xend=rda1,yend=rda2,colour=fact.all.treat.proc),arrow=arrow(length=unit(0.2,"cm")), show.legend = FALSE)+
  scale_color_manual(values=cost.col.line, name = "Treatment", labels=c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  scale_fill_manual(values=cost.col.fill, name = "Treatment", guide = guide_legend(override.aes = list(shape = 21)), labels=c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))+
  scale_shape_manual(values = c(21, 22), name = "Dataset", labels = c("16S", "Metabolome"))+
  scale_x_continuous(name = "Dimension 1")+
  scale_y_continuous(name = "Dimension 2")+
  coord_fixed(ratio=1.48)+
  theme_bw()
fig4B

write_csv(ctest, file.path(dirOutput, "Fig4B_data.csv"))

fig4<-plot_grid(fig4A, fig4B, labels="AUTO", nrow =1, rel_heights = c(1, 1), align = c("hv"), axis = "l")
fig4
ggsave("fig4.jpg", plot = fig4, path = dirFigs, width = 14, height = 5, units = "in", dpi = 320)

#Next, run a mantel test between pro.16 and pro.Meta
mantel.results <- mantel(pro.16, pro.Meta.dist)
mantel.results
