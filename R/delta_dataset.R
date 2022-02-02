#create a metabolite delta dataset

df.area.ABC_0<- df.area %>% dplyr::filter(Experiment == "ABCDOM") %>%
  dplyr::filter(Timepoint_char == "T0") %>%
  dplyr::filter(File.Name != "Mo_2019_324_B") %>%
  dplyr::filter(File.Name != "Mo_2019_543") %>%
  dplyr::arrange(Bottle_NR)

df.area.ABC_36<- df.area %>% dplyr::filter(Experiment == "ABCDOM") %>%
  dplyr::filter(Timepoint_char == "Tend") %>%
  dplyr::arrange(Bottle_NR)

df.delta <- df.area.ABC_36[M:ncol(df.area.ABC_36)] - df.area.ABC_0[M:ncol(df.area.ABC_0)]
df.delta <- bind_cols(df.area.ABC_0[1:(M-1)], df.delta)
df.delta$Timepoint_char <- "Delta"

write_csv(df.delta, paste0(dirOutput, "\\df.delta.csv"))

# make NMDS plot
#non metric multidimentional scaling, all data
ord.mod.area<- metaMDS(df.delta[M:ncol(df.delta)], distance = 'euclidean', k=3, autotransform = FALSE, noshare = FALSE, wascores = FALSE)
ordi.plot1.1<-ordiplot(ord.mod.area, choices = c(1,2))
sites.long1<-cbind(df.delta[1:(M-1)], ord.mod.area$points)

#create object/list to store the sites.lists in.
ordi.data$NMDS.delta <- sites.long1

#make the treatment a factor to be ordered, and define the color scheme
cost.col.fill<-c("blue","red", "white", "white", "grey70", "grey70", "grey70")
cost.col.line<-c("dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "lightgrey")
fact.all.treat<-factor(ordi.data$NMDS.delta$Treatment, levels = c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", 'Bleached + Heated', "Ambient Water Control", "Heated Water Control"))
cost.shape <- c(21,24)

ggplot() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  xlab("NMDS3") +
  ylab("NMDS2") +
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data=ordi.data$NMDS.delta,
             aes(x=MDS3, y=MDS2, fill = fact.all.treat, colour=fact.all.treat, shape= Experiment),
             size=5, stroke=1.5) +
  scale_shape_manual(values = cost.shape, name = "Experiment")+
  scale_color_manual(values = cost.col.line, name = "Treatment")+
  scale_fill_manual(values = cost.col.fill, name = "Treatment", guide = guide_legend(override.aes = list(shape = 21)))+
  ggtitle("Manhattan Delta")+
  coord_fixed(ratio=1)
ggsave("NMDS_manhattan_dela.jpg", path = dirFigs, width = 6.75, height = 5, units = "in", dpi = 320)

scatterplot3d(ordi.data$NMDS.delta$MDS1, ordi.data$NMDS.delta$MDS2, ordi.data$NMDS.delta$MDS3,
              pch = 21, color = cost.col.line[fact.all.treat],
              bg = cost.col.fill[fact.all.treat], type = "h", angle = 80)
