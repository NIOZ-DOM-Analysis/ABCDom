'Diversity.R

calculates the H diversity, evenness and number of metabolites from the cleaned data.
'
########
# define some things first:
# determine the dataset you want to use
ls.diversity <- list()
ls.diversity$df.div.T0 <- df.area %>% dplyr::filter(Timepoint_char == "T0" & Experiment == "ABCDOM" & Treatment != "Inoculum")
ls.diversity$df.div.Tend <- df.area %>% dplyr::filter(Timepoint_char == "Tend" & Experiment == "ABCDOM")

# what dataset did you use?
O <- 'area'

#today's date
D<-format(Sys.time(), '%Y%m%d')
# define the features
ls.diversity$features_T0 <- ls.diversity$df.div.T0[,M:ncol(ls.diversity$df.div.T0)]
ls.diversity$features_Tend <- ls.diversity$df.div.Tend[,M:ncol(ls.diversity$df.div.Tend)]

#first T0
# define the grouping
Treatment <- factor(ls.diversity$df.div.T0$Treatment, levels = c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", 'Bleached + Heated', "Ambient Water Control", "Heated Water Control", "not applicable"))
U <- Treatment
N <- 'Treatment'


ls.diversity$H1_T0 <- diversity(ls.diversity$features_T0, index = 'shannon', MARGIN = 1, base = exp(1))
ls.diversity$S1_T0 <- specnumber(ls.diversity$features_T0)
ls.diversity$J1_T0 <- ls.diversity$H1_T0/log(ls.diversity$S1_T0)

# create new dataframe to perform the statistics
ls.diversity$H1_T0 <- bind_cols(U, ls.diversity$H1_T0) %>% rename('U' = '...1') %>% rename('H' = '...2')
ls.diversity$S1_T0 <- bind_cols(U, ls.diversity$S1_T0) %>% rename('U' = '...1') %>% rename('S' = '...2')
ls.diversity$J1_T0 <- bind_cols(U, ls.diversity$J1_T0) %>% rename('U' = '...1') %>% rename('J' = '...2')

cost.col.fill<-c("blue","red", "white", "white", "grey70", "grey70", "grey70")
cost.col.line<-c("dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "lightgrey")

ggplot(ls.diversity$H1_T0, aes(x=U, y=H, group=U, color = U, fill = U))+
  stat_boxplot(geom = 'errorbar')+
  geom_boxplot(aes(), show.legend = TRUE, size = 1.2)+
  theme_classic()+
  theme(legend.position = "right")+
  scale_x_discrete(name="", guide = guide_axis(n.dodge = 2))+
  scale_fill_manual(values = cost.col.fill, name = "Treatment")+
  scale_color_manual(values = cost.col.line, name = "Treatment")+
  scale_y_continuous(name = "H' Diveristy")
ggsave("Metabolome_diversity_T0.jpeg" , path = dirFigs, dpi = 300, height=5, width=7)

ggplot(ls.diversity$S1_T0, aes(x=U, y=S, group=U, color = U, fill = U))+
  stat_boxplot(geom = 'errorbar')+
  geom_boxplot(aes(), show.legend = TRUE, size = 1.2)+
  theme_classic()+
  theme(legend.position = "right")+
  scale_x_discrete(name="", guide = guide_axis(n.dodge = 2))+
  scale_fill_manual(values = cost.col.fill, name = "Treatment")+
  scale_color_manual(values = cost.col.line, name = "Treatment")+
  scale_y_continuous(name = "Number of Features")
ggsave("Metabolome_number_of_features_T0.jpeg" , path = dirFigs, dpi = 300, height=5, width=7)

ggplot(ls.diversity$J1_T0, aes(x=U, y=J, group=U, color = U, fill = U))+
  stat_boxplot(geom = 'errorbar')+
  geom_boxplot(aes(), show.legend = TRUE, size = 1.2)+
  theme_classic()+
  theme(legend.position = "right")+
  scale_x_discrete(name="", guide = guide_axis(n.dodge = 2))+
  scale_fill_manual(values = cost.col.fill, name = "Treatment")+
  scale_color_manual(values = cost.col.line, name = "Treatment")+
  scale_y_continuous(name = "Pilou's evenness J=H'/log(S)")
ggsave("Metabolome_evenness_T0.jpeg" , path = dirFigs, dpi = 300, height=5, width=7)


# Now Tend
# define the grouping
Treatment <- factor(ls.diversity$df.div.Tend$Treatment, levels = c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", 'Bleached + Heated', "Ambient Water Control", "Heated Water Control", "not applicable"))
U <- Treatment
N <- 'Treatment'


ls.diversity$H1_Tend <- diversity(ls.diversity$features_Tend, index = 'shannon', MARGIN = 1, base = exp(1))
ls.diversity$S1_Tend <- specnumber(ls.diversity$features_Tend)
ls.diversity$J1_Tend <- ls.diversity$H1_Tend/log(ls.diversity$S1_Tend)

# create new dataframe to perform the statistics
ls.diversity$H1_Tend <- bind_cols(U, ls.diversity$H1_Tend) %>% rename('U' = '...1') %>% rename('H' = '...2')
ls.diversity$S1_Tend <- bind_cols(U, ls.diversity$S1_Tend) %>% rename('U' = '...1') %>% rename('S' = '...2')
ls.diversity$J1_Tend <- bind_cols(U, ls.diversity$J1_Tend) %>% rename('U' = '...1') %>% rename('J' = '...2')

cost.col.fill<-c("blue","red", "white", "white", "grey70", "grey70", "grey70")
cost.col.line<-c("dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "lightgrey")

ggplot(ls.diversity$H1_Tend, aes(x=U, y=H, group=U, color = U, fill = U))+
  stat_boxplot(geom = 'errorbar')+
  geom_boxplot(aes(), show.legend = TRUE, size = 1.2)+
  theme_classic()+
  theme(legend.position = "right")+
  scale_x_discrete(name="", guide = guide_axis(n.dodge = 2))+
  scale_fill_manual(values = cost.col.fill, name = "Treatment")+
  scale_color_manual(values = cost.col.line, name = "Treatment")+
  scale_y_continuous(name = "H' Diveristy")
ggsave("Metabolome_diversity_Tend.jpeg" , path = dirFigs, dpi = 300, height=5, width=7)

ggplot(ls.diversity$S1_Tend, aes(x=U, y=S, group=U, color = U, fill = U))+
  stat_boxplot(geom = 'errorbar')+
  geom_boxplot(aes(), show.legend = TRUE, size = 1.2)+
  theme_classic()+
  theme(legend.position = "right")+
  scale_x_discrete(name="", guide = guide_axis(n.dodge = 2))+
  scale_fill_manual(values = cost.col.fill, name = "Treatment")+
  scale_color_manual(values = cost.col.line, name = "Treatment")+
  scale_y_continuous(name = "Number of features")
ggsave("Metabolome_number_of_features_Tend.jpeg" , path = dirFigs, dpi = 300, height=5, width=7)

ggplot(ls.diversity$J1_Tend, aes(x=U, y=J, group=U, color = U, fill = U))+
  stat_boxplot(geom = 'errorbar')+
  geom_boxplot(aes(), show.legend = TRUE, size = 1.2)+
  theme_classic()+
  theme(legend.position = "right")+
  scale_x_discrete(name="", guide = guide_axis(n.dodge = 2))+
  scale_fill_manual(values = cost.col.fill, name = "Treatment")+
  scale_color_manual(values = cost.col.line, name = "Treatment")+
  scale_y_continuous(name = "Pilou's evenness J=H'/log(S)")
ggsave("Metabolome_evenness_Tend.jpeg" , path = dirFigs, dpi = 300, height=5, width=7)
