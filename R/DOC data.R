# DOC data (run symbiont data first!)
library(tidyverse)

DOC_raw <- read_csv(file.path(dirRAW, "DOC", "DOC_with_SampleName.csv"))

# remove the following samples based on Carlson lab instructions because septa upside down
filter.out <- c("ABC_037", "ABC_048", "ABC_051", "ABC_054", "ABC_092", "ABC_097", "ABC_108")

DOC_dat <- right_join( full_metadata, DOC_raw, by = "Sample Name" )%>%
  filter(!`Sample Name`%in% filter.out)

tmp<-DOC_dat %>%
  filter(Timepoint_char == "T0") %>%
  filter(is.na(Bleaching_Status)) %>%
  group_by(Treatment) %>%
  summarise(background_DOC = mean(uMC))

DOC_dat <- DOC_dat %>%
  mutate(Background_DOC = case_when(Stress_status == "Thermal Stress" ~ tmp$background_DOC[2],
                                    Stress_status == "Ambient" ~ tmp$background_DOC[1])) %>%
  mutate(Control_Corrected_DOC = uMC - Background_DOC)%>%
  mutate(Control_Corrected_DOC_flux = Control_Corrected_DOC / 3.5)


# put together with symbiont data
DOC_dat <- left_join(DOC_dat, sym_dat_aquaria, by = "PLANC_aquaria") %>%
  mutate(DOC_SA_Normalized = uMC / sum_SA_no_outliers) %>%
  mutate(Control_Corrected_DOC_SA_Normalized = Control_Corrected_DOC /sum_SA_no_outliers) %>%
  mutate(Control_Corrected_DOC_flux_SA_Normalized = Control_Corrected_DOC_flux / sum_SA_no_outliers) %>%
  mutate(Control_Corrected_DOC_flux_SA_Normalized_dm2 = Control_Corrected_DOC_flux_SA_Normalized*100)

DOC_dat$Treatment.x <- factor(DOC_dat$Treatment.x, levels = c("Control", "Heated", "Bleached", "Bleached + Heated", "Negative Control", "Negative Control + Heated"))

#create manual
cost.col.fill<-c("dodgerblue1","firebrick1", "white", "white", "grey70", "grey70", "grey70")
cost.col.line<-c("dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "lightgrey")
cost.shape <- c(21,24)

#select T0 data

tmp <- DOC_dat %>%
  filter(Timepoint_char == "T0")

tmp %>%
  mutate(coral_water = case_when(is.na(Bleaching_Status)~ "control",
                                 !is.na(Bleaching_Status)~ "Control")) %>%
  kruskal_test(uMC~coral_water)
# water controls significantly lower than corals

tmp %>%
  kruskal_test(uMC~Treatment.x)
# no significant difference between all treatments

#visualize
figS2<-ggplot(tmp, aes(x=Treatment.x,y=uMC,color=Treatment.x,fill=Treatment.x))+
  stat_boxplot(geom = 'errorbar', linewidth = 2)+
  geom_boxplot(linewidth = 1.2)+
  # geom_point(size = 3)+
  scale_color_manual(values=cost.col.line)+
  # scale_fill_manual(values=cost.col.fill, guide = guide_legend(override.aes = list(size = 1)))+
  scale_fill_manual(values=cost.col.fill)+
  theme_classic()+
  # theme(legend.key.height=unit(0.5,"in"))+
  theme(legend.position = "none")+
  ylab("DOC (µM)")+
  xlab("")
#use this figure for a potential supplement
#ggsave('FigS2_DOC_per treatment.jpeg', path = dirFigs, width = 9, height = 5.5, units = "in", dpi = 320)

#write_csv(tmp, file.path(dirOutput, "FigS2_data.csv"))


tmp <- DOC_dat %>%
  filter(Timepoint_char == "T0")%>%
  filter(!is.na(Bleaching_Status))

fig2A<-ggplot(tmp,aes(x=Treatment.x,y=Control_Corrected_DOC_flux_SA_Normalized_dm2 ,color=Treatment.x,fill=Treatment.x))+
  stat_summary(geom="bar", fun="mean", size=1.5)+
  stat_summary(fun.data = mean_se, geom = "linerange")+
  geom_point(pch=21, size=2, stroke= 2)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill)+
  theme_classic()+
  theme(legend.position = "none")+
  ylab (expression("Surface Area Normalized DOC flux (µM (dm"^2*")"^-1*"h"^-1*")"))+
  xlab("")
#ggsave('fig2A_DOC_flux Surface area normalized_control corrected_per treatment barplot.jpeg', path = dirFigs, width = 7.5, height = 5.5, dpi = 300)
#Use as panel A for Figure 2

#write_csv(tmp, file = file.path(dirOutput, "Fig2A_data.csv"))

tmp <- DOC_dat %>%
  filter(Timepoint_char == "T0")

#Check the distribution of raw DOC values.
hist(tmp$uMC) #not very normal
hist(log10(tmp$uMC)) #still not normal

# Since raw DOC values don't seem to be normal, test the effect of Treatment using a K-W test.
tmp %>%
  kruskal_test(uMC~Treatment.x)
 #marginally significant

#double check with an ANOVA
tmp %>%
  anova_test(uMC~Treatment.x)

# now do the stats on the flux
tmp %>%
  kruskal_test(Control_Corrected_DOC_flux_SA_Normalized_dm2~Treatment.x)
#marginally significant

#summary results
tmp %>%
  group_by(Treatment.x)%>%
  summarise(mean(Control_Corrected_DOC_flux_SA_Normalized_dm2))


