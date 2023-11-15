#symbiont analysis
library(vegan)
library(rstatix)
library(tidyverse)
library(multcomp)
library(dplyr)

mutate <- dplyr::mutate
select <- dplyr::select

# Read in data.
PLANC_Blastate_FCM <- read_csv(file.path(dirRAW, "Sym_Counts", "PLANC_Blastate_FCM.csv"), ) %>% #load the raw FCM data
  mutate(Well = replace(Well, Well == "H8(1)", "H8"))
PLANC_Blastate_FCM_Sample_Map <- read_csv(file.path(dirRAW, "Sym_Counts", "PLANC_Blastate_FCM_Sample_Map.csv"), ) #load the FCM sample map data
PLANC_nubbin_sizing_notes <- read_csv(file.path(dirRAW, "Sym_Counts", "PLANC_nubbin_sizing_notes.csv"), ) #load the nubbin size notes
PLANC_metadata <- read_csv(file.path(dirRAW,"Sym_Counts","Symbiont_metadata.csv"))

# Merge the three data frames.

sym_dat<- full_join(PLANC_Blastate_FCM_Sample_Map, PLANC_Blastate_FCM, by = c("Plate", "Well"))
sym_dat<- merge(sym_dat, PLANC_nubbin_sizing_notes, by = join_by("Sample"== "Code"))

#make two columns concentration and events based on  if we have to use Sym-SCC data
sym_dat <- sym_dat %>%
  mutate(sym.con = case_when(is.na(Notes.x) ~  `Sym-FSC Events/μL(V)`,
                             .default = `Sym-SSC Events/μL(V)`), .before = Notes.x) %>%
  mutate(sym.events = case_when(is.na(Notes.x) ~  `Sym-FSC Events`,
                             .default = `Sym-SSC Events`), .before = Notes.x) %>%
  mutate(sym.total = sym.con*Vol*1000) %>%
  mutate(sym.SA = sym.total/SA) %>%
  mutate(log10.sym.SA = log10(sym.SA))


# sum the sym counts and SA values by aquaria and timepoint, without outlier1
tmp <- sym_dat %>%
  filter(!is.na(log10.sym.SA)) %>%
  filter(is.na(Outlier1)) %>%
  group_by(Aquaria_Timepoint_Heat) %>%
  summarise(sum_sym_total = sum(sym.total),
            sum_SA= sum(SA),
            sym_SA_norm = sum_sym_total / sum_SA)

# sum the sym counts and SA values by aquaria and timepoint, without outlier
tmp2 <- sym_dat %>%
  filter(!is.na(log10.sym.SA)) %>%
  filter(is.na(Outlier)) %>%
  group_by(Aquaria_Timepoint_Heat) %>%
  summarise(sum_sym_total = sum(sym.total),
            sum_SA= sum(SA),
            sym_SA_norm = sum_sym_total / sum_SA)

sym_dat_aquaria <- full_join(tmp, tmp2, by = "Aquaria_Timepoint_Heat", suffix = c("_no_outliers1", "_no_outliers"))

sym_dat_aquaria <- right_join(PLANC_metadata, sym_dat_aquaria, by = join_by(PLANC_aquaria == Aquaria_Timepoint_Heat))

# visualize Symbiont at collection and do stats with stringent outliers (outlier)
tmp <- sym_dat %>%
  filter(!is.na(log10.sym.SA)) %>%
  filter(is.na(Outlier)) %>%
  filter(Timepoint == "T0") %>%
  unite(Species_Bl, c("Collection_Bleaching_Level1", "Species"), remove = F)

tmp$Collection_Bleaching_Level1 <- factor(tmp$Collection_Bleaching_Level1, levels =  c("HE", "BL"))
tmp$Species_Bl <- factor(tmp$Species_Bl)

Fig1C <- ggplot(tmp,(aes(x=Collection_Bleaching_Level1,y=log10.sym.SA,color=Collection_Bleaching_Level1, fill=Collection_Bleaching_Level1)))+
  facet_wrap(.~Species, scales="fixed")+
  stat_boxplot(geom = 'errorbar', linewidth = 0.81)+
  geom_boxplot(size=1)+
  scale_color_manual(name = "Bleaching Status\nat Collection", labels = c("Unbleached", "Bleached"), values=c("dodgerblue3", "dodgerblue3"))+
  scale_fill_manual(name = "Bleaching Status\nat Collection", labels = c("Unbleached", "Bleached"), values=c("dodgerblue1", "white"))+
  scale_x_discrete(name = NULL, labels=NULL)+
  scale_y_continuous(limits = c(3,6), name = expression("Log"[10]*" Symbiodiniaceae cells per cm"^2))+
  # theme(text=element_text(size=45), plot.margin=unit(c(1,1,1,1), units = "in"), strip.text.x=element_text(size=45))+
  theme_classic()
Fig1C
#ggsave('Fig1C_Symbiont cells per cm2_per treatment vMARCH barplot.jpeg', path = dirFigs, dpi = 300, width=6, height=3.5, units = "in")

#write_csv(tmp, file.path(dirOutput,  "Fig1C_data.csv"))

tmp %>%
  anova_test(log10.sym.SA ~ Collection_Bleaching_Level1 * Species)

tmp %>%
  tukey_hsd(log10.sym.SA ~ Collection_Bleaching_Level1 * Species)

tmp %>%
  anova_test(log10.sym.SA ~ Species_Bl)

tmp %>%
  tukey_hsd(log10.sym.SA ~ Species_Bl)


amod <- aov(log10.sym.SA ~ Species_Bl, data = tmp)
tuk <- glht(amod, linfct = mcp(Species_Bl = "Tukey"))
### extract information
tidy(cld(tuk)) #reorder de letters


# visualize Symbiont after pre-treatment #also with most stringent outlier removal
tmp <- sym_dat_aquaria %>%
  filter(Timepoint_char_PlanC == "T7")
tmp$Treatment <- factor(tmp$Treatment, levels = c("Control", "Heated", "Bleached", "Bleached + Heated"))

cost.col.fill<-c("dodgerblue1","firebrick1", "white", "white", "grey70", "grey70")
cost.col.line<-c("dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3")


Fig1D <-ggplot(tmp,aes(x=Treatment,y=log10(sym_SA_norm_no_outliers),color=Treatment,fill=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 0.8)+
  geom_boxplot(size=1)+
  # geom_point(size=3, pch=21)+
  scale_color_manual(values=cost.col.line)+
  scale_fill_manual(values=cost.col.fill)+
  scale_x_discrete(name = "")+
  theme_classic()+
  # theme(text=element_text(size=34),legend.key.height=unit(2,"cm"),axis.text.x=element_blank(), plot.margin=unit(c(1,1,1,1.5), "cm"))+
  ylab(expression(paste("Log10 Symbiodiniaceae\ncells per cm^2")))+
  labs(color="Treatment",fill="Treatment")
Fig1D
#ggsave('fig1D_Symbiont cells per cm2_per treatment vMARCH barplot.jpeg', path = dirFigs, dpi = 300, width=15, height=9, units = "in")

#write_csv(tmp, file.path(dirOutput,  "Fig1D_data.csv"))

tmp %>%
  anova_test(sym_SA_norm_no_outliers ~ Treatment)
tmp %>%
 tukey_hsd(sym_SA_norm_no_outliers ~ Treatment)
amod <- aov(sym_SA_norm_no_outliers ~ Treatment, data = tmp)
tuk <- glht(amod, linfct = mcp(Treatment = "Tukey"))
### extract information
tidy(cld(tuk))

# seperate for heating and treatment
tmp %>%
  anova_test(sym_SA_norm_no_outliers ~ Temperature*Stress_status)

