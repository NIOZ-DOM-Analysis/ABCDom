"FCM microbial growth data prep and analysis"
library(vegan)
library(rstatix)
library(multcomp)
library(magrittr)
library(tidyverse)

mutate <- dplyr::mutate
select <- dplyr::select


PlateII.Map <- read_csv(file.path(dirRAW, "FCM", "abcDOM.SYBR.NOV.15.Plate.II.Map.csv")) %>%
  select(SampleID, Well) %>%
  rename("Sample Name" = "SampleID")
PlateI.Map <- read_csv(file.path(dirRAW, "FCM", "abcDOM.SYBR.NOV.15.Plate.I.Map.csv")) %>%
  select(SampleID, Well) %>%
  rename("Sample Name" = "SampleID")
PlateII <- read_csv(file.path(dirRAW, "FCM", "abcDOM.SYBR.NOV.15.Plate.II.csv")) %>% rename ("Well" = "Sample Name")
PlateI <- read_csv(file.path(dirRAW, "FCM", "abcDOM.SYBR.NOV.15.Plate.I.csv")) %>% rename ("Well" = "Sample Name")

# merge the plate and map, then filter for polygon

FCMI <- PlateI %>% select (Well, `Gate Name`, Concentration) %>%
  left_join(., PlateI.Map, by = c("Well")) %>% filter( `Gate Name` == "SYBR Polygon")

FCMII <- PlateII %>% select (Well, `Gate Name`, Concentration) %>%
  left_join(., PlateII.Map, by = c("Well")) %>% filter( `Gate Name` == "SYBR Polygon")

FCM.dat<-bind_rows(FCMI, FCMII) %>%
  right_join(full_metadata, ., by = c("Sample Name")) %>%
  filter(!is.na(Treatment)) %>%
  mutate(Concentration = as.numeric(as.character(sub(",","", Concentration))))

rm(PlateI, PlateI.Map, PlateII, PlateII.Map, FCMI, FCMII)



# # Check distribution of concentration data.
hist(FCM.dat$Concentration) #not normal, strong right skew
hist(sqrt(FCM.dat$Concentration)) #more normal

#calculate the mean and SE for bleaching susceptability treatments and timepoint

FCM.summ <- FCM.dat %>% group_by(Treatment, Timepoint, Timepoint_char) %>%
  summarise( mean = mean(Concentration), sd = sd(Concentration), n = n(), se = sd/sqrt(n))

Treatment <- factor(FCM.summ$Treatment, levels = c("Control", "Heated", "Bleached","Bleached + Heated","Negative Control", "Negative Control + Heated"))
cost.col.fill<-c("dodgerblue1","firebrick1", "white", "white", "grey70", "grey70")
cost.col.line<-c("dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3")


fig2B<-ggplot(FCM.summ,aes(x=Timepoint,y=mean, fill=Treatment,shape=Treatment,color=Treatment,group=Treatment))+
  geom_point(size=8)+
  geom_line(linewidth=3)+
  scale_shape_manual(values=c(21,21,21,21,21,21), labels = c("Control [C]", "Heated [A]", "Bleached [B/C]", "Bleached + Heated [B]", "Negative Control [D]", "Negative Control + Heated [E]"), name = "Treatment")+
  geom_pointrange(aes(ymin=mean-se,ymax=mean+se))+
  scale_color_manual(values=cost.col.line, labels = c("Control [C]", "Heated [A]", "Bleached [B/C]", "Bleached + Heated [B]", "Negative Control [D]", "Negative Control + Heated [E]"), name = "Treatment")+
  scale_fill_manual(values=cost.col.fill, labels = c("Control [C]", "Heated [A]", "Bleached [B/C]", "Bleached + Heated [B]", "Negative Control [D]", "Negative Control + Heated [E]"), name = "Treatment")+
  theme_classic()+
  theme(legend.position="right")+
  ylab(label="Concentration (cells/µL)")+
  xlab(label="Time (hours)")
fig2B
#NEED TO ADD FIGURE OUTPUT

#write_csv(FCM.summ, file = file.path(dirOutput, "fig2B_summ_data.csv"))
#write_csv(FCM.dat, file = file.path(dirOutput, "fig2B_data.csv"))


# Stats

tmp <- FCM.dat %>% filter(Timepoint_char == "T24")

#check distribution
hist(tmp$Concentration) #not super normal
hist(sqrt(tmp$Concentration)) #more normal


#ANOVA

FCM.dat %>%
  mutate("Conc_sqrt" = sqrt(Concentration)) %>%
  anova_test(Conc_sqrt ~ Treatment * as.factor(Timepoint))

FCM.dat %>%
  mutate("Conc_sqrt" = sqrt(Concentration)) %>%
  tukey_hsd(Conc_sqrt ~ Treatment * as.factor(Timepoint)) %>% view()

#or the same but as written in the paper
tmp <- FCM.dat %>%
  mutate("Conc_sqrt" = sqrt(Concentration))
tmp$Treatment <- factor(tmp$Treatment, levels = c("Control", "Heated", "Bleached","Bleached + Heated","Negative Control", "Negative Control + Heated"))

tmp2<- aov(Conc_sqrt ~ Treatment * as.factor(Timepoint), tmp)
anova(tmp2)
TukeyHSD(tmp2, "Treatment")
TukeyHSD(tmp2, "as.factor(Timepoint)")
TukeyHSD(tmp2, "Treatment:as.factor(Timepoint)")

#only 0 and 24 hours
#ANOVA
FCM.dat %>%
  filter(Timepoint_char == "T0" | Timepoint_char == "T24")  %>%
  group_by(Timepoint_char, Treatment) %>%
  summarise(mean(Concentration), sd(Concentration), )

FCM.dat %>%
  filter(Timepoint_char == "T0" | Timepoint_char == "T24")  %>%
  anova_test(Concentration ~ Treatment * as.factor(Timepoint) )

FCM.dat %>%
  filter(Timepoint_char == "T0" | Timepoint_char == "T24")%>%
  anova_test(Concentration ~ as.factor(Timepoint) )


#tukey HSD
FCM.dat %>%
  filter(Timepoint_char == "T0" | Timepoint_char == "T24")  %>%
  tukey_hsd(Concentration ~Treatment * as.factor(Timepoint))

# now without watercontrols
FCM.dat %>%
  filter(Timepoint_char == "T0" | Timepoint_char == "T24")  %>%
  filter(!is.na(Bleaching_Status)) %>%
  anova_test(Concentration ~ Treatment * as.factor(Timepoint))

FCM.dat %>%
  filter(Timepoint_char == "T0" | Timepoint_char == "T24")  %>%
  filter(!is.na(Bleaching_Status)) %>%
  anova_test(Concentration ~ Treatment)

#tukey HSD
FCM.dat %>%
  filter(Timepoint_char == "T0" | Timepoint_char == "T24")  %>%
  filter(!is.na(Bleaching_Status)) %>%
  tukey_hsd(Concentration ~ Treatment )


###############


#ANOVA sqrt

FCM.dat %>%
  filter(Timepoint_char == "T0" | Timepoint_char == "T24")%>%
  mutate("Conc_sqrt" = sqrt(Concentration)) %>%
  anova_test(Conc_sqrt ~ Treatment)
# p =0.002

#tukey HSD
FCM.dat %>%
  filter(Timepoint_char == "T0" | Timepoint_char == "T24")%>%
  mutate("Conc_sqrt" = sqrt(Concentration)) %>%
  tukey_hsd(Conc_sqrt ~ Treatment )


# now without watercontrols
FCM.dat %>%
  filter(Timepoint_char == "T0" | Timepoint_char == "T24")%>%
  filter(!is.na(Bleaching_Status)) %>%
  mutate("Conc_sqrt" = sqrt(Concentration)) %>% view()
  anova_test(Conc_sqrt ~ Treatment)
# p = 0.075

#tukey HSD
FCM.dat %>%
  filter(Timepoint_char == "T0" | Timepoint_char == "T24")%>%
  filter(!is.na(Bleaching_Status)) %>%
  mutate("Conc_sqrt" = sqrt(Concentration)) %>%
  tukey_hsd(Conc_sqrt ~ Treatment )



# continue with specific growth rate

FCM.dat_growth <- FCM.dat %>%
  filter (Timepoint == 0 | Timepoint == 24) %>%
  select (Treatment, Timepoint_char, Bottle_NR, Concentration) %>%
  pivot_wider(values_from = Concentration, names_from = Timepoint_char) %>%
  mutate (log10_T0 = log10(T0), log10_T24 = log10(T24)) %>%
  mutate (Specific_Growth_rate = ((log10_T24 - log10_T0) *2.303)/24)

FCM.dat_growth$Treatment <- factor(FCM.dat_growth$Treatment, levels = c("Control", "Heated", "Bleached","Bleached + Heated","Negative Control", "Negative Control + Heated"))

# in paper:
FCM.dat_growth %>%
  anova_test(Specific_Growth_rate ~ Treatment )

FCM.dat_growth %>%
  tukey_hsd(Specific_Growth_rate ~ Treatment )

tmp<- aov(Specific_Growth_rate ~ Treatment, data = FCM.dat_growth)
tmp<-glht(tmp, linfct = mcp(Treatment = "Tukey"))
tmp<-tidy(cld(tmp, decreasing = TRUE)) %>%
  mutate(letters= str_to_upper(letters))


figS3<-ggplot(FCM.dat_growth, aes(Treatment, Specific_Growth_rate, color=Treatment, fill=Treatment))+
  stat_boxplot(geom = 'errorbar', size = 2.5)+
  geom_boxplot(linewidth = 2)+
  # geom_point(size = 3)+
  scale_color_manual(values=cost.col.line )+
  scale_fill_manual(values=cost.col.fill)+
  theme_classic()+
  ylab ("Specific Growth Rate (log10 cells per hour)")+
  xlab("")+
  geom_text(data=tmp, aes(x=Treatment, y=0.118, label=letters), color="black", size=5)
#NEED TO EXPORT GRAPH

write_csv(FCM.dat_growth, file.path(dirOutput, "figS3_data.csv"))



