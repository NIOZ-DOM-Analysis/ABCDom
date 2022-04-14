'Diversity.R

calculates the H diversity, evenness and number of metabolites from the cleaned data.
'
########
# define some things first:
# determine the dataset you want to use
ls.diversity <- list()
ls.diversity$df.div.T0 <- df.area %>% dplyr::filter(Timepoint_char == "T0" & Treatment != "Inoculum")
ls.diversity$df.div.Tend <- df.area %>% dplyr::filter(Timepoint_char == "Tend")
ls.diversity$df.all_features<- df.area %>% dplyr::filter(Treatment != "Inoculum")
# what dataset did you use?
O <- 'area'

#today's date
D<-format(Sys.time(), '%Y%m%d')
# define the features
ls.diversity$features_T0 <- ls.diversity$df.div.T0[,M:ncol(ls.diversity$df.div.T0)]
ls.diversity$features_Tend <- ls.diversity$df.div.Tend[,M:ncol(ls.diversity$df.div.Tend)]
ls.diversity$all_features <- ls.diversity$df.all_features[,M:ncol(ls.diversity$df.all_features)]

#first T0
# define the grouping
Treatment <- factor(ls.diversity$df.div.T0$Treatment, levels = c("Non-bleached + Ambient", "Non-bleached + Heated", "Bleached + Ambient", 'Bleached + Heated', "Ambient Water Control", "Heated Water Control"))
U <- Treatment
N <- 'Treatment'


ls.diversity$H1_T0 <- diversity(ls.diversity$features_T0, index = 'shannon', MARGIN = 1, base = exp(1))
ls.diversity$S1_T0 <- specnumber(ls.diversity$features_T0)
ls.diversity$J1_T0 <- ls.diversity$H1_T0/log(ls.diversity$S1_T0)

# create new dataframe to perform the statistics
ls.diversity$H1_T0 <- bind_cols(U, ls.diversity$H1_T0) %>% dplyr::rename('U' = '...1') %>% dplyr::rename('H' = '...2')
ls.diversity$S1_T0 <- bind_cols(U, ls.diversity$S1_T0) %>% dplyr::rename('U' = '...1') %>% dplyr::rename('S' = '...2')
ls.diversity$J1_T0 <- bind_cols(U, ls.diversity$J1_T0) %>% dplyr::rename('U' = '...1') %>% dplyr::rename('J' = '...2')

cost.col.fill<-c("blue","red", "white", "white", "grey70", "grey70", "grey70")
cost.col.line<-c("dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "dodgerblue3", "firebrick3", "lightgrey")

plot1<-ggplot(ls.diversity$H1_T0, aes(x=U, y=H, group=U, color = U, fill = U))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(aes(), show.legend = TRUE, size = 1.5)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_x_discrete(name="", guide = guide_axis(angle = 45))+
  scale_fill_manual(values = cost.col.fill, name = "Treatment")+
  scale_color_manual(values = cost.col.line, name = "Treatment")+
  scale_y_continuous(name = "H' Diveristy")
ggsave("Metabolome_diversity_T0.jpeg" , path = dirFigs, dpi = 300, height=5, width=7)

plot2<-ggplot(ls.diversity$S1_T0, aes(x=U, y=S, group=U, color = U, fill = U))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(aes(), show.legend = TRUE, size = 1.5)+
  theme_classic()+
  theme(legend.position = "none")+
  scale_x_discrete(name="", guide = guide_axis(angle = 45))+
  scale_fill_manual(values = cost.col.fill, name = "Treatment")+
  scale_color_manual(values = cost.col.line, name = "Treatment")+
  scale_y_continuous(name = "Richness (S)")
ggsave("Metabolome_number_of_features_T0.jpeg" , path = dirFigs, dpi = 300, height=5, width=7)

plot3<-ggplot(ls.diversity$J1_T0, aes(x=U, y=J, group=U, color = U, fill = U))+
  stat_boxplot(geom = 'errorbar', size = 2)+
  geom_boxplot(aes(), show.legend = TRUE, size = 1.5)+
  theme_classic()+
  scale_x_discrete(name="", guide = guide_axis(angle = 45))+
  scale_fill_manual(values = cost.col.fill, name = "Treatment")+
  scale_color_manual(values = cost.col.line, name = "Treatment")+
  scale_y_continuous(name = "Pilou's evenness J=H'/log(S)")
ggsave("Metabolome_evenness_T0.jpeg" , path = dirFigs, dpi = 300, height=5, width=7)


#make it together in a plot
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend<-get_legend(plot3)
plot3 <- plot3 + theme(legend.position="none")


plot4<-plot_grid(plot1, plot2, plot3, ncol = 3, nrow = 1, rel_widths=c(1,1,1))
plot4
ggsave(paste0("Metabolome_Diversity_Richness_Evenness.png"), plot4, path =dirFigs, width = 16, height = 6)


#now do stats

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~ H ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#normality
set.seed(25)

ls.diversity$H1_T0.shap<-shapiro.test(ls.diversity$H1_T0$H)
tmp <- factor(ls.diversity$H1_T0$U)

fout3<-file(paste0(dirOutput,'\\Df1', O, N, 'Diversity stats.txt'), 'w')
writeLines(c('Normality'), fout3)
writeLines(c("shapiro test for normality"), fout3)
writeLines(c("data: ", ls.diversity$H1_T0.shap$data.name), fout3)
writeLines(c("W: ", ls.diversity$H1_T0.shap$statistic), fout3)
writeLines(c("p.value: ", ls.diversity$H1_T0.shap$p.value), fout3)
writeLines(c(" "), fout3)
writeLines(c(" "), fout3)

for (i in levels(tmp)){
  tmp1<-subset(ls.diversity$H1_T0, ls.diversity$H1_T0$U == i )
  tryCatch(tmp2<-shapiro.test(tmp1$H),
            error = function (e){
              writeLines(c("There was an error here, probably too few samples in ", i), fout3)
            })

    writeLines(c("data: ", tmp2$data.name), fout3)
    writeLines(c("subset: ", i), fout3)
    writeLines(c("W: ", tmp2$statistic), fout3)
    writeLines(c("p.value: ", tmp2$p.value), fout3)
    writeLines(c(" "), fout3)
    writeLines(c(" "), fout3)

}


#ANOVA on H'(diversity) of Coral - Algae - Interface
tmp <- tapply(ls.diversity$H1_T0$H, ls.diversity$H1_T0$U, mean)
tmp1 <- tapply(ls.diversity$H1_T0$H, ls.diversity$H1_T0$U, sd)
tmp2 <- tapply(ls.diversity$H1_T0$H, ls.diversity$H1_T0$U, length)
tmp3 <- tmp1 / sqrt(tmp2)
tmp4<-data.frame(mean = tmp, std.dev = tmp1, n = tmp2)

writeLines(c('Diversity'), fout3)
write.table(ls.diversity$H1_T0, fout3, sep='\t')
writeLines(c(' '), fout3)
writeLines(c(' '), fout3)

writeLines(c('H1 descriptives'), fout3)
write.table(tmp4, fout3, sep='\t')
writeLines(c(' '), fout3)
writeLines(c(' '), fout3)

AnovaH1 <- lm(H ~ U, data = ls.diversity$H1_T0)
AnovaH1.table<-anova(AnovaH1)  #this creates an ANOVA table
writeLines(c('H1 Anova'), fout3)
write.table(AnovaH1.table, fout3, sep='\t')
writeLines(c(' '), fout3)
writeLines(c(' '), fout3)

#R^2 of the ANOVA
AnovaSummaryH1 <- summary(AnovaH1)
AnovaH1.Rsq<-AnovaSummaryH1$r.squared
writeLines(c('R square'), fout3)
write.table(AnovaH1.Rsq, fout3, sep='\t')
writeLines(c(' '), fout3)


#TukeyHSD
AoVH1<-aov(AnovaH1)
TukH1<-(TukeyHSD(AoVH1, ordered = TRUE))
TukH1.1<-as.data.frame(TukH1$U)
writeLines(c('Tukey HSD'), fout3)
write.table(TukH1.1, fout3, sep='\t')
writeLines(c(' '), fout3)

#Kruskal Wallis rank sum-test
KrukH1<-as.list.data.frame(kruskal.test(H ~ U, data = ls.diversity$H1_T0))
writeLines(c('Kruskal Wallis rank sum-test'), fout3)
write.table(KrukH1, fout3, sep='\t')
writeLines(c(' '), fout3)

#Dunn's test on Kruskal wallis
dunnH1<-as.list.data.frame(dunn.test(ls.diversity$H1_T0$H, g=ls.diversity$H1_T0$U, table = FALSE, list = TRUE))
dunntable1<-as.data.frame(dunnH1$comparisons)
dunntable2<-as.data.frame(dunnH1$Z)
dunntable3<-as.data.frame(dunnH1$P)
dunntable4<-as.data.frame(dunnH1$P.adjusted)
dunntable<-cbind(dunntable1,dunntable2,dunntable3,dunntable4)

writeLines(c('Dunn s test'), fout3)
writeLines(c('Kruskal-Wallis chi-squared =',dunnH1$chi2), fout3)
write.table(dunntable, fout3, sep='\t')
writeLines(c(' '), fout3)

flush(fout3)
close(fout3)

#~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~ J ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(25)

ls.diversity$J1_T0.shap<-shapiro.test(ls.diversity$J1_T0$J)
tmp <- factor(ls.diversity$J1_T0$U)

fout3<-file(paste0(dirOutput,'\\Df1', O, N, 'evenness stats.txt'), 'w')
writeLines(c('Normality'), fout3)
writeLines(c("shapiro test for normality"), fout3)
writeLines(c("data: ", ls.diversity$J1_T0.shap$data.name), fout3)
writeLines(c("W: ", ls.diversity$J1_T0.shap$statistic), fout3)
writeLines(c("p.value: ", ls.diversity$J1_T0.shap$p.value), fout3)
writeLines(c(" "), fout3)
writeLines(c(" "), fout3)

for (i in levels(tmp)){
  tmp1<-subset(ls.diversity$J1_T0, ls.diversity$J1_T0$U == i )
  tryCatch(tmp2<-shapiro.test(tmp1$J),
           error = function (e){
             writeLines(c("There was an error here, probably too few samples in ", i), fout3)
           })
  writeLines(c("data: ", tmp2$data.name), fout3)
  writeLines(c("subset: ", i), fout3)
  writeLines(c("W: ", tmp2$statistic), fout3)
  writeLines(c("p.value: ", tmp2$p.value), fout3)
  writeLines(c(" "), fout3)
  writeLines(c(" "), fout3)
}


#ANOVA on J (evenness)
tmp <- tapply(ls.diversity$J1_T0$J, ls.diversity$J1_T0$U, mean)
tmp1 <- tapply(ls.diversity$J1_T0$J, ls.diversity$J1_T0$U, sd)
tmp2 <- tapply(ls.diversity$J1_T0$J, ls.diversity$J1_T0$U, length)
tmp3 <- tmp1 / sqrt(tmp2)
tmp4<-data.frame(mean = tmp, std.dev = tmp1, n = tmp2)

writeLines(c('evenness'), fout3)
write.table(ls.diversity$J1_T0, fout3, sep='\t')
writeLines(c(' '), fout3)
writeLines(c(' '), fout3)

writeLines(c('J1 descriptives'), fout3)
write.table(tmp4, fout3, sep='\t')
writeLines(c(' '), fout3)
writeLines(c(' '), fout3)

AnovaJ1 <- lm(J ~ U, data = ls.diversity$J1_T0)
AnovaJ1.table<-anova(AnovaJ1)  #this creates an ANOVA table
writeLines(c('J1 Anova'), fout3)
write.table(AnovaJ1.table, fout3, sep='\t')
writeLines(c(' '), fout3)
writeLines(c(' '), fout3)

#R^2 of the ANOVA
AnovaSummaryJ1 <- summary(AnovaJ1)
AnovaJ1.Rsq<-AnovaSummaryJ1$r.squared
writeLines(c('R square'), fout3)
write.table(AnovaJ1.Rsq, fout3, sep='\t')
writeLines(c(' '), fout3)

#TukeyHSD
AoVJ1<-aov(AnovaJ1)
TukJ1<-(TukeyHSD(AoVJ1, ordered = TRUE))
TukJ1.1<-as.data.frame(TukJ1$U)
writeLines(c('Tukey HSD'), fout3)
write.table(TukJ1.1, fout3, sep='\t')

#Kruskal Wallis rank sum-test
KrukJ1<-as.list.data.frame(kruskal.test(J ~ U, data = ls.diversity$J1_T0))
writeLines(c('Kruskal Wallis rank sum-test'), fout3)
write.table(KrukJ1, fout3, sep='\t')
writeLines(c(' '), fout3)

#Dunn's test on Kruskal wallis
dunnJ1<-as.list.data.frame(dunn.test(ls.diversity$J1_T0$J, g=ls.diversity$J1_T0$U, table = FALSE, list = TRUE))
dunntable1<-as.data.frame(dunnJ1$comparisons)
dunntable2<-as.data.frame(dunnJ1$Z)
dunntable3<-as.data.frame(dunnJ1$P)
dunntable4<-as.data.frame(dunnJ1$P.adjusted)
dunntable<-cbind(dunntable1,dunntable2,dunntable3,dunntable4)

writeLines(c('Dunn s test'), fout3)
writeLines(c('Kruskal-Wallis chi-squared =',dunnJ1$chi2), fout3)
write.table(dunntable, fout3, sep='\t')
writeLines(c(' '), fout3)

flush(fout3)
close(fout3)

#~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~ S ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(25)

ls.diversity$S1_T0.shap<-shapiro.test(ls.diversity$S1_T0$S)
tmp <- factor(ls.diversity$S1_T0$U)

fout4<-file(paste0(dirOutput,'\\Df1', O, N, 'Richness stats.txt'), 'w')
writeLines(c('Normality'), fout4)
writeLines(c("shapiro test for normality"), fout4)
writeLines(c("data: ", ls.diversity$S1_T0.shap$data.name), fout4)
writeLines(c("W: ", ls.diversity$S1_T0.shap$statistic), fout4)
writeLines(c("p.value: ", ls.diversity$S1_T0.shap$p.value), fout4)
writeLines(c(" "), fout4)
writeLines(c(" "), fout4)

for (i in levels(tmp)){
  tmp1<-subset(ls.diversity$S1_T0, ls.diversity$S1_T0$U == i )
  tryCatch(tmp2<-shapiro.test(tmp1$S),
           error = function (e){
             writeLines(c("There was an error here, probably too few samples in ", i), fout3)
           })

  writeLines(c("data: ", tmp2$data.name), fout4)
  writeLines(c("subset: ", i), fout4)
  writeLines(c("W: ", tmp2$statistic), fout4)
  writeLines(c("p.value: ", tmp2$p.value), fout4)
  writeLines(c(" "), fout4)
  writeLines(c(" "), fout4)
  }


#ANOVA on S specnumber (richness)
tmp <- tapply(ls.diversity$S1_T0$S, ls.diversity$S1_T0$U, mean)
tmp1 <- tapply(ls.diversity$S1_T0$S, ls.diversity$S1_T0$U, sd)
tmp2 <- tapply(ls.diversity$S1_T0$S, ls.diversity$S1_T0$U, length)
tmp3 <- tmp1 / sqrt(tmp2)
tmp4<-data.frame(mean = tmp, std.dev = tmp1, n = tmp2)


writeLines(c('richness'), fout3)
write.table(ls.diversity$S1_T0, fout4, sep='\t')
writeLines(c(' '), fout4)
writeLines(c(' '), fout4)

writeLines(c('S1 descriptives'), fout4)
write.table(tmp4, fout4, sep='\t')
writeLines(c(' '), fout4)
writeLines(c(' '), fout4)

AnovaS1 <- lm(S ~ ls.diversity$S1_T0$U, data = ls.diversity$S1_T0)
AnovaS1.table<-anova(AnovaS1)  #this creates an ANOVA table
writeLines(c('S1 Anova'), fout4)
write.table(AnovaS1.table, fout4, sep='\t')
writeLines(c(' '), fout4)
writeLines(c(' '), fout4)

#R^2 of the ANOVA
AnovaSummaryS1 <- summary(AnovaS1)
AnovaS1.Rsq<-AnovaSummaryS1$r.squared
writeLines(c('R square'), fout4)
write.table(AnovaS1.Rsq, fout4, sep='\t')
writeLines(c(' '), fout4)

#TukeyHSD
AoVS1<-aov(AnovaS1)
TukS1<-(TukeyHSD(AoVS1, ordered = TRUE))
TukS1.1<-as.data.frame(TukS1$U)
writeLines(c('Tukey HSD'), fout4)
write.table(TukS1.1, fout4, sep='\t')
writeLines(c(' '), fout3)

#Kruskal Wallis rank sum-test
KrukS1<-as.list.data.frame(kruskal.test(S ~ U, data = ls.diversity$S1_T0))
writeLines(c('Kruskal Wallis rank sum-test'), fout4)
write.table(KrukS1, fout4, sep='\t')
writeLines(c(' '), fout4)

#Dunn's test on Kruskal wallis
dunnS1<-as.list.data.frame(dunn.test(ls.diversity$S1_T0$S, g=ls.diversity$S1_T0$U, table = FALSE, list = TRUE))
dunntable1<-as.data.frame(dunnS1$comparisons)
dunntable2<-as.data.frame(dunnS1$Z)
dunntable3<-as.data.frame(dunnS1$P)
dunntable4<-as.data.frame(dunnS1$P.adjusted)
dunntable<-cbind(dunntable1,dunntable2,dunntable3,dunntable4)

writeLines(c('Dunn s test'), fout4)
writeLines(c('Kruskal-Wallis chi-squared =',dunnS1$chi2), fout4)
write.table(dunntable, fout4, sep='\t')
writeLines(c(' '), fout4)


flush(fout4)
close(fout4)

#~~~~~~~~~~~~~~~~Clear workspace~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


rm(tmp)
rm(tmp1)
rm(tmp2)
rm(tmp3)
rm(tmp4)
rm(AnovaH1)
rm(AnovaH1.table)
rm(AnovaH1.Rsq)
rm(AnovaSummaryH1)
rm(KrukH1)
rm(AoVH1)
rm(TukH1)
rm(TukH1.1)

rm(fout3)


rm(AnovaJ1)
rm(AnovaJ1.table)
rm(AnovaJ1.Rsq)
rm(AnovaSummaryJ1)
rm(KrukJ1)
rm(AoVJ1)
rm(TukJ1)
rm(TukJ1.1)



rm(AnovaS1)
rm(AnovaS1.table)
rm(AnovaS1.Rsq)
rm(AnovaSummaryS1)
rm(KrukS1)
rm(AoVS1)
rm(TukS1)
rm(TukS1.1)


rm(dunnH1)
rm(dunnJ1)
rm(dunnS1)
rm(dunntable)
rm(dunntable1)
rm(dunntable2)
rm(dunntable3)
rm(dunntable4)
rm(fout4)
