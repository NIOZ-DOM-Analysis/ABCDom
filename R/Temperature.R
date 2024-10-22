# temperature data
#make temperature overview to show bleaching period
#data comes from LTER site.
library(tidyverse)


temp <- list.files(file.path(dirRAW, "Temperature"), pattern = "Daily_average", full.names = TRUE)
temp.data <- lapply(temp, read_csv)

# lets first select the right period (1 jan '18 - 31 dec '19)

temp.data[[1]] <- temp.data[[1]] %>%
  dplyr::filter(time_local > '2018-01-01 00:00:00') %>%
  dplyr::filter(time_local < '2020-01-01 00:00:00') %>%
  dplyr::select(time_local, contains("Temperature")) %>%
  dplyr::select(time_local,!contains("data_coverage"))

temp.data[[1]][temp.data[[1]] == -999] <- NA

temp.data[[2]] <- temp.data[[2]] %>%
  dplyr::filter(time_local > '2018-01-01 00:00:00') %>%
  dplyr::filter(time_local < '2020-01-01 00:00:00') %>%
  dplyr::select(time_local, contains("Temperature")) %>%
  dplyr::select(time_local,!contains("data_coverage"))

temp.data[[2]][temp.data[[2]] == -999] <- NA

temp.data[[3]] <- temp.data[[3]] %>%
  dplyr::filter(time_local > '2018-01-01 00:00:00') %>%
  dplyr::filter(time_local < '2020-01-01 00:00:00') %>%
  dplyr::select(time_local, contains("Temperature")) %>%
  dplyr::select(time_local,!contains("data_coverage"))

temp.data[[3]][temp.data[[3]] == -999] <- NA

temp<-full_join(temp.data[[1]], temp.data[[2]], by = "time_local", suffix = c("_FOR01", ""))
temp<-full_join(temp, temp.data[[3]], by = "time_local", suffix = c("_FOR4", "_FOR5"))

# lets make the average temp of all
tmp<-temp
tmp$mean_temp<- apply(temp[,-1],1, mean,na.rm = TRUE)
tmp$sd_temp<- apply(temp[,-1],1, sd,na.rm = TRUE)

threshold <- 29 # set colour-transition threshold
fieldwork.start<-as.POSIXct("2019-05-04 00:00:00", tz = "UTC")
fieldwork.end<-as.POSIXct("2019-05-25 00:00:00", tz = "UTC")
bleaching.start<-as.POSIXct("2019-04-01 00:00:00", tz = "UTC")
bleaching.end<-as.POSIXct("2019-05-30 00:00:00", tz = "UTC")


Fig1B<-ggplot(tmp, aes(x=time_local, y=mean_temp))+
  geom_line(y = threshold, color = "orange", size = 2.3)+
  annotate("rect", xmin = fieldwork.start, xmax = fieldwork.end, ymin = c(25), ymax = c(30.5), alpha= 0.2, color = "purple", fill = "purple")+
  geom_ribbon(aes(ymin = mean_temp -sd_temp, ymax = mean_temp + sd_temp), fill = "lightblue", alpha = 0.8)+
  geom_line(color = "black", linewidth=3)+
  geom_line(data = subset(tmp, mean_temp >= 28.9 & time_local > bleaching.start & time_local < bleaching.end),
            color = "red", linewidth = 3)+
  theme_bw(base_size = 26)+
  scale_x_datetime(name = "Date", date_breaks = "2 months", date_labels = "%b %Y")+
  scale_y_continuous(name = "mean seawater \n Temperature °C ")
Fig1B
#ggsave("Fig1B_temperature2018-2019.png", path = dirFigs, units = "in", dpi = 320, width = 20, height = 8.5)

#write_csv(tmp, file.path(dirOutput, "fig1B_data.csv"))
