'Stuff for figure 1'
#make temperature overview to show bleaching period
#data comes from LTER site.

temp <- list.files(paste0(dirFigs, "\\Stuff for figure 1"), pattern = ".csv", full.names = TRUE)
temp.data <- lapply(temp, read_csv)

# lets first select the right period (1 jan '18 - 31 dec '19)

temp.data[[1]] <- temp.data[[1]] %>%
  dplyr::filter(time_local > '2018-01-01 00:00:00') %>%
  dplyr::filter(time_local < '2020-01-01 00:00:00') %>%
  dplyr::select(time_local, contains("Temperature")) %>%
  dplyr::select(time_local,!contains("data_coverage"))

temp.data[[1]][temp.data[[1]] == -999] <- NA

ggplot(temp.data[[1]])+
  geom_line(aes(x=time_local, y=upper_watercolumn_temperature_sbe39_C), color = "green")+
  geom_line(aes(x=time_local, y=bottom_watercolumn_temperature_sbe39_C), color = "red")+
  geom_line(aes(x=time_local, y=middle_watercolumn_temperature_sbe39_C), color = "blue")+
  geom_line(aes(x=time_local, y=temperature_shallow_sbe37_C), color = "yellow")+
  geom_line(aes(x=time_local, y=temperature_deeper_sbe37_C), color = "purple")

# lets make the average temp of all
temp.data[[1]] <- temp.data[[1]] %>% dplyr::mutate(mean_col = rowMeans(select(temp.data[[1]],
             upper_watercolumn_temperature_sbe39_C:temperature_deeper_sbe37_C), na.rm = TRUE))
temp.data[[1]] <- temp.data[[1]] %>% mutate(Color = ifelse(mean_col > 29.5, "blue", "red"))

threshold <- 29.5 # set colour-transition threshold
fieldwork.start<-as.POSIXct("2019-05-04 00:00:00", tz = "UTC")
fieldwork.end<-as.POSIXct("2019-05-25 00:00:00", tz = "UTC")

Temp.plot<-ggplot(temp.data[[1]], aes(x=time_local, y=mean_col))+
  geom_line(color = "black", size=2)+
  geom_line(data = subset(temp.data[[1]], mean_col >= 29.5),
            color = "red", size=2)+
  theme_bw()+
  annotate("rect", xmin = fieldwork.start, xmax = fieldwork.end, ymin = c(26), ymax = c(30.5), alpha= 0.2, color = "purple", fill = "purple")+
  scale_x_datetime(name = "Date", date_breaks = "2 months", date_labels = "%b %Y")+
  scale_y_continuous(name = "mean seawater \n Temperature °C ")

Temp.plot
