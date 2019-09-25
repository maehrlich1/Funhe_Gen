####HOBO Data Wrangling####
library(dplyr)
library(ggplot2)
library(RColorBrewer)
setwd("/Users/Moritz/Documents/Academic/RSMAS/PhD/SF16_GBS/Plots/HOBO/")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols=gg_color_hue(2)

hobo <- read.delim("O2_Basin_P1_prelim.txt", header=T)

png("basin_P1_temp.png", width = 3600, height = 2400, res = 300)

ggplot(hobo)+
  geom_line(aes(x=as.POSIXct(strptime(hobo$Time,"%d/%m/%Y %H:%M",tz="")), y=Temp..C, colour=Microhabitat), alpha=.9)+
  geom_smooth(aes(x=as.POSIXct(strptime(hobo$Time,"%d/%m/%Y %H:%M",tz="")), y=Temp..C, colour=Microhabitat), method="loess", span=0.1)+
  labs(x="Month", y="Temperature (C)")+
  theme_bw()+
  theme(text = element_text(size=24), legend.position = c(.8,.1))+
  scale_color_manual(values=c(twoggcols[2],twoggcols[1]))

dev.off()

#Boxplot of Temp data
png("basin_P1_temp_box.png", width = 2400, height = 2400, res = 300)
ggplot(hobo)+
  geom_boxplot(aes(x=Microhabitat, y=Temp..C, fill=Microhabitat))+
  labs(y="Temperature (C)")+
  theme_bw()+
  theme(text = element_text(size=24))+
  scale_fill_manual(values=c(cols[2],cols[1]), guide=F)
dev.off()

#Line graph over time - O2
png("Ba_P1_O2_1week.png", width = 1200, height = 600, res = 300)
ggplot(hobo)+
  geom_line(aes(x=as.POSIXct(strptime(hobo$Time,"%d/%m/%Y %H:%M",tz="")), y=DO.Adj.Conc..mg.L, colour=Microhabitat), alpha=.9)+
  #geom_smooth(aes(x=as.POSIXct(strptime(hobo$Time,"%d/%m/%Y %H:%M",tz="")), y=DO.Adj.Conc..mg.L, colour=Microhabitat), method="loess", span=0.1)+
  geom_hline(yintercept = 0.5, col= "black", lty="twodash")+
  labs(x="Day", y="Dissolved Oxygen (mg/L)")+
  scale_x_datetime(limits=c(as.POSIXct("0018-6-28"),as.POSIXct("0018-7-05")), date_minor_breaks="1 day")+
  coord_cartesian(y=c(0,15))+
  scale_color_manual(values=c(cols[2],cols[1]))+
  theme_bw()+
  theme(text = element_text(size=9),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = c(.15,.8),
        legend.title = element_blank(),
        legend.text=element_text(size=7),
        legend.key.height = unit(0.2, "cm"), 
        legend.key.width = unit(0.2, "cm"),
        legend.box.background = element_rect(colour = "black"))
dev.off()

#Line graph over time - Temp
png("Ba_P1_Temp_1week.png", width = 1200, height = 600, res = 300)
ggplot(hobo)+
  geom_line(aes(x=as.POSIXct(strptime(hobo$Time,"%d/%m/%Y %H:%M",tz="")), y=Temp..C, colour=Microhabitat), alpha=.9)+
  #geom_smooth(aes(x=as.POSIXct(strptime(hobo$Time,"%d/%m/%Y %H:%M",tz="")), y=DO.Adj.Conc..mg.L, colour=Microhabitat), method="loess", span=0.1)+
  geom_hline(yintercept = 0.5, col= "black", lty="twodash")+
  labs(x="Day", y="Temperature (C)")+
  scale_x_datetime(limits=c(as.POSIXct("0018-6-28"),as.POSIXct("0018-7-05")), date_minor_breaks="1 day")+
  coord_cartesian(y=c(17,39))+
  scale_color_manual(values=c(cols[2],cols[1]))+
  theme_bw()+
  theme(text = element_text(size=9),
        legend.position = c(.1,.86),
        legend.title = element_blank(),
        legend.text=element_text(size=7),
        legend.key.height = unit(0.2, "cm"), 
        legend.key.width = unit(0.2, "cm"),
        legend.box.background = element_rect(colour = "black"))
dev.off()
