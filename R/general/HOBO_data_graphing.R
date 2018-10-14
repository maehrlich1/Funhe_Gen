####HOBO Data Wrangling####
library(dplyr)
library(ggplot2)
library(RColorBrewer)
setwd("/Users/Moritz/Documents/Academic/RSMAS/PhD/SF16_GBS/Plots/")
cols <- brewer.pal(8,"Paired")
myPng <- function(..., width=6, height=6, res=300, ps=12) {png(..., width=width*res, height=height*res, res=res, pointsize=ps)}
threeggcols <- c("#F8766D","#00BA38","#619CFF")
twoggcols <- c("#F8766D","#00BFC4")

hobo <- read.delim("/Users/Moritz/Desktop/O2_Basin_P1_prelim.txt", header=T)

ggplot(hobo)+
  geom_line(aes(x=as.POSIXct(strptime(hobo$Time,"%d/%m/%Y %H:%M",tz="")), y=Temp..C, colour=Microhabitat), alpha=.7)+
  geom_smooth(aes(x=as.POSIXct(strptime(hobo$Time,"%d/%m/%Y %H:%M",tz="")), y=Temp..C, colour=Microhabitat), method="loess", span=0.1)+
  labs(x="Month", y="Temperature (C)")+
  theme_bw()+
  theme(text = element_text(size=24), legend.position = c(.8,.1))



ggplot(hobo)+
  geom_line(aes(x=as.POSIXct(strptime(hobo$Time,"%d/%m/%Y %H:%M",tz="")), y=DO.Adj.Conc..mg.L, colour=Microhabitat), alpha=.9)+
  #geom_smooth(aes(x=as.POSIXct(strptime(hobo$Time,"%d/%m/%Y %H:%M",tz="")), y=DO.Adj.Conc..mg.L, colour=Microhabitat), method="loess", span=0.1)+
  labs(x="Month", y="Dissolved Oxygen (mg/L)")+
  scale_x_datetime(limits=c(as.POSIXct("0018-6-20"),as.POSIXct("0018-7-09")), date_minor_breaks="1 day")+
  coord_cartesian(y=c(0,20))+
  theme_bw()+
  theme(text = element_text(size=24), legend.position = c(.2,.8))