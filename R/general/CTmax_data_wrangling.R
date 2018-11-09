####CTmax data wrangling####
library(dplyr)
library(ggplot2)
library(RColorBrewer)
setwd("/Users/Moritz/Desktop/")
cols <- brewer.pal(8,"Paired")
threeggcols <- c("#F8766D","#00BA38","#619CFF")
twoggcols <- c("#F8766D","#00BFC4")

#Import
master <- read.delim("~/Desktop/RUMFS_May18_CTmax_data_WORKING.txt", header=T)

png("prelim_CTmax_hist.png", width = 1600, height = 1600, res = 300)

ggplot(master, aes(x=CTmax))+
  geom_histogram(fill=twoggcols[2], alpha=0.8, binwidth = 0.25)+
  theme_bw()+
  theme(text = element_text(size=18))

dev.off()

#####################

png("prelim_CTmax_team.png", width = 1600, height = 1600, res = 300)

ggplot(master, aes(x=Experimenter.s, y=CTmax))+
  geom_jitter(aes(col=Experimenter.s), position=position_jitter(0.2))+
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="black", alpha=0.5)+
  theme(legend.position="none", text = element_text(size=18))

dev.off()

#######################

png("prelim_CTmax_length.png", width = 2000, height = 1600, res = 300)

ggplot(master, aes(x=Std..Length..mm., y=CTmax))+
  geom_point(aes(col=Sex))+
  stat_smooth(method=lm, col="black", alpha=0.2)+
  theme_bw()+
  theme(legend.position=c(0.8,0.2), text = element_text(size=18))

dev.off()

cor.test(master$CTmax, master$Std..Length..mm., method="pearson")

########################

png("prelim_CTmax_recovtime.png", width = 1600, height = 2400, res = 300)

ggplot(master, aes(x=Recovery.Time*Timepoint, y=CTmax))+
  geom_boxplot(aes(group=Recovery.Time*Timepoint, fill=as.factor(Recovery.Time)))+
  theme_bw()+
  theme(legend.position="none", text = element_text(size=18))+
  labs(x="Days")+
  facet_wrap(~Recovery.Time, nrow=3)

dev.off()

############################

png("prelim_CTmax_recovtime.png", width = 1600, height = 2400, res = 300)

ggplot(master, aes(x=Recovery.Time*Timepoint, y=CTmax))+
  stat_summary(aes(y=CTmax, group=Group, col=Group), fun.y=mean, geom="line")+
  theme_bw()+
  theme(text = element_text(size=18))+
  labs(x="Days")

dev.off()