########Plotting VCF File INFO field data for hard filtering############
library(ggplot2)
library(dplyr)
source("~/Documents/Coding/R/General Scripts/multiplot_fn.R")
myPng <- function(..., width=6, height=6, res=300) {png(..., width=width*res, height=height*res, res=res)}

########SNPS##########
###Import data###
#can downsample data beforehand to make read-in and plotting quicker. in the end we care about the distribution which shouldnt change too much with this many data points
data <- read.delim("~/NEC_WORK/variants/merged/master_snps_raw_info_red.table", header=T)

###Plot data###

QD_plot <- ggplot(data, aes(x=QD)) +
  geom_density(color="red", fill="red", alpha=.5)+
  geom_vline(xintercept=quantile(data$QD, 0.05, na.rm = T), linetype="twodash")+
  geom_vline(xintercept=7)+
  labs(x="Quality by Depth", y="Density")

FS_plot <- ggplot(data, aes(x=FS)) +
  geom_density(color="blue", fill="blue", alpha=.5) +
  scale_x_log10(breaks=c(1,10,20,50,100))+
  geom_vline(xintercept=quantile(data$FS, 0.95, na.rm = T), linetype="twodash")+
  geom_vline(xintercept=32)+
  labs(x="Fisher Strand Bias", y="Density")

SOR_plot <- ggplot(data, aes(x=SOR)) +
  geom_density(color="darkgreen", fill="darkgreen", alpha=.5, bw=.1) +
  geom_vline(xintercept=quantile(data$SOR, 0.95, na.rm = T), linetype="twodash")+
  geom_vline(xintercept=2)+
  xlim(0,6)+
  labs(x="Strand Odds Ratio", y="")

MQ_plot <- ggplot(data, aes(x=MQ)) +
  geom_density(color="darkblue", fill="darkblue", alpha=.5, bw=.2) +
  geom_vline(xintercept=quantile(data$MQ, 0.05, na.rm = T), linetype="twodash")+
  geom_vline(xintercept=58)+
  xlim(0,75)+
  labs(x="RMS Mapping Quality", y="")

MQRankSum_plot <- ggplot(data, aes(x=MQRankSum)) +
  geom_density(color="orange", fill="orange", alpha=.5) +
  geom_vline(xintercept=quantile(data$MQRankSum, 0.05, na.rm = T), linetype="twodash")+
  geom_vline(xintercept=-0.05)+
  xlim(-1.5,1.5)+
  labs(x="Mapping Quality Rank Sum", y="")

ReadPosRankSum_plot <- ggplot(data, aes(x=ReadPosRankSum)) +
  geom_density(color="purple", fill="purple", alpha=.5, bw=.075) +
  geom_vline(xintercept=c(quantile(data$ReadPosRankSum, 0.05, na.rm = T),quantile(data$ReadPosRankSum, 0.95, na.rm = T)), linetype="twodash")+
  geom_vline(xintercept=c(-1.2,1.75))+
  xlim(-3,3)+
  labs(x="Read Position Rank Sum", y="")

myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/snp_stats/raw_snps_info.png", res=300, width = 18, height = 12)

multiplot(QD_plot, FS_plot, SOR_plot, MQ_plot, MQRankSum_plot, ReadPosRankSum_plot, cols=3)

dev.off()

#######INDELS########

data <- read.table("~/NEC_WORK/variants/merged/master_raw_indels_info.table", header=T)

str(data)

###Reduce dataset to make plotting quicker###

data_red <- data %>% sample_frac(0.1)

###Plot data###

QD_plot <- ggplot(data, aes(x=QD)) +
  geom_density(color="red", fill="red", alpha=.3, bw=.75)+
  geom_vline(xintercept=quantile(data$QD, 0.05, na.rm = T), linetype="twodash")+
  geom_vline(xintercept=7)+
  xlab("Quality by Depth")

FS_plot <- ggplot(data, aes(x=FS)) +
  geom_density(color="blue", fill="blue", alpha=.3) +
  scale_x_log10(breaks=c(10,20,25,30,35,40,45,50,60,70,80,90,100))+
  geom_vline(xintercept=quantile(data$FS, 0.95, na.rm = T), linetype="twodash")+
  geom_vline(xintercept=32)+
  xlab("Fisher's Strand Bias")

SOR_plot <- ggplot(data, aes(x=SOR)) +
  geom_density(color="darkgreen", fill="darkgreen", alpha=.3, bw=.1) +
  geom_vline(xintercept=quantile(data$SOR, 0.95, na.rm = T), linetype="twodash")+
  geom_vline(xintercept=2)+
  xlim(0,6)+
  xlab("Strand Odds Ratio")

ReadPosRankSum_plot <- ggplot(data, aes(x=ReadPosRankSum)) +
  geom_density(color="purple", fill="purple", alpha=.3, bw=.1) +
  geom_vline(xintercept=c(quantile(data$ReadPosRankSum, 0.05, na.rm = T), quantile(data$ReadPosRankSum, 0.95, na.rm = T)), linetype="twodash")+
  geom_vline(xintercept=c(-1.8,2.2))+
  xlim(-4,4)+
  xlab("Within-Read Position Rank Sum")

InbreedingCoeff_plot <- ggplot(data, aes(x=InbreedingCoeff)) +
  geom_density(color="orange", fill="orange", alpha=.3) +
  geom_vline(xintercept=quantile(data$InbreedingCoeff, 0.05, na.rm = T), linetype="twodash")+
  geom_vline(xintercept=-0.25)+
  xlim(-1,1)+
  xlab("Inbreeding Coefficient")

png(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/variant_stats/raw_indels_info.png", res=300, width = 3600, height = 1600)

multiplot(QD_plot, FS_plot, SOR_plot, ReadPosRankSum_plot, InbreedingCoeff_plot, cols=3)

dev.off()
