####AF Plots
library(abind)
library(reshape2)
library(dplyr)
library(ggplot2)
library(cowplot)
library(scales)
library(akima)
library(lme4)
library(viridis)
setwd("/Users/Moritz/Documents/Academic/RSMAS/PhD/SF16_GBS/Plots/spatiotemporal_corr/")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols=gg_color_hue(4)

####DEFINE PARAMETERS####
#Populations
populations=c("Basin","Pond1","Pond2","Pond3")
#Read in raw data
rawAF <- list(read.delim("/Users/Moritz/Fuse_vols/STORAGE/SF16_GBS/Results/BaSU.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq")),
              read.delim("/Users/Moritz/Fuse_vols/STORAGE/SF16_GBS/Results/P1SU.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq")),
              read.delim("/Users/Moritz/Fuse_vols/STORAGE/SF16_GBS/Results/P2SU.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq")),
              read.delim("/Users/Moritz/Fuse_vols/STORAGE/SF16_GBS/Results/P3SU.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq")),
              read.delim("/Users/Moritz/Fuse_vols/STORAGE/SF16_GBS/Results/BaFT.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq")),
              read.delim("/Users/Moritz/Fuse_vols/STORAGE/SF16_GBS/Results/P1FT.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq")),
              read.delim("/Users/Moritz/Fuse_vols/STORAGE/SF16_GBS/Results/P2FT.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq")),
              read.delim("/Users/Moritz/Fuse_vols/STORAGE/SF16_GBS/Results/P3FT.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq"))
)

#####COMPUTATION#####
#SNPcount
snpcount=length(rawAF[[1]]$ref_freq)
#Number of populations
numpop=length(populations)
#Empirical reference AFs, array with dimensions = (SNP:Population:Season)
AF <- abind(cbind(
  rawAF[[1]]$ref_freq,
  rawAF[[2]]$ref_freq,
  rawAF[[3]]$ref_freq,
  rawAF[[4]]$ref_freq
), cbind(
  rawAF[[5]]$ref_freq,
  rawAF[[6]]$ref_freq,
  rawAF[[7]]$ref_freq,
  rawAF[[8]]$ref_freq
), along=3)
#Sample sizes per snp, population and season, array with dimensions = (SNP:Population:Season)
ss <- abind(cbind(
  rawAF[[1]]$N_chrom/2,
  rawAF[[2]]$N_chrom/2,
  rawAF[[3]]$N_chrom/2,
  rawAF[[4]]$N_chrom/2
), cbind(
  rawAF[[5]]$N_chrom/2,
  rawAF[[6]]$N_chrom/2,
  rawAF[[7]]$N_chrom/2,
  rawAF[[8]]$N_chrom/2
), along=3)

#Calculating temporal AF changes
dAF <- AF[,,2] - AF[,,1]
colnames(dAF) <- populations
dAF.melt <- melt(dAF, varnames = c("SNP","Population"), value.name = "dAF")

#Calculating spatial AF differences in Fall
diffAF_fall <- cbind(AF[,1,2] - AF[,2,2],
                     AF[,1,2] - AF[,3,2],
                     AF[,1,2] - AF[,4,2])
colnames(diffAF_fall) <- c("BvsP1","BvsP2","BvsP3")

diffAF_fall.melt <- melt(diffAF_fall, varnames=c("SNP","BP_pair"), value.name = "diff_AF")

#Only partially melt the dAF data to keep the Basin dAF column conserved
dAF.melt.special <- melt(as.data.frame(dAF), id="Basin", variable.name = "Pond", value.name = "dAF_pond")

#Join the temporal with the spatial differences
st <- cbind(dAF.melt.special, diffAF_fall.melt$diff_AF)
colnames(st) <- c("dAF_Basin","Pond","dAF_Pond","diffAF_fall")

#Plot the spatiotemporal data
#png("spatiotemp_point.png", width = 2400, height = 1200, res = 300)
ggplot(st %>% group_by(Pond) %>% arrange(abs(diffAF_fall)), aes(dAF_Basin, dAF_Pond))+
  geom_vline(aes(xintercept=0), lty="twodash", col="black")+
  geom_hline(aes(yintercept=0), lty="twodash", col="black")+
  geom_point(aes(color=abs(diffAF_fall), fill=diffAF_fall, alpha=abs(diffAF_fall)), shape=21, size=1.5)+
  labs(y=expression("Pond Allele Frequency Shift, "*italic(Delta*AF)),
       x=expression("Basin Allele Frequency Shift, "*italic(Delta*AF)))+
  theme_bw()+
  scale_x_continuous(breaks = c(-0.4,-0.2,0,0.2,0.4))+
  scale_y_continuous(breaks = c(-0.4,-0.2,0,0.2,0.4))+
  #coord_fixed()+
  scale_color_gradient(low="gray", high="black")+
  scale_fill_gradient2(low="red", mid="gray", high="blue",
                       limits=c(-0.25,0.25),
                       oob=squish,
                       name="Spatial AF Difference")+
  guides(color=F, size=F, alpha=F, fill=guide_colorbar(frame.colour = "black"))+
  theme(text = element_text(size=10),
        legend.direction = "horizontal",
        legend.position = "top",
        legend.title.align = 1)+
  facet_wrap(~Pond, ncol = 3)
#dev.off()

#Plot a correlation between spatial and temporal data in 2 dimensions. Have to come up with a composite AF change (Basin and Pond). Use pythagorean theorem! That is, the distance from the origin in the previous graph.

st <- st %>% mutate(comp_dAF = sqrt(dAF_Basin^2 + dAF_Pond^2))

png("spatiotemp_composite.png", width = 2400, height = 900, res = 300)
ggplot(st, aes(x=comp_dAF, y=abs(diffAF_fall)))+
  geom_point(alpha=.8)+
  geom_smooth(method="lm", se=F, color=cols[1])+
  labs(y="Spatial Allele Frequency Difference", x="Absolute Composite Allele Frequency Shift")+
  facet_wrap(~Pond)+
  theme_bw()
dev.off()

spatiotemp.model <- lmList(data=st, abs(diffAF_fall) ~ comp_dAF | Pond)

sink("spatiotemp_composite_regression.data")
summary(spatiotemp.model)
sink()

#Try to plot a raster/heatmap of the data to avoid overcrowding
#To do this we need to generate a grid of our data, i.e. interpolate between data points. Use the akima package for this.

grid.all <- data.frame(NULL)
for(i in levels(st$Pond)){
  cut <- subset(st, Pond==i)
  grid.data <- interp(x=cut$dAF_Basin, y=cut$dAF_Pond, z=cut$diffAF_fall,
                      nx = 26, ny = 26,
                      duplicate = "mean")
  grid.all <- rbind(grid.all,
                    data.frame(expand.grid(x = grid.data$x, y = grid.data$y), z = c(grid.data$z), Pond=i))
}

png("spatiotemp_heatmap.png", width = 2400, height = 1400, res = 300)
ggplot(grid.all, aes(x=x, y=y))+
  geom_raster(aes(fill=z))+
  geom_vline(aes(xintercept=0), lty="twodash", col="black")+
  geom_hline(aes(yintercept=0), lty="twodash", col="black")+
  labs(y=expression("Pond Allele Frequency Shift, "*italic(Delta*AF)),
       x=expression("Basin Allele Frequency Shift, "*italic(Delta*AF)))+
  scale_x_continuous(breaks = c(-0.4,-0.2,0,0.2,0.4))+
  scale_y_continuous(breaks = c(-0.4,-0.2,0,0.2,0.4))+
  guides(fill=guide_colorbar(frame.colour = "black"))+
  facet_wrap(~Pond, ncol = 3)+
  scale_fill_viridis(option="inferno",
                     na.value=NA,
                     name="Spatial AF Difference")+
  theme_cowplot()+
  theme(text = element_text(size=10),
        legend.direction = "horizontal",
        legend.position = "top")
dev.off()