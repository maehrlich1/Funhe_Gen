###Simple Fst distribution plots###
library(dplyr)
library(qvalue)
library(ggplot2)
library(RColorBrewer)
setwd("/Users/Moritz/Documents/Academic/RSMAS/PhD/SF16_GBS/Plots/")
#cols <- brewer.pal(8,"Paired")
#pops <- c(rep("Basin Fall",26), rep("Basin Spring",25), rep("Pond1 Fall",27), rep("Pond1 Spring",19), rep("Pond2 Fall",30), rep("Pond2 Spring",19), rep("Pond3 Fall",24), rep("Pond3 Spring",23))
#bp_pops <- c(rep("Basin",51), rep("Pond",142))
#sf_pops <- c(rep("Fall",26), rep("Spring",25), rep("Fall",27), rep("Spring",19), rep("Fall",30), rep("Spring",19), rep("Fall",24), rep("Spring",23))
#myPng <- function(..., width=6, height=6, res=300, ps=12) {png(..., width=width*res, height=height*res, res=res, pointsize=ps)}
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols=gg_color_hue(4)

#create vector of fst sets
comp <- c("BaSUvsBaFT","P1SUvsP1FT","P2SUvsP2FT","P3SUvsP3FT")

#create list to store data and for outliers
fst <-list()

#Read in SNP position data because the bash script drops the chrom field for some reason
#chrom <- read.delim(pipe("cut -f 1 /Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/re_Q30_DP5_bial_hwe_CR0.93.pos"), header = T)

#Read in the pvalues from the permutation analysis
for(i in comp) {
  
  fst[[i]] <- read.delim(paste("~/Fuse_vols/SCRATCH/SF16_GBS/Results/", i, "_maf5p.weir.fst", sep=""), header=T)
  
  #stat[[i]] <- mutate(stat[[i]], ID=paste(Chrom,Position,sep="_"))
  
  #stat[[i]] <- mutate(stat[[i]], q_val=qvalue(stat[[i]]$p_val)$qvalues)
  
  #outliers[[i]] <- filter(stat[[i]], p_val<=0.001)
  
}

#Convert the list to a dataframe with a new variable

fst_df <- bind_rows("Basin"=fst[[1]], "Pond1"=fst[[2]], "Pond2"=fst[[3]], "Pond3"=fst[[4]], .id="Habitat")
fst_df$Habitat <- as.factor(fst_df$Habitat)

#make density graphs

png("sf_fst_histogram_zoom_basin.png", width = 2400, height = 1600, res = 300)

ggplot(subset(fst_df, WEIR_AND_COCKERHAM_FST>0.25 & Habitat=="Basin"))+
  #geom_density(aes(x=WEIR_AND_COCKERHAM_FST, col=Habitat), alpha=0.8)+
  geom_histogram(aes(x=WEIR_AND_COCKERHAM_FST), fill=cols[1], bins=50, alpha=0.9)+
  #geom_point(data=subset(sf_stat, p_val<=0.001), aes(x=Fst, y=-log10(p_val)), col=twoggcols[2], alpha=0.8)+
  geom_vline(xintercept = 0, lty="twodash", alpha=0.8)+
  #facet_wrap(~Habitat)+
  labs(y="Count", x=expression(italic(F[ST])))+
  coord_cartesian(x=c(0.25,0.5), y=c(0,20))+
  theme_bw()+
  scale_fill_discrete(guide=F)+
  theme(legend.position = c(.8,.8), text = element_text(size=24), axis.title.x = element_blank())

dev.off()
