#Compare Tajima's D of a subset of SNPs to the genome-wide distribution
library(dplyr)
library(ggplot2)
library(lme4)
setwd("/Users/Moritz/Documents/Academic/RSMAS/PhD/SF16_GBS/Plots/TajimasD")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols=gg_color_hue(4)
populations = c("Basin", "Pond1", "Pond2", "Pond3")
snpnum=10861
winsize = 300 #window size of NucDiv and TajD bins

#Get the p-values for temporally changing SNPs
pvals.all.melt <- read.delim("/Users/Moritz/Documents/Academic/RSMAS/PhD/SF16_GBS/Plots/pval_comparison/temporal_pvals.no_arl.melt.txt")

#Add locus data, create a combined p_value by taking the geometric mean and calculate the start position of the tajD bins
pos <- read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/r_Q30_DP5_bial_hwe_CR90_maf5p_thin300.pos")
pvals.all.melt <- cbind(pos,pvals.all.melt) %>% mutate(geom_mean_pval = (barnard*perms*sims)^(1/3),
                                                       lgm_pval = -log10(geom_mean_pval),
                                                       outlier = ifelse(geom_mean_pval <= 0.01, T,F),
                                                       BIN_START = POS - POS%%winsize)

#Read in the Tajima's D values for each population
tajD <- rbind(cbind(data.frame(Habitat="Basin"), read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/BaFT.300bpwin.TajimaD")),
              cbind(data.frame(Habitat="Pond1"), read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P1FT.300bpwin.TajimaD")),
              cbind(data.frame(Habitat="Pond2"), read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P2FT.300bpwin.TajimaD")),
              cbind(data.frame(Habitat="Pond3"), read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P3FT.300bpwin.TajimaD"))
              )

#Join the two data sets to merge individual SNPs with their TajD bins
tajD_pvals <- left_join(pvals.all.melt, tajD, by = c("Habitat","CHROM","BIN_START"))

#Need to construct little helper tables to plot mean TajDs
mean_tajD <- data.frame(Habitat = populations, TajimaD = c(mean(subset(tajD_pvals, Habitat=="Basin")$TajimaD, na.rm=T),
                                                           mean(subset(tajD_pvals, Habitat=="Pond1")$TajimaD, na.rm=T),
                                                           mean(subset(tajD_pvals, Habitat=="Pond2")$TajimaD, na.rm=T),
                                                           mean(subset(tajD_pvals, Habitat=="Pond3")$TajimaD, na.rm=T)
                                                           )
                        )

outlier_mean_tajD <- data.frame(Habitat = populations, TajimaD = c(mean(subset(tajD_pvals, Habitat=="Basin" & geom_mean_pval <= 0.01)$TajimaD, na.rm=T),
                                                                  mean(subset(tajD_pvals, Habitat=="Pond1" & geom_mean_pval <= 0.01)$TajimaD, na.rm=T),
                                                                  mean(subset(tajD_pvals, Habitat=="Pond2" & geom_mean_pval <= 0.01)$TajimaD, na.rm=T),
                                                                  mean(subset(tajD_pvals, Habitat=="Pond3" & geom_mean_pval <= 0.01)$TajimaD, na.rm=T)
                                                                  )
                                )

#Plot Tajima's D as a function of temporal significance
png("tajD_vs_tempPval.png", width = 2760, height = 1800, res = 300)
ggplot(tajD_pvals, aes(x = lgm_pval, y = TajimaD))+
  geom_point(alpha=0.7)+
  geom_hline(data=mean_tajD, aes(yintercept = TajimaD))+
  geom_smooth(aes(col=Habitat), method = "lm", span=0.5)+
  labs(y="Tajima's D", x=expression(-log[10]*"(Mean Temporal "*italic(p)*"-value"*")"))+
  theme_bw(base_size = 12)+
  theme(legend.background=element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  facet_wrap(~Habitat)
dev.off()

#Test for a significant slope

tajD.model <- lmList(data=tajD_pvals, TajimaD ~ lgm_pval | Habitat)
sink("LM_tajD_temporalPvals.stat")
summary(tajD.model)
sink()

#Plot distributions of outliers vs all
png("tadJ_outliersvsALL.png", width = 2760, height = 1200, res = 300)
ggplot()+
  geom_density(data = tajD_pvals, aes(x = TajimaD), fill = "gray")+
  geom_density(data = subset(tajD_pvals, outlier == T), aes(x = TajimaD), fill = cols[1], alpha=0.7)+
  geom_vline(data=mean_tajD, aes(xintercept = TajimaD))+
  geom_vline(data=outlier_mean_tajD, aes(xintercept = TajimaD), lty="twodash")+
  labs(x="Tajima's D", y="Density")+
  theme_bw(base_size = 12)+
  theme(legend.background=element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  facet_wrap(~Habitat)
dev.off()

#Testing outlier loci mean TajD vs general mean TajD

sink("tajD_outliers_vs_ALL.stat")
t.test(subset(tajD_pvals, Habitat == "Basin" & outlier == T)$TajimaD,
       subset(tajD_pvals, Habitat == "Basin")$TajimaD,
       alternative = "less")
t.test(subset(tajD_pvals, Habitat == "Pond1" & outlier == T)$TajimaD,
       subset(tajD_pvals, Habitat == "Pond1")$TajimaD,
       alternative = "less")
t.test(subset(tajD_pvals, Habitat == "Pond2" & outlier == T)$TajimaD,
       subset(tajD_pvals, Habitat == "Pond2")$TajimaD,
       alternative = "less")
t.test(subset(tajD_pvals, Habitat == "Pond3" & outlier == T)$TajimaD,
       subset(tajD_pvals, Habitat == "Pond3")$TajimaD,
       alternative = "less")
sink()

#####NUCLEOTIDE DIVERSITY#######
#Import the pi data
pi_fall <- rbind(cbind(data.frame(Habitat="Basin"), read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/BaFT.300bpwin.pi")),
              cbind(data.frame(Habitat="Pond1"), read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P1FT.300bpwin.pi")),
              cbind(data.frame(Habitat="Pond2"), read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P2FT.300bpwin.pi")),
              cbind(data.frame(Habitat="Pond3"), read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P3FT.300bpwin.pi"))
)

pi_spring <- rbind(cbind(data.frame(Habitat="Basin"), read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/BaSU.300bpwin.pi")),
              cbind(data.frame(Habitat="Pond1"), read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P1SU.300bpwin.pi")),
              cbind(data.frame(Habitat="Pond2"), read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P2SU.300bpwin.pi")),
              cbind(data.frame(Habitat="Pond3"), read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P3SU.300bpwin.pi"))
)

pi_diff <- full_join(pi_spring, pi_fall, by = c("Habitat", "CHROM", "BIN_START", "BIN_END")) %>%
            rename(N_VARIANTS_spring = N_VARIANTS.x, N_VARIANTS_fall = N_VARIANTS.y, PI_spring = PI.x, PI_fall = PI.y) %>%
            mutate(PI_diff = PI_fall - PI_spring) %>%
            mutate(BIN_START = BIN_START -1)

#Join pi and temporal pvalue data
pi_pvals <- left_join(pvals.all.melt, pi_diff, by = c("Habitat","CHROM","BIN_START"))

#Contruct helper tables
mean_pi_fall <- data.frame(Habitat = populations, PI_fall = c(mean(subset(pi_pvals, Habitat=="Basin")$PI_fall, na.rm=T),
                                                           mean(subset(pi_pvals, Habitat=="Pond1")$PI_fall, na.rm=T),
                                                           mean(subset(pi_pvals, Habitat=="Pond2")$PI_fall, na.rm=T),
                                                           mean(subset(pi_pvals, Habitat=="Pond3")$PI_fall, na.rm=T)
                                                           )
                           )

mean_pi_diff <- data.frame(Habitat = populations, PI_diff = c(mean(subset(pi_pvals, Habitat=="Basin")$PI_diff, na.rm=T),
                                                              mean(subset(pi_pvals, Habitat=="Pond1")$PI_diff, na.rm=T),
                                                              mean(subset(pi_pvals, Habitat=="Pond2")$PI_diff, na.rm=T),
                                                              mean(subset(pi_pvals, Habitat=="Pond3")$PI_diff, na.rm=T)
                                                              )
                           )

outlier_mean_pi_fall <- data.frame(Habitat = populations, PI_fall = c(mean(subset(pi_pvals, Habitat=="Basin" & geom_mean_pval <= 0.01)$PI_fall, na.rm=T),
                                                                   mean(subset(pi_pvals, Habitat=="Pond1" & geom_mean_pval <= 0.01)$PI_fall, na.rm=T),
                                                                   mean(subset(pi_pvals, Habitat=="Pond2" & geom_mean_pval <= 0.01)$PI_fall, na.rm=T),
                                                                   mean(subset(pi_pvals, Habitat=="Pond3" & geom_mean_pval <= 0.01)$PI_fall, na.rm=T)
                                                                   )
                                   )

outlier_mean_pi_diff <- data.frame(Habitat = populations, PI_diff = c(mean(subset(pi_pvals, Habitat=="Basin" & geom_mean_pval <= 0.01)$PI_diff, na.rm=T),
                                                                      mean(subset(pi_pvals, Habitat=="Pond1" & geom_mean_pval <= 0.01)$PI_diff, na.rm=T),
                                                                      mean(subset(pi_pvals, Habitat=="Pond2" & geom_mean_pval <= 0.01)$PI_diff, na.rm=T),
                                                                      mean(subset(pi_pvals, Habitat=="Pond3" & geom_mean_pval <= 0.01)$PI_diff, na.rm=T)
                                                                      )
                                   )

#Plot nucleotide diversity in fall
png("pi_vs_tempPval_fall.png", width = 2760, height = 1800, res = 300)
ggplot(pi_pvals, aes(x = lgm_pval, y = PI_fall))+
  geom_point(alpha=0.7)+
  geom_hline(data=mean_pi_fall, aes(yintercept = PI_fall))+
  geom_smooth(aes(col=Habitat), method = "lm", span=0.5)+
  labs(y="Nucleotide Diversity", x=expression(-log[10]*"(Mean Temporal "*italic(p)*"-value"*")"))+
  theme_bw(base_size = 12)+
  theme(legend.background=element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  facet_wrap(~Habitat)
dev.off()

#Test for a significant slope

pi.fall.model <- lmList(data=pi_pvals, PI_fall ~ lgm_pval | Habitat)
sink("LM_pi_fall_temporalPvals.stat")
summary(pi.fall.model)
sink()

#Plot nucleotide diversity changes
png("piDIFF_vs_tempPval.png", width = 2760, height = 1800, res = 300)
ggplot(pi_pvals, aes(x = lgm_pval, y = PI_diff))+
  geom_point(alpha=0.7)+
  geom_hline(data=mean_pi_diff, aes(yintercept = PI_diff))+
  geom_smooth(aes(col=Habitat), method = "lm", span=0.5)+
  labs(y="Change in Nucleotide Diversity", x=expression(-log[10]*"(Mean Temporal "*italic(p)*"-value"*")"))+
  theme_bw(base_size = 12)+
  theme(legend.background=element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  facet_wrap(~Habitat)
dev.off()

#Test for a significant slope

pi.diff.model <- lmList(data=pi_pvals, PI_diff ~ lgm_pval | Habitat)
sink("LM_piDIFF_temporalPvals.stat")
summary(pi.diff.model)
sink()

#Plot distributions of outliers vs all in fall
png("pi_outliersvsALL_fall.png", width = 2760, height = 1200, res = 300)
ggplot()+
  geom_density(data = pi_pvals, aes(x = PI_fall), fill = "gray")+
  geom_density(data = subset(pi_pvals, outlier == T), aes(x = PI_fall), fill = cols[1], alpha=0.7)+
  geom_vline(data=mean_pi_fall, aes(xintercept = PI_fall))+
  geom_vline(data=outlier_mean_pi_fall, aes(xintercept = PI_fall), lty="twodash")+
  labs(x="Nucleotide Diversity", y="Density")+
  theme_bw(base_size = 12)+
  theme(legend.background=element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  facet_wrap(~Habitat)
dev.off()

#Testing outlier loci mean pi in fall vs general mean pi in fall
sink("pi_outliers_vs_ALL_fall.stat")
t.test(subset(pi_pvals, Habitat == "Basin" & outlier == T)$PI_fall,
       subset(pi_pvals, Habitat == "Basin")$PI_fall,
       alternative = "less")
t.test(subset(pi_pvals, Habitat == "Pond1" & outlier == T)$PI_fall,
       subset(pi_pvals, Habitat == "Pond1")$PI_fall,
       alternative = "less")
t.test(subset(pi_pvals, Habitat == "Pond2" & outlier == T)$PI_fall,
       subset(pi_pvals, Habitat == "Pond2")$PI_fall,
       alternative = "less")
t.test(subset(pi_pvals, Habitat == "Pond3" & outlier == T)$PI_fall,
       subset(pi_pvals, Habitat == "Pond3")$PI_fall,
       alternative = "less")
sink()

#Plot distributions of temporal outliers vs all regarding change in pi
png("piDIFF_outliersvsALL.png", width = 2760, height = 1200, res = 300)
ggplot()+
  geom_density(data = pi_pvals, aes(x = PI_diff), fill = "gray")+
  geom_density(data = subset(pi_pvals, outlier == T), aes(x = PI_diff), fill = cols[1], alpha=0.7)+
  geom_vline(data=mean_pi_diff, aes(xintercept = PI_diff))+
  geom_vline(data=outlier_mean_pi_diff, aes(xintercept = PI_diff), lty="twodash")+
  labs(x="Nucleotide Diversity", y="Density")+
  theme_bw(base_size = 12)+
  theme(legend.background=element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  facet_wrap(~Habitat)
dev.off()

#Testing outlier loci mean pi in fall vs general mean pi in fall
sink("piDIFF_outliers_vs_ALL.stat")
t.test(subset(pi_pvals, Habitat == "Basin" & outlier == T)$PI_diff,
       subset(pi_pvals, Habitat == "Basin")$PI_diff,
       alternative = "less")
t.test(subset(pi_pvals, Habitat == "Pond1" & outlier == T)$PI_diff,
       subset(pi_pvals, Habitat == "Pond1")$PI_diff,
       alternative = "less")
t.test(subset(pi_pvals, Habitat == "Pond2" & outlier == T)$PI_diff,
       subset(pi_pvals, Habitat == "Pond2")$PI_diff,
       alternative = "less")
t.test(subset(pi_pvals, Habitat == "Pond3" & outlier == T)$PI_diff,
       subset(pi_pvals, Habitat == "Pond3")$PI_diff,
       alternative = "less")
sink()
