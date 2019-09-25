####Looking for concordance/discordance among populations
library(abind)
library(reshape2)
library(dplyr)
library(ggplot2)
setwd("/Users/Moritz/Documents/Academic/RSMAS/PhD/SF16_GBS/Plots/concordance")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols=gg_color_hue(2)

####DEFINE PARAMETERS####
#Populations
populations=c("Basin","Pond1","Pond2","Pond3")
#Number of populations
numpop=length(populations)
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

#####WRANGLING#####
#SNPcount
snpcount=length(rawAF[[1]]$ref_freq)
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

#Calculating temporal AF changes
dAF <- AF[,,2] - AF[,,1] #matrix with dimensions = (SNP:Population)
colnames(dAF) <- populations
dAF.melt <- melt(dAF, varnames = c("SNP","Population"), value.name = "dAF")

#Need to add p-values to these dAFs, use previously saved data
pvals.all.melt <- read.delim("/Users/Moritz/Documents/Academic/RSMAS/PhD/SF16_GBS/Plots/pval_comparison/temporal_pvals.all.melt.txt")

#Make combined composite p-values (several methods but geometric mean works well)
pvals.all.melt <- pvals.all.melt %>% mutate(geom_mean_pval = (barnard*perms*sims)^(1/3))
#Give the SNPs IDs
pvals.all.melt <- cbind(pvals.all.melt, SNP_ID = rep(1:snpcount,numpop))
#Cast the data frame and generate a cross-pond p-value and a pond-basin p-value
geom_mean_pval <- dcast(pvals.all.melt, SNP_ID~Habitat, value.var = "geom_mean_pval") %>%
                  mutate(Xpond_geom_mean_pval = Pond1*Pond2*Pond3,
                         Ba_P1_geom_mean_pval = Basin*Pond1,
                         Ba_P2_geom_mean_pval = Basin*Pond2,
                         Ba_P3_geom_mean_pval = Basin*Pond3)

#Join dAF with pvalue data
dAF_pvals <- cbind(dAF,geom_mean_pval[,6:9])

#Find number of outliers for various p-value thresholds
thresholds <- 10^-(seq(0,9,0.01))
conc_prop <- data.frame(thresholds=thresholds, num_snps=NA, num_conc=NA, conc_prop=NA)#Dimensions: (Threshold:Conc_snps)

for(th in 1:length(thresholds)){
  
  conc_prop[th,2] <- nrow(subset(dAF_pvals, Xpond_geom_mean_pval <= thresholds[th]))
  conc_prop[th,3] <- length(
                            which(
                                  abs(
                                      rowMeans(
                                              sign(
                                                    subset(dAF_pvals, Xpond_geom_mean_pval <= thresholds[th])[2:4]
                                                    )
                                              )
                                      )
                                  ==1) 
                            )
  conc_prop[th,4] <- conc_prop[th,3]/conc_prop[th,2]
      
}

###Pond concordance
#Access simulated and permuted data to have a null expectation to compare against
#simulated AF change data from simulation script




png("conc_prop_smooth.png", width = 1000, height = 1000, res = 300)
ggplot(conc_prop, aes(x=-log10(thresholds), y=conc_prop))+
  geom_hline(aes(yintercept=0.25), col="black", lty="twodash")+
  #geom_line(col=cols[1])+
  geom_smooth(span=0.2, col=cols[1])+
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1), breaks = seq(0,1,0.25))+
  scale_x_continuous(breaks = seq(0,10,2))+
  labs(x=expression(Compound~-log[10]*"("*italic(p)~value*")"~"Threshold"), y="Proportion of concordant AF shifts among ponds")+
  theme_bw()+
  theme(text = element_text(size=9))
dev.off()

#Histogram of SNP numbers
png("conc_num_hist.png", width = 1000, height = 300, res = 300)
ggplot(conc_prop, aes(x=-log10(thresholds), y=num_conc))+
  geom_histogram(stat="identity")+
  scale_x_continuous(breaks = seq(0,10,2))+
  #labs(x=expression(Compound~-log[10]*"("*italic(p)~value*")"~"Threshold"), y="No. of SNPs below threshold")+
  theme_bw()+
  theme(text = element_text(size=9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

#Basin discordance for each pond
thresholds <- 10^-(seq(0,6,0.01))
discor_prop <- data.frame(thresholds=thresholds,
                          Pond1=NA,
                          Pond2=NA,
                          Pond3=NA)

for(th in 1:length(thresholds)){

  discor_prop[th,2] <- length(
                          which(
                            abs(
                              rowMeans(
                                sign(
                                  subset(dAF_pvals, Ba_P1_geom_mean_pval <= thresholds[th])[,c(1,2)]
                                )
                              )
                            )
                          ==1) 
                        )/nrow(subset(dAF_pvals, Ba_P1_geom_mean_pval <= thresholds[th]))

  discor_prop[th,3] <- length(
                          which(
                            abs(
                              rowMeans(
                                sign(
                                  subset(dAF_pvals, Ba_P2_geom_mean_pval <= thresholds[th])[,c(1,3)]
                                )
                              )
                            )
                          ==1) 
                        )/nrow(subset(dAF_pvals, Ba_P2_geom_mean_pval <= thresholds[th]))
  
  discor_prop[th,4] <- length(
                          which(
                            abs(
                              rowMeans(
                                sign(
                                  subset(dAF_pvals, Ba_P3_geom_mean_pval <= thresholds[th])[,c(1,4)]
                                )
                              )
                            )
                          ==1)
                        )/nrow(subset(dAF_pvals, Ba_P3_geom_mean_pval <= thresholds[th]))
  
}

#melt to plot
discor_prop.melt <- melt(discor_prop, id.vars = "thresholds", variable.name = "Pond", value.name = "conc_prop")

#Plot Basin-Pond concordance
#Plot Pond concordance
png("discord_prop_smooth.png", width = 1000, height = 1000, res = 300)
ggplot(discor_prop.melt, aes(x=-log10(thresholds), y=conc_prop))+
  geom_hline(aes(yintercept=0.5), col="black", lty="twodash")+
  #geom_line(aes(col=Pond))+
  geom_smooth(aes(col=Pond), span = 0.2)+
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1), breaks = seq(0,1,0.25))+
  scale_x_continuous(breaks = seq(0,6,2))+
  labs(x=expression(Compound~-log[10]*"("*italic(p)~value*")"~"Threshold"), y="Proportion of concordant AF shifts (Pond vs Basin)")+
  theme_bw()+
  theme(text = element_text(size=9),
        legend.title = element_blank(),
        #legend.text = element_text(size=9),
        legend.position = c(0.175,0.2),
        legend.background=element_blank(),
        legend.box.background = element_rect(colour = "black"))
dev.off()
