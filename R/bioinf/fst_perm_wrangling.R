####FST Permutation Wrangling####
library(dplyr)
library(qvalue)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
setwd("/Users/Moritz/Documents/Academic/RSMAS/PhD/SF16_GBS/Plots/")
threeggcols <- c("#F8766D","#00BA38","#619CFF") # red,green,blue in that order
twoggcols <- c("#F8766D","#00BFC4") #coral and blue in that order

####SIMPLE FDR APPLICATION (NOT OPTIMAL!!)###
#The issue is that the minimum p-value achieved through permutations is pval = 1/N_PERM.
#This number is usually not very low, hence multiple test corrections (especially for 100,000 SNPs)
#will readily throw out these values. VERY CONSERVATIVE!

#create vector of your permutation sets
comp_mh <- c("BaFT_vs_P1FT","BaFT_vs_P2FT","BaFT_vs_P3FT")
comp_sf <- c("BaSU_vs_BaFT","P1SU_vs_P1FT","P2SU_vs_P2FT","P3SU_vs_P3FT")

#create list to store data and for outliers
stat_mh <-list()
stat_sf <- list()
outliers_mh <- list()
outliers_sf <- list()

#Read in the pvalues from the permutation analysis
for(i in comp_mh) {
  
  stat_mh[[i]] <- read.delim(paste("/Users/Moritz/Fuse_vols/STORAGE/SF16_GBS/perms/perm_thin/", i, "/pval/master.pval", sep=""), col.names = c("Chrom","Pos","He","FST","perms_in_He_range","perms_with_higher_Fst","p_val"))

  stat_mh[[i]] <- mutate(stat_mh[[i]], coord=paste(Chrom,Pos,sep="_"))
  
  #stat_mh[[i]] <- mutate(stat_mh[[i]], q_val=qvalue(stat_mh[[i]]$p_val)$qvalues)
  
  stat_mh[[i]] <- mutate(stat_mh[[i]], p_adj=p.adjust(stat_mh[[i]]$p_val, method="BH"))
  
  stat_mh[[i]] <- mutate(stat_mh[[i]], p_bonf=p.adjust(stat_mh[[i]]$p_val, method="bonferroni"))
  
  #stat_mh[[i]] <- inner_join(stat_mh[[i]],common_dnw4741_mae84616,by=c("CHROM","POS"))
  
  outliers[[i]] <- filter(stat_mh[[i]], p_val<=0.01)
  
}

for(i in comp_sf) {
  
  stat_sf[[i]] <- read.delim(paste("/Users/Moritz/Fuse_vols/STORAGE/SF16_GBS/perms/perm_thin/", i, "/pval/master.pval", sep=""), col.names = c("Chrom","Pos","He","FST","perms_in_He_range","perms_with_higher_Fst","p_val"))
  
  stat_sf[[i]] <- mutate(stat_sf[[i]], coord=paste(Chrom,Pos,sep="_"))
  
  #stat_sf[[i]] <- mutate(stat_sf[[i]], q_val=qvalue(stat_sf[[i]]$p_val)$qvalues)
  
  stat_sf[[i]] <- mutate(stat_sf[[i]], p_adj=p.adjust(stat_sf[[i]]$p_val, method="BH"))
  
  stat_sf[[i]] <- mutate(stat_sf[[i]], p_bonf=p.adjust(stat_sf[[i]]$p_val, method="bonferroni"))
  
  #stat_sf[[i]] <- inner_join(stat_sf[[i]],common_dnw4741_mae84616,by=c("CHROM","POS"))
  
  outliers[[i]] <- filter(stat_sf[[i]], p_val<=0.01)
  
}


#Save outliers for later analysis
sf_outliers <- bind_rows("Basin"=outliers[[1]], "Pond1"=outliers[[2]], "Pond2"=outliers[[3]], "Pond3"=outliers[[4]], .id="Habitat")
sf_outliers$Habitat <- as.factor(sf_outliers$Habitat)

bp_outliers <- bind_rows("Spring"=outliers[[5]], "Fall"=outliers[[6]], .id="Season")
bp_outliers$Season <- as.factor(bp_outliers$Season)

#After generating outlier sets, can now compare if some outliers are found in more than one set
#Use package VennDiagram. Since this takes a list of vectors as input you need to select the "Position" column from each dataframe in the outliers list.

venn_sf <- lapply(outliers[c(1:4)],function (x) as.vector(x$coord))
names(venn_sf) <- levels(sf_outliers$Habitat)

#venn_bp <- lapply(outliers[c(5:6)],function (x) as.vector(x$coord))

venn_overlap <- list("Basin vs. Ponds - Spring" = subset(bp_outliers, Season=="Spring")$coord, "Basin vs. Ponds - Fall" = subset(bp_outliers, Season=="Fall")$coord, "Allele frequency shifted"=sf_outliers$coord)



venn.diagram(venn_sf, filename = "test.png", fill=c(threeggcols[2],threeggcols[3],threeggcols[1],twoggcols[1]), col=NA, alpha=0.5, height=1800, width=1800)#, cat.pos=c(340,20,180), cat.dist=c(.07,.07,.05), margin=.1)



#Make a plot of Fst vs p-value. First make a dataframe from the list

stat_mh.df <- bind_rows("Basin_vs_Pond1"=stat_mh[[1]], "Basin_vs_Pond2"=stat_mh[[2]], "Basin_vs_Pond3"=stat_mh[[3]], .id="Habitat")
stat_mh.df$Habitat <- as.factor(stat_mh.df$Habitat)

stat_sf.df <- bind_rows("Basin"=stat_sf[[1]], "Pond1"=stat_sf[[2]], "Pond2"=stat_sf[[3]], "Pond3"=stat_sf[[4]], .id="Habitat")
stat_sf.df$Habitat <- as.factor(stat_sf.df$Habitat)

#########???????????/###############

test_sf <- subset(stat_sf.df, Habitat=="Pond1")
test_mh <- subset(stat_mh.df, Habitat=="Basin_vs_Pond1")

test <- full_join(test_sf,test_mh, by="coord")

ggplot(test, aes(x=-log10(p_val.x), y=-log10(p_val.y)))+
  geom_point()+
  geom_smooth(method = "lm")

##############??????????###############


png("sf_fst_pval.png", width = 2400, height = 2400, res = 300)

ggplot(sf_stat)+
  geom_point(aes(x=FST, y=-log10(p_val)), alpha=0.8)+
  geom_point(data=subset(sf_stat, p_val<=0.01), aes(x=FST, y=-log10(p_val)), col=twoggcols[2], alpha=0.8, size=2)+
  geom_point(data=subset(sf_stat, p_adj<=0.1), aes(x=FST, y=-log10(p_val)), col=twoggcols[1], alpha=0.8, size=2)+
  #geom_point(data=subset(sf_stat, p_bonf<=0.01), aes(x=FST, y=-log10(p_val)), col=threeggcols[1], alpha=0.8, size=2)+
  #geom_hline(yintercept = -log10(0.01), col=twoggcols[2], lty="twodash")+
  labs(y=expression("-log("*italic(p)~value*")"), x=expression(italic(F[ST])))+
  theme_bw()+ 
  theme(text = element_text(size=24))+
  facet_wrap(~Habitat, ncol=2)

dev.off()

#Check the pvalue distribution

png("sf_pval_dist.png", width = 2400, height = 2400, res = 300)

ggplot(sf_stat)+
  geom_histogram(aes(x=p_val), bins=50)+
  labs(x=expression(italic(p)~value), y="Count")+
  theme_bw()+
  theme(text = element_text(size=18))+
  facet_wrap(~Habitat, ncol=2)

dev.off()

#If the p value distributions look good, can now FDR correct and plot q-values

png("sf_fst_qval.png", width = 2400, height = 2400, res = 300)

ggplot(sf_stat)+
  geom_point(aes(x=FST, y=-log10(q_val)), alpha=0.8)+
  geom_point(data=subset(sf_stat, q_val<=0.1), aes(x=FST, y=-log10(q_val)), col=twoggcols[2], alpha=0.8)+
  geom_hline(yintercept = -log10(0.1), col=twoggcols[2], lty="twodash")+
  labs(y=expression("-log("*italic(q)~value*")"), x=expression(italic(F[ST])))+
  theme_bw()+
  theme(text = element_text(size=18))+
  facet_wrap(~Habitat, ncol=2)

dev.off()

#And a plot of BH-corrected p-values

png("sf_fst_padj.png", width = 2400, height = 2400, res = 300)

ggplot(sf_stat)+
  geom_point(aes(x=FST, y=-log10(p_adj)), alpha=0.8)+
  geom_point(data=subset(sf_stat, p_adj<=0.1), aes(x=FST, y=-log10(p_adj)), col=twoggcols[2], alpha=0.8)+
  geom_hline(yintercept = -log10(0.1), col=twoggcols[2], lty="twodash")+
  labs(y=expression("-log("*italic(p)~adjusted*")"), x=expression(italic(F[ST])))+
  theme_bw()+
  theme(text = element_text(size=18))+
  facet_wrap(~Habitat, ncol=2)

dev.off()

#Finally a LOSITAN style plot using the He values

png("sf_losi.png", width = 2400, height = 2400, res = 300)

ggplot(sf_stat)+
  geom_point(aes(x=He, y=FST), alpha=0.6)+
  geom_point(data=subset(sf_stat, q_val<=0.1), aes(x=He, y=FST), col=twoggcols[2], alpha=0.8)+
  labs(y=expression(italic(F[ST])), x=expression(italic(H[E])))+
  theme_bw()+
  theme(text = element_text(size=18))+
  facet_wrap(~Habitat, ncol=2)

dev.off()

##########################
###ACCURATE CUMULATIVE DISTRIBUTION FUNCTION (PNORM or PGAMMA)###
#Fit a CDF to the permuted Fst values (Difficult!). Like this we get a more accurate idea of low p-values
#i.e. we are no longer limited by the minimum p-value = 1/N_PERM

fst_perm <- read.delim(pipe("head -n 50001 /Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/perm/BaSU_vs_BaFT/fst/chunk_008.vcf.weir.fst"))

test <- subset(fst_perm, POS==2345134)$WEIR_AND_COCKERHAM_FST
test[which(test<0)] <- 0


fit <- fitdistr(test, "f", start=list(df1=1, df2=1))
hist(rgamma(10000, fit$estimate[1], fit$estimate[2]) , col = "red")
hist(test, binwidth)

qqplot(test,rgamma(length(test), fit$estimate[1], fit$estimate[2]))

pgamma(0.1, fit$estimate[1],fit$estimate[2], lower.tail = F)

gg <- as.data.frame(test)

ggplot(stat$P1SU_vs_P1FT)+
  geom_histogram(aes(x=stat$P1SU_vs_P1FT$p_val), binwidth = 0.01, alpha=0.7)+
  geom_histogram(aes(x=rexp(length(test), fit$estimate[1])), binwidth = 0.01, fill="red", alpha=0.5)
  #geom_density(aes(x=test), adjust=8)+
  #geom_density(aes(x=rpois(length(test), fit$estimate[1])), col="red")
  

