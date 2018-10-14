####FST Permutation Wrangling####
library(dplyr)
library(qvalue)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
setwd("/Users/Moritz/Documents/Academic/RSMAS/PhD/SF16_GBS/Plots/")
cols <- brewer.pal(8,"Paired")
pops <- c(rep("Basin Fall",26), rep("Basin Spring",25), rep("Pond1 Fall",27), rep("Pond1 Spring",19), rep("Pond2 Fall",30), rep("Pond2 Spring",19), rep("Pond3 Fall",24), rep("Pond3 Spring",23))
bp_pops <- c(rep("Basin",51), rep("Pond",142))
sf_pops <- c(rep("Fall",26), rep("Spring",25), rep("Fall",27), rep("Spring",19), rep("Fall",30), rep("Spring",19), rep("Fall",24), rep("Spring",23))
myPng <- function(..., width=6, height=6, res=300, ps=12) {png(..., width=width*res, height=height*res, res=res, pointsize=ps)}
threeggcols <- c("#F8766D","#00BA38","#619CFF")
twoggcols <- c("#F8766D","#00BFC4")

####SIMPLE FDR APPLICATION (NOT OPTIMAL!!)###
#The issue is that the minimum p-value achieved through permutations is pval = 1/N_PERM.
#This number is usually not very low, hence multiple test corrections (especially for 100,000 SNPs)
#will readily throw out these values. VERY CONSERVATIVE!

#create vector of your permutation sets
comp <- c("BaSU_vs_BaFT","P1SU_vs_P1FT","P2SU_vs_P2FT","P3SU_vs_P3FT")

#create list to store data and for outliers
stat <-list()
outliers <- list()

#Read in SNP position data because the bash script drops the chrom field for some reason
#chrom <- read.delim(pipe("cut -f 1 /Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/re_Q30_DP5_bial_hwe_CR0.93.pos"), header = T)

#Read in the pvalues from the permutation analysis
for(i in comp) {
  
  #stat[[i]] <- read.delim(paste("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/perm/", i, "/pval/master.pval", sep=""), header=T)

  #stat[[i]] <- mutate(stat[[i]], ID=paste(CHROM,Pos,sep="_"))
  
  stat[[i]] <- mutate(stat[[i]], q_val=qvalue(stat[[i]]$p_val, lambda=0)$qvalues) #the lambda dictates that BH correction is used
  
  outliers[[i]] <- filter(stat[[i]], q_val<=0.05)
  
}

#Save outliers for later analysis
sf_outliers <- bind_rows("Basin"=outliers[[1]], "Pond1"=outliers[[2]], "Pond2"=outliers[[3]], "Pond3"=outliers[[4]], .id="Habitat")
sf_outliers$Habitat <- as.factor(sf_outliers$Habitat)

bp_outliers <- bind_rows("Spring"=outliers[[5]], "Fall"=outliers[[6]], .id="Season")
bp_outliers$Season <- as.factor(bp_outliers$Season)

#After generating outlier sets, can now compare if some outliers are found in more than one set
#Use package VennDiagram. Since this takes a list of vectors as input you need to select the "Position" column from each dataframe in the outliers list.

venn_sf <- lapply(outliers[c(1:4)],function (x) as.vector(x$ID))

venn_bp <- lapply(outliers[c(5:6)],function (x) as.vector(x$ID))

venn_overlap <- list("Basin vs. Ponds - Spring" = subset(bp_outliers, Season=="Spring")$ID, "Basin vs. Ponds - Fall" = subset(bp_outliers, Season=="Fall")$ID, "Allele frequency shifted"=sf_outliers$ID)

venn.diagram(venn_overlap, filename = "venn_overlap.png", fill=c(threeggcols[2],threeggcols[3],threeggcols[1]), col=NA, alpha=0.5, height=1800, width=1800, cat.pos=c(340,20,180), cat.dist=c(.07,.07,.05), margin=.1)



#Make a plot of Fst vs p-value. First make a dataframe from the list

sf_stat <- bind_rows("Basin"=stat[[1]], "Pond1"=stat[[2]], "Pond2"=stat[[3]], "Pond3"=stat[[4]], .id="Microhabitat")
sf_stat$Microhabitat <- as.factor(sf_stat$Microhabitat)

png("fst_pval_SF_allhabs.png", width = 2400, height = 1800, res = 300)

ggplot(sf_stat)+
  geom_point(aes(x=Fst, y=-log10(p_val)), alpha=0.6)+
  geom_point(data=subset(sf_stat, p_val<=0.05/28000), aes(x=Fst, y=-log10(p_val)), col=twoggcols[1], alpha=0.8)+
  geom_hline(yintercept = -log10(0.05/28000), col=twoggcols[1])+
  labs(y=expression("-log("*italic(p)~value*")"), x=expression(italic(F[ST])))+
  theme_bw()+
  theme(text = element_text(size=18))+
  facet_wrap(~Microhabitat, ncol=2)

dev.off()

#Can also plot q-values

png("fst_qval_SF _allhabs.png", width = 2400, height = 1800, res = 300)

ggplot(sf_stat)+
  geom_point(aes(x=Fst, y=-log10(q_val)), alpha=0.6)+
  geom_point(data=subset(sf_stat, q_val<=0.05), aes(x=Fst, y=-log10(q_val)), col=twoggcols[2], alpha=0.8)+
  geom_hline(yintercept = -log10(0.05), col=twoggcols[2])+
  labs(y=expression("-log("*italic(q)~value*")"), x=expression(italic(F[ST])))+
  theme_bw()+
  theme(text = element_text(size=18))+
  facet_wrap(~Microhabitat, ncol=2)

dev.off()


#Finally a LOSITAN style plot using the He values

png("he_fst_SF_allhabs.png", width = 2400, height = 1800, res = 300)

ggplot(sf_stat)+
  geom_point(aes(x=He, y=Fst), alpha=0.6)+
  geom_point(data=subset(sf_stat, q_val<=0.05), aes(x=He, y=Fst), col=twoggcols[2], alpha=0.8)+
  labs(y=expression(italic(F[ST])), x=expression(italic(H[E])))+
  theme_bw()+
  theme(text = element_text(size=18))+
  facet_wrap(~Microhabitat, ncol=2)

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
  

