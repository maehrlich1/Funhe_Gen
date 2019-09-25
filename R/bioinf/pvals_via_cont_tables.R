#Get pvalues with a Barnard's Exact Test (or Chi2 test of equal proportions MUCH QUICKER!)
#Like this we can get p-values for the simulations too, then compare p-value distributions!
#Note that the simulations and these tests are VERY similar!
####Simulating AF changes due to random death and sampling
library(abind)
library(reshape2)
library(dplyr)
library(ggplot2)
library(Exact)
setwd("/Users/Moritz/Documents/Academic/RSMAS/PhD/SF16_GBS/Plots/Barnards Test/")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols=gg_color_hue(2)

#Populations
populations=c("Basin","Pond1","Pond2","Pond3")

#Read in Alelle count data
rawAC <- list(read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/BaSU.frq.count", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_count","alt_count")),
              read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P1SU.frq.count", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_count","alt_count")),
              read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P2SU.frq.count", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_count","alt_count")),
              read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P3SU.frq.count", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_count","alt_count")),
              read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/BaFT.frq.count", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_count","alt_count")),
              read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P1FT.frq.count", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_count","alt_count")),
              read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P2FT.frq.count", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_count","alt_count")),
              read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P3FT.frq.count", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_count","alt_count"))
)

#Generate array of AC data, dimensions: (SNP:Allele:Population:Season)

empAC <- abind(abind(rawAC[[1]][,5:6],
                     rawAC[[2]][,5:6],
                     rawAC[[3]][,5:6],
                     rawAC[[4]][,5:6],
                     along=3),
               abind(rawAC[[5]][,5:6],
                     rawAC[[6]][,5:6],
                     rawAC[[7]][,5:6],
                     rawAC[[8]][,5:6],
                     along=3),
               along = 4)

#Calculate per SNP p-values using Barnard's test. Fisher test does not work well since it is Conditional! i.e. row and column margins of the contingency table are taken as fixed!
ptm <- proc.time()
#Generate empty matrix to hold pvalues
pvals.barnard <- matrix(NA, nrow=snpcount, ncol=numpop)
#Loop over each SNP and population and skip SNPs where column margins are zero
for(pop in 1:numpop){
  for(snp in 1:snpcount){
    
    cont.tbl <- rbind(empAC[snp,,pop,2],empAC[snp,,pop,1])
    
    if(min(colSums(cont.tbl))==0){
      pvals.barnard[snp,pop] <- NA
      next
    }
    
    pvals.barnard[snp,pop] <- exact.test(rbind(empAC[snp,,pop,2],empAC[snp,,pop,1]), model="Binomial", to.plot=F)$p.value
    
  }
}

proc.time() - ptm

#Export p-values for later use
colnames(pvals.barnard) <- populations
write.table(pvals.barnard, "Barnard_pvals.txt", quote=F, sep="\t", row.names = F)

#Plot Barnards p-values as a function of dAF
colnames(pvals.barnard) <- populations
pvals.barnard.melt <- melt(pvals.barnard, varnames = c("SNP","Population"), value.name = "pval")
BpvalsempdAF <- full_join(pvals.barnard.melt, empdAF.melt, by=c("SNP","Population"))

png("dAFvspvalBarnard.png", width = 2400, height = 2400, res = 300)
ggplot(BpvalsempdAF)+
  geom_point(aes(x=dAF, y=-log10(pval)), alpha=0.8)+
  geom_point(data=subset(BpvalsempdAF, pval<=0.01), aes(x=dAF, y=-log10(pval)), col=cols[1], alpha=0.8, size=2)+
  geom_vline(xintercept = 0, col="darkgray", lty="twodash")+
  scale_x_continuous(breaks=seq(-0.4,0.4,0.2))+
  labs(y=expression(-log[10]*"("*italic(p)~value*")"), x=expression(italic(Delta*AF)))+
  theme_bw()+ 
  theme(text = element_text(size=24))+
  facet_wrap(~Population, ncol=2)
dev.off()

#Plotting pvalue distribution
png("pvalBarnard_dist.png", width = 2400, height = 2400, res = 300)
ggplot(pvals.barnard.melt)+
  geom_histogram(aes(x=pval, fill=Population), color="darkgrey", binwidth=0.05, boundary=0.05)+
  scale_y_continuous(expand = c(0,0), limits = c(0,800))+
  labs(x=expression(italic(p)~value), y="Count")+
  theme_bw()+
  scale_fill_manual(values=gg_color_hue(numpop))+
  theme(text = element_text(size=24), legend.position = "none")+
  facet_wrap(~Population, ncol=2)
dev.off()
