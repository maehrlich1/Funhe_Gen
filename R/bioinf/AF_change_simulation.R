####Simulating AF changes due to random death and sampling
library(abind)
library(reshape2)
library(dplyr)
library(ggplot2)
setwd("/Users/Moritz/Documents/Academic/RSMAS/PhD/SF16_GBS/Plots/simulations")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols=gg_color_hue(2)

####DEFINE PARAMETERS####
#Populations
populations=c("Basin","Pond1","Pond2","Pond3")
#Estimated effective population sizes
Ne=c(1300,400,400,400)
#iterations
iter=10000
#Window size for grouping simulation data. Empirical data will be compared to simulations that have a starting AF of +-windowsize/2.
winsize <- 0.01
#Read in raw data
rawAF <- list(read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/BaSU.raw.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq")),
                read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P1SU.raw.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq")),
                read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P2SU.raw.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq")),
                read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P3SU.raw.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq")),
                read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/BaFT.raw.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq")),
                read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P1FT.raw.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq")),
                read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P2FT.raw.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq")),
                read.delim("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/P3FT.raw.frq", col.names = c("CHROM","POS","N_alleles","N_chrom","ref_freq","alt_freq"))
                )

#####COMPUTATION#####
#SNPcount
snpcount=length(rawAF[[1]]$ref_freq)
#Number of populations
numpop=length(populations)
#Empirical reference AFs, array with dimensions = (SNP:Population:Season)
empAF <- abind(cbind(
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

#Population AFs under the null hypothesis are a weighted mean of the both the spring and fall sample, matrix with dimensions = (SNP:Population)
nullAF <- apply(empAF*ss, c(1,2), sum) / apply(ss, c(1,2), sum)

#Function generating HWE population of a given AF and containing Ne individuals for each snp
#Genotype encoding
#Homozygote ref. = 2, heterozygote = 1, homozygote alt. = 0
generateHWEPop <- function(p,Ne) {
  q=1-p
  c(rep(2,round((p**2)*Ne)),
    rep(1,round(2*p*q*Ne)),
    rep(0,round((q**2)*Ne)))
}

#Generate NA list to hold simulated populations with dimensions = [[population]][SNP:individuals]
simpop <- vector("list", numpop)
for (pop in 1:numpop){
  simpop[[pop]] <- lapply(nullAF[,pop], generateHWEPop, Ne=Ne[pop])
}

###Iterative subsampling of the simulated populations loop. Avoid "growing" objects in R (c, cbind, rbind). This is super slow! Instead, initialize an output matrix and fill it in!

ptm <- proc.time() #Just to keep track of computation time

#Function to sample the simulated HWE populations with a specific sample size
subsampleHWEPop <- function(simpop,ss) {
  sum(sample(simpop,ss,replace = F))/(2*ss)
}

#Generate NA array to hold simulated AFs, dimensions = (SNP:Simulation:Population:Season)
simAF <- array(NA, dim =c(snpcount,iter,numpop,2))

for (pop in 1:numpop){
  for (i in 1:iter){
    for (season in 1:2){
      
      simAF[,i,pop,season] <- mapply(subsampleHWEPop, simpop=simpop[[pop]], ss=ss[,pop,season])
      
    }
  }
}

proc.time() - ptm #Calculate elepsed computing time

#######AF ANALYSIS#########
#Calculating the AF change for both empirical and simulated data
empdAF <- empAF[,,2] - empAF[,,1]
colnames(empdAF) <- populations
simdAF <- simAF[,,,2] - simAF[,,,1]
dimnames(simdAF)[[3]] <- populations

#Plot empirical vs simulation data
#Need to put the array data into a long-format data frame for ggplot to work with.
empdAF.melt <- melt(empdAF, varnames = c("SNP","Population"), value.name = "dAF")
simdAF.melt <- melt(simdAF, varnames = c("SNP","Simulation","Population"), value.name = "dAF")

#OPTIONAL: Reduce number of simulations to experiment with plotting etc.
#simdAF.test <- simdAF[,1:100,]
#simdAF.test.melt <- melt(simdAF.test, varnames = c("SNP","Simulation","Population"), value.name = "dAF")

#Density plot of delta AFs
png("simVSempdAF.png", width = 2400, height = 2400, res = 300)

ggplot()+
  geom_line(data=simdAF.test.melt, aes(x=dAF, group=Simulation), stat='density', color="grey", lwd=0.1)+
  geom_line(data=empdAF.melt, aes(x=dAF), stat='density', color="red")+
  geom_vline(xintercept=0, color="darkgrey", lty="twodash")+
  facet_wrap(~Population)+
  scale_x_continuous(breaks=c(-0.5,-0.3,-0.1,0.1,0.3,0.5))+
  theme_bw()+
  theme(text = element_text(size=18))+
  labs(y="Density", x=expression("Temporal Change in Allele Frequency, "*italic(Delta*AF)))
  
dev.off()

#dAF as a function of initial AF
colnames(nullAF) <- populations
nullAF.melt <- melt(nullAF, varnames = c("SNP","Population"), value.name = "AF")

nullAFvsempdAF <- full_join(nullAF.melt, empdAF.melt, by=c("SNP","Population"))
nullAFvssimdAF <- full_join(nullAF.melt, simdAF.test.melt, by=c("SNP","Population"))

png("dAFvsnullAF.png", width = 2400, height = 2400, res = 300)

ggplot()+
  geom_point(data=nullAFvssimdAF, aes(x=AF, y=abs(dAF)), color="grey", size=0.01)+
  geom_point(data=nullAFvsempdAF, aes(x=AF, y=abs(dAF)), color="red", size=0.1)+
  facet_wrap(~Population)+
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1))+
  theme_bw()+
  theme(text = element_text(size=18))+
  labs(x="Null AF", y=expression("Temporal Change in Allele Frequency, "*italic(Delta*AF)))

dev.off()

####P-VALUE CALCULATION####
#In order to increase the number of simulations being used for p-value calculation, simulations with similar starting null AFs can be grouped. This is a QC to check the size of the simulation sets/bins being used to generate p-values. 
setsize <- matrix(NA, nrow=snpcount, ncol=numpop)

for(pop in 1:numpop){
  for(snp in 1:snpcount){
    
    setsize[snp,pop] <- length(which(abs(nullAF[,pop]-nullAF[snp,pop])<=winsize/2))
    
  }
}

hist(-log10(1/(setsize*iter)), main="-log10 Minimum attainable p-value")

#Define how p-values should be calculated
getpval <- function(x,set) {
  
  size <- length(set)
  p <- ifelse(x!=0,
              sum(abs(set) >= abs(x))/size,
              1)
  ifelse(p!=0,
         p,
         1/size)

}

ptm <- proc.time() #Time the pvalue calculation

#Generate empty matrix to hold pvalues, dimensions = (SNP:Population)
pvals <- matrix(NA, nrow=snpcount, ncol=numpop)
#Loop over every SNP. Get the set of simulations that have similar starting null AF and compare the empirical dAF to these.
for(pop in 1:numpop){
  for(snp in 1:snpcount){
    
    pvals[snp,pop] <- getpval(empdAF[snp,pop], simdAF[which(abs(nullAF[,pop]-nullAF[snp,pop])<=winsize/2), ,pop])
    
  }
}

proc.time() - ptm #Output the p-value calculation time.

#Export p-values for later use
colnames(pvals) <- populations
write.table(pvals, "simulation_pvals.txt", quote=F, sep="\t", row.names = F)

#Plotting pvalue distribution
colnames(pvals) <- populations
pvals.melt <- melt(pvals, varnames = c("SNP","Population"), value.name = "pval")


png("pval_dist.png", width = 2400, height = 2400, res = 300)
ggplot(pvals.melt)+
  geom_histogram(aes(x=pval, fill=Population), color="darkgrey", binwidth=0.05, boundary=0.05)+
  scale_y_continuous(expand = c(0,0), limits = c(0,800))+
  labs(x=expression(italic(p)~value), y="Count")+
  theme_bw()+
  scale_fill_manual(values=gg_color_hue(numpop))+
  theme(text = element_text(size=24), legend.position = "none")+
  facet_wrap(~Population, ncol=2)
dev.off()

#Plotting pvalue as a function of dAF
pvalsempdAF <- full_join(pvals.melt, empdAF.melt, by=c("SNP","Population"))

png("dAFvspval.png", width = 2400, height = 2400, res = 300)
ggplot(pvalsempdAF)+
  geom_point(aes(x=dAF, y=-log10(pval)), alpha=0.8)+
  geom_point(data=subset(pvalsempdAF, pval<=0.01), aes(x=dAF, y=-log10(pval)), col=cols[1], alpha=0.8, size=2)+
  geom_vline(xintercept = 0, col="darkgray", lty="twodash")+
  scale_x_continuous(breaks=seq(-0.4,0.4,0.2))+
  labs(y=expression(-log[10]*"("*italic(p)~value*")"), x=expression(italic(Delta*AF)))+
  theme_bw()+ 
  theme(text = element_text(size=24))+
  facet_wrap(~Population, ncol=2)
dev.off()
