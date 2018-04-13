###########PCAdapt package##########
#for creating PCAs and finding outlier SNPs based on PC loadings
setwd("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/r_files")
library(pcadapt)
library(ggplot2)
library(dplyr)
library(qvalue)
library(tidyr)
cols <- c("#E41A1C","#984EA3","#FF7F00","#377EB8","#4DAF4A")
cols_south <- c("#E41A1C","#FF7F00","#377EB8")
cols_north <- c("#984EA3", "#4DAF4A")
pops <- c(rep("Miami",15), rep("Sylt",16), rep("Varna",16), rep("Villefranche",9), rep("Woods Hole",16))
pops_south <- c(rep("Miami",15), rep("Varna",16), rep("Villefranche",9))
pops_north <- c(rep("Sylt",16), rep("Woods Hole",16))
myPng <- function(..., width=6, height=6, res=300) {png(..., width=width*res, height=height*res, res=res)}
threeggcols <- c("#F8766D","#00BA38","#619CFF")
twoggcols <- c("#F8766D","#00BFC4")
set.seed(1) #the PCA starts from a random point in order to search for PCs. The result might be different depending on this starting point. By setting the seed we can ensure reproducibility.

####Convert the vcf file into something readable by PCAdapt########
#this may take some time ~15mins
read.pcadapt("/Users/Moritz/NEC_WORK/variants/merged/no_mito/master_snps_ss_GQ13_1miss_maf5p_bial.vcf", type="vcf")
#this creates a stupid tmp.pcadapt file in your working directory. BUG!! this file should work but needs to be renamed manually
file.rename("tmp.pcadapt","master_ss_GQ13_1miss_maf5p_bial.pcadapt")

#########Choosing number of PCs###########
#Here we just call the file created in the step before and run the actual PCA
pcadapt1 <- pcadapt("master_ss_GQ13_1miss_maf5p_bial.pcadapt", K=3)

#check out the screeplot of variance explained by different PCs
#plot(pcadapt2, option="screeplot")

#ggplot alternative also works but we have to do some calculations since SINGULAR VALUES are stored in the pcadapt object NOT eigenvalues
#the proportion of explained variance can be calculated by squaring the SV and dividing by the number of original variables i.e. SNPs
pcadapt1_scree <- cbind.data.frame(c(1:length(pcadapt1$singular.values)), pcadapt1$singular.values, (pcadapt1$singular.values^2)/length(pcadapt1$maf))
colnames(pcadapt1_scree) <- c("PC","sv","var_expl")

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/south_1miss/varexpl.png", res=300, width = 6, height = 6)

ggplot(pcadapt1_scree)+
  geom_bar(aes(x=PC, y=var_expl*100), stat="identity", colour=c(rep("red",2),rep("gray20",13)))+
  scale_x_continuous(breaks= c(1:length(pcadapt1_scree$PC)))+
  xlab("Principal Component")+
  ylab("Variance Explained (%)")

dev.off()

###BEFORE MOVING ON REDO THE PCA WITH THE K VALUE YOU JUST FOUND!!!! THE NUMBER OF PCs INFLUENCES THE P-VALUES STRONGLY!!!!
###ALSO RERUN THE pcadapt1_scree data WITH THE NEW K...WILL BE USED LATER
#########Plotting the PCA#########
#plot the actual pca use i and j to plot other PCs
#plot(pcadapt1, option="scores", i=1, j=2, pop=pops_south)

#Alternatively plot using ggplot. scores are found in $scores but need to be multiplied by the proportion of variance they explain in order to give a visually accurate representation
pcadapt1_scores <- cbind.data.frame(pops_north, t(t(pcadapt1$scores)*pcadapt1_scree$var_expl))
colnames(pcadapt1_scores) <- c("Population", paste("PC",seq(1:(length(pcadapt1$singular.values))), sep=""))

#PC1&2
myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/south_1miss/pca_PC1+2.png", res=300, height = 6,width = 9)

ggplot(pcadapt1_scores, aes(x=PC1, y=PC2))+
  geom_point(aes(col=Population), size=3, alpha=.8)+
  coord_fixed()+
  geom_hline(yintercept=0, alpha=.5)+
  geom_vline(xintercept=0, alpha=.5)+
  scale_color_manual(values=cols_south)+
  labs(x="Principal Component 1", y="Principal Component 2")+
  theme(legend.position = c(0.6,0.3), axis.text = element_blank(), axis.ticks = element_blank())

dev.off()

#ZOOM into northern cluster

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/pca_PC1+2_ZOOMNORTH.png", res=300, width = 12)

ggplot(subset(pcadapt1_scores, Population=="Sylt" | Population=="Woods Hole"), aes(x=PC1, y=PC2))+
  geom_point(aes(col=Population), size=3, alpha=.8)+
  coord_fixed()+
  geom_hline(yintercept=0, alpha=.5)+
  scale_color_manual(values=cols_north)+
  labs(x="Principal Component 1", y="Principal Component 2")+
  theme(legend.position = c(0.1,0.75), axis.text = element_blank(), axis.ticks = element_blank())

dev.off()

#PC2&3
myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/pca_PC2+3.png", res=300, width = 6, height = 12)

ggplot(pcadapt1_scores, aes(x=PC3, y=PC2))+
  geom_point(aes(col=Population), size=3, alpha=.8)+
  coord_fixed()+
  geom_hline(yintercept=0, alpha=.5)+
  geom_vline(xintercept=0, alpha=.5)+
  scale_color_manual(values=cols)+
  labs(x="Principal Component 3", y="Principal Component 2")+
  theme(legend.position = c(0.8,0.2), axis.text = element_blank(), axis.ticks = element_blank())

dev.off()

#ZOOM into northern cluster
myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/pca_PC2+3_ZOOMNORTH.png", res=300, width = 6, height = 12)

ggplot(subset(pcadapt1_scores, Population=="Sylt" | Population=="Woods Hole"), aes(x=PC3, y=PC2))+
  geom_point(aes(col=Population), size=3, alpha=.8)+
  coord_fixed()+
  geom_hline(yintercept=0, alpha=.5)+
  geom_vline(xintercept=0, alpha=.5)+
  scale_color_manual(values=cols_north)+
  labs(x="Principal Component 3", y="Principal Component 2")+
  theme(legend.position = c(0.2,0.2), axis.text = element_blank(), axis.ticks = element_blank())

dev.off()

#PC1&3
myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/pca_PC1+3.png", res=300, width = 12, height = 4)

ggplot(pcadapt1_scores, aes(x=PC1, y=PC3))+
  geom_point(aes(col=Population), size=3, alpha=.8)+
  coord_fixed()+
  geom_hline(yintercept=0, alpha=.5)+
  geom_vline(xintercept=0, alpha=.5)+
  scale_color_manual(values=cols)+
  labs(x="Principal Component 3", y="Principal Component 2")+
  theme(legend.position = c(0.3,0.275), axis.text = element_blank(), axis.ticks = element_blank())

dev.off()

#ZOOM into northern cluster
myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/pca_PC1+3_ZOOMNORTH.png", res=300, width = 6, height = 12)

ggplot(subset(pcadapt1_scores, Population=="Sylt" | Population=="Woods Hole"), aes(x=PC1, y=PC3))+
  geom_point(aes(col=Population), size=3, alpha=.8)+
  coord_fixed()+
  geom_hline(yintercept=0, alpha=.5)+
  geom_vline(xintercept=0, alpha=.5)+
  scale_color_manual(values=cols_north)+
  labs(x="Principal Component 3", y="Principal Component 2")+
  theme(legend.position = c(0.2,0.2), axis.text = element_blank(), axis.ticks = element_blank())

dev.off()

###can also plot a GRID MATRIX using the pairs function
#pairs(pcadapt1_scores[,2:4], col=pcadapt1_scores$Population, pch=20)

#GGally package extends ggplot2 to add a function called ggpairs() for these plot matrices which is better
library(GGally)

#3 PCs
#choose the variables with "columns" / "mapping" can be used to add aes() like in ggplot / "upper" "diag" "lower" define the plot regions, by supplying "continuous", "discrete" or "combo" the type of plot which needs to be altered can be called. the type and options can be supplied with wrap.
ggp <- ggpairs(pcadapt1_scores, columns=c("PC1","PC2","PC3"),
               mapping = aes(color=Population),
               upper = NULL,
               diag = list(continuous=wrap("densityDiag", adjust=2, alpha=.6)),
               lower = list(continuous=wrap("points", size=3, alpha=.6)))

#Unfortunately cant supply themes like in ggplot. BUT plots are saved in a matrix form e.g. ggpairs[a,b]. if we supply themes to each plot in the matrix separately it does work. use a loop.
for (row in seq_len(ggp$nrow))
  for (col in seq_len(ggp$ncol))
    ggp[row, col] <- ggp[row, col]+
                      geom_hline(yintercept=0, alpha=.5)+
                      geom_vline(xintercept=0, alpha=.5)+
                      scale_fill_manual(values = cols_north)+
                      scale_color_manual(values=cols_north)+
                      theme(axis.ticks = element_blank(), axis.text=element_blank())

#we only want 1 legend rather that 6 so we add it manually to just one plot and move it to a position we like
ggp[2,2] <- ggp[2,2] + theme(legend.position=c(1.5,1.5), legend.title = element_text(size=18), legend.text = element_text(size=14), legend.key.size = unit(0.8,"cm"))

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/north_1miss/pca_matrix_3PCs.png", res=300, width = 20, height = 10)

ggp

dev.off()

#2 PCs
#choose the variables with "columns" / "mapping" can be used to add aes() like in ggplot / "upper" "diag" "lower" define the plot regions, by supplying "continuous", "discrete" or "combo" the type of plot which needs to be altered can be called. the type and options can be supplied with wrap.
ggp <- ggpairs(pcadapt1_scores, columns=c("PC1","PC2"),
               mapping = aes(color=Population),
               upper = NULL,
               diag = list(continuous=wrap("densityDiag", adjust=2, alpha=.6)),
               lower = list(continuous=wrap("points", size=3, alpha=.6)))

#Unfortunately cant supply themes like in ggplot. BUT plots are saved in a matrix form e.g. ggpairs[a,b]. if we supply themes to each plot in the matrix separately it does work. use a loop.
for (row in seq_len(ggp$nrow))
  for (col in seq_len(ggp$ncol))
    ggp[row, col] <- ggp[row, col]+
                        geom_hline(yintercept=0, alpha=.5)+
                        geom_vline(xintercept=0, alpha=.5)+
                        scale_fill_manual(values = cols_north)+
                        scale_color_manual(values=cols_north)+
                        theme(axis.ticks = element_blank(), axis.text=element_blank())

#we only want 1 legend rather that 6 so we add it manually to just one plot and move it to a position we like
ggp[2,2] <- ggp[2,2] + theme(legend.position=c(0.5,1.5), legend.title = element_text(size=20), legend.text = element_text(size=16), legend.key.size = unit(1,"cm"))

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/north_1miss/pca_matrix_2PCs.png", res=300, width = 9, height = 11)

ggp

dev.off()

#########PROJECTING NEW/OTHER SAMPLES ONTO EXISTING PCs#############
#This works by multiplying the genotypes of the new samples by the loadings of the PCs previously found. We have to do some work beforehand though...
#NOTE: The calculations here do NOT replicate those of PCAdapt, but they are pretty close...To be consistent, use the loadings from PCAdapt but recalculate the scores of both old and new samples yourself.
#First calculate the loadings using PCAdapt and just a subset of your samples (those that should determine the PCs)
G_master_1miss <- read.table("master_ss_GQ13_1miss_maf5p_bial.pcadapt")#we need a common SNP set later on so we have to use this set and just take a subset of it
G_master_1miss_subset_north <- bind_cols(G_master_1miss[,16:31], G_master_1miss[,57:72])#this takes Sylt and Woods Hole
G_master_1miss_subset_south <- bind_cols(G_master_1miss[,1:15], G_master_1miss[,32:56])#this takes Miami, Varna and Villefranche

write.table(G_master_1miss_subset_north, "master_ss_GQ13_1miss_maf5p_bial_subsnorth.pcadapt", row.names = F, col.names = F)#need to write a file, otherwise it wont work
write.table(G_master_1miss_subset_south, "master_ss_GQ13_1miss_maf5p_bial_subssouth.pcadapt", row.names = F, col.names = F)#need to write a file, otherwise it wont work

#run the PCA
pcadapt_justnorth <- pcadapt("master_ss_GQ13_1miss_maf5p_bial_subsnorth.pcadapt", K=2)#use K=2 from before
pcadapt_justsouth <- pcadapt("master_ss_GQ13_1miss_maf5p_bial_subssouth.pcadapt", K=2)#use K=2 from before

#check out how it looks quickly
#plot(pcadapt_justnorth, option="scores", i=1, j=2, pop=pops_north)

#extract the loadings to calculate our new scores with
loadings_north <- pcadapt_justnorth$loadings
loadings_south <- pcadapt_justsouth$loadings

#and set NA loci to zero. these were loci that had a maf below 0.05 after subsetting the genotype matrix
loadings_north[which(is.na(loadings_north))] <- 0
loadings_south[which(is.na(loadings_south))] <- 0

#Back to the genotype matrix. Remember: we want to plot all to be consistent in the calculations.
#its still in pcadapt format where NAs are represented as 9s, need to change this

G_master_1miss[G_master_1miss==9] <- NA

#normalize the genotype matrix. this normalization is based on the MAF and is commonly used in population genomics (Patterson 2006).
#create our own function for this which first retrieves the alternative AFs, convert them to mafs and then uses them to normalize
geno.normalize <- function(G) {
  alt_afs = apply(G, 1, function(m) sum(m, na.rm=T)/(2*sum(is.na(m) == F)))
  mafs <- ifelse(alt_afs>0.5, 1-alt_afs, alt_afs)
  as.matrix((G-mafs)/sqrt(2*mafs*(1-mafs)))
}

G_norm <- geno.normalize(G_master_1miss)

#next we need to deal with missing data. here we use the median genotype value of that SNP to replace NAs with. This is not ideal but its the best i could do.
#again write a function that loops through every SNP and replaces any NA values it finds with the median value.
geno.impute <- function(G) {
  for(i in 1:nrow(G)) {
    G[i,is.na(G[i,])] <- median(G[i,], na.rm = T)
  }
  G
}

G_norm_imp <- geno.impute(G_norm)

#MAYBE STILL NA VALUES LEFT!
#remaining NAs are caused by invariant sites where the maf was calculated as 0. this gives a division by zero during the normalization and NA values for all samples at that locus.
#if we used the entire dataset we shouldnt have any fixed sites. nonetheless, by setting NAs to zero we cause them to have no influence on the scores.
#G_norm_imp[which(is.na(G_norm_imp))] <- 0

#Finally we can calculate the scores of our new samples by multiplying the normalized and imputed genotype matrix by the loadings. Need to transpose the genotype matrix t for it to work the right way.
new_scores_northPC <- as.data.frame(t(G_norm_imp) %*% loadings_north) %>% mutate(Population=pops)
new_scores_southPC <- as.data.frame(t(G_norm_imp) %*% loadings_south) %>% mutate(Population=pops)

colnames(new_scores_northPC) <- c(paste("PC",seq(1:(ncol(new_scores_northPC)-1)), sep=""), "Population")
colnames(new_scores_southPC) <- c(paste("PC",seq(1:(ncol(new_scores_southPC)-1)), sep=""), "Population")

#plot the newly calculated scores as projections on the PCs defined earlier
#Projection onto Northern PCs
myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/projections/all_onto_northPCs.png", res=300, width = 6, height = 8.5)

ggplot(new_scores_northPC, aes(x=-PC1, y=-PC2))+
  geom_point(aes(col=Population), size=3, alpha=.8)+
  coord_fixed()+
  geom_hline(yintercept=0, alpha=.5)+
  geom_vline(xintercept=0, alpha=.5)+
  scale_color_manual(values=cols)+
  labs(x="Principal Component 1", y="Principal Component 2")+
  theme(legend.position = c(0.5,0.2), axis.text = element_blank(), axis.ticks = element_blank())

dev.off()

##ZOOM if needed
myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/projections/all_onto_northPCs_ZOOM.png", res=300, width = 9, height = 6)

ggplot(subset(new_scores_northPC, Population=="Miami" | Population=="Varna" | Population=="Villefranche"), aes(x=PC1, y=PC2))+
  geom_point(aes(col=Population), size=3, alpha=.8)+
  scale_color_manual(values=cols_south)+
  labs(x="Principal Component 1", y="Principal Component 2")+
  coord_fixed()+
  theme(legend.position = c(0.85,0.3), axis.text = element_blank(), axis.ticks = element_blank())

dev.off()

#Projection onto Southern PCs
myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/projections/all_onto_southPCs.png", res=300, width = 8, height = 6)

ggplot(new_scores_southPC, aes(x=PC1, y=PC2))+
  geom_point(aes(col=Population), size=3, alpha=.8)+
  coord_fixed()+
  geom_hline(yintercept=0, alpha=.5)+
  geom_vline(xintercept=0, alpha=.5)+
  scale_color_manual(values=cols)+
  labs(x="Principal Component 1", y="Principal Component 2")+
  theme(legend.position = c(0.6,0.2), axis.text = element_blank(), axis.ticks = element_blank())

dev.off()

##ZOOM if needed
myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/projections/all_onto_southPCs_ZOOM.png", res=300, width = 8, height = 6)

ggplot(subset(new_scores_southPC, Population=="Sylt" | Population=="Woods Hole"), aes(x=PC1, y=PC2))+
  geom_point(aes(col=Population), size=3, alpha=.8)+
  geom_hline(yintercept=0, alpha=.5)+
  geom_vline(xintercept=0, alpha=.5)+
  scale_color_manual(values=cols_north)+
  labs(x="Principal Component 1", y="Principal Component 2")+
  coord_fixed()+
  theme(legend.position = c(0.2,0.8), axis.text = element_blank(), axis.ticks = element_blank())

dev.off()

###########################################
#######SNP loadings and statistics##########
#############################################

######check out the structure of the pcadapt object
summary(pcadapt1)

#######histogram of the p-values DO THIS FIRST!!!
#look at the shape of the histogram to check how your test went!
hist(pcadapt1$pvalues, xlab="p-values", main=NULL, breaks=50)

#import SNP location data for later comparison of outliers. PCAdapt only uses biallelic SNPs so make sure to get only this data
pos <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/master_snps_ss_GQ13_1miss_maf5p_bial.pos")

#ggplot, $pvalues. extract p-values. values of zero just means that R cannot handle such low values, they should be real though.
pcadapt1_pval <- cbind.data.frame(pos, c(1:NROW(pcadapt1$pvalues)), pcadapt1$pvalues)
colnames(pcadapt1_pval) <- c("chrom","pos","snp_id","p_val")
pcadapt1_pval <- mutate(pcadapt1_pval, log_p_val= -log(p_val))

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/pvalhist.png", res=300, width = 6, height = 6)

ggplot(pcadapt1_pval)+
  geom_density(aes(x=p_val), adjust=.25, fill="coral", colour="gray20", alpha=.8)+
  labs(x=expression(italic(p)*"-value"), y="Density")

dev.off()

#########manattan plot of p-values associated with snps. a p-value of zero means R could not handle such a low value. they should be included though
#plot(pcadapt1, option="manhattan")

#ggplot
myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/pvalmanha.png", res=300, width = 12, height = 6)

ggplot(pcadapt1_pval)+
  geom_point(aes(x=snp_id, y=log_p_val), alpha=.6, color="coral")+
  geom_hline(yintercept = -log(0.001), colour="gray20", lty="twodash")+
  labs(x="SNP ID", y=expression("-log10("*italic(p)*"-value)"))+
  annotate("text", x=2e5 ,y=-log(0.001)*3, label="p-value = 0.001", colour="gray20")

dev.off()

##########qqplot to check distribution of p-values. with the threshold command (vertical line) the top x percent of data are separated. Any point above the diagonal are smaller than the expcted p-values!
#here the pcadapt function is actually more useful but we have to get rid of zeroes first so we'll set them to na
myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/pvalqq.png", res=300, width = 6, height = 6)

plot(pcadapt1, option="qqplot", threshold=0.1)

dev.off()

###########Plot the Mahalonobis distance test statistic directly (should be chi-squared with K degrees of freedom)
#plot(pcadapt1, option="stat.distribution")

#ggplot, $chi2.stat
pcadapt1_md2 <- cbind.data.frame(c(1:NROW(pcadapt1$stat)), pcadapt1$chi2.stat)
colnames(pcadapt1_md2) <- c("snp_id", "md2_gif_adj")

myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/md2adjhist_ZOOM.png", res=300, width = 8, height = 6)

#CAREFUL: need to adjust the number of interpolation points (n) in the overlaid function to be larger than the number of bins used in the histogram!!
ggplot(pcadapt1_md2, aes(x=md2_gif_adj))+
  geom_histogram(aes(y=..density..), color="coral", fill="coral", alpha=.8, bins=7000)+
  stat_function(fun=dchisq, color="black", args=list(df=3), n = 100000, lty="twodash")+
  labs(x=expression("GIF-corrected Mahalanobis Distance,"~italic(D^2)), y="Density")+
  coord_cartesian(x=c(0,40))

dev.off()

##################Thresholding to define outliers################
#for this use the package "qvalue" which transforms p to q-values, we can add these qvalue calcs straight to the pvals
#Also calculate the -log10 of the qvals since we are interested in very small values
pcadapt1_pq <- pcadapt1_pval %>% mutate(q_val=qvalue(pcadapt1_pval$p_val)$qvalues) %>% mutate(log_q_val=-log(q_val))

#save this data
#save(pcadapt1_pq, file="/Users/Moritz/Documents/Academic/GEOMAR/Thesis/r_files/master_snps_ss_GQ13_1miss_maf5p_bial.pqvals")

#plot the qvalues and see which values lie below alpha where alpha is the false discovery rate i.e. for alpha=0.1 10% of outlier SNPs are false positive
#have to plot -log(qval) in order get outliers to be large rather than small numbers, then use -log(alpha) as the threshold 

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/qvalmanha.png", res=300, width = 8, height = 4)

ggplot(pcadapt1_pq)+
  geom_point(aes(x=snp_id, y=log_q_val), color="coral", alpha=.6)+
  geom_hline(yintercept = -log(0.01), color="gray20", lty="twodash")+
  annotate("text", label="alpha = 0.01", x=1e4, y=-log(0.01)*2)+
  labs(x="SNP ID", y=expression("-log10("*italic(q)*"-value)"))

dev.off()

myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/qvalhist_ZOOM.png", res=300, width = 8, height = 6)

ggplot(pcadapt1_pq)+
  geom_histogram(aes(x=log_q_val, y=..density..), fill="coral", color="coral", alpha=.8, binwidth=.25)+
  geom_vline(xintercept = -log(0.01), colour="gray20")+
  annotate("text", label="alpha = 0.01", x=-log(0.01)*0.95, y=0.5, angle=90)+
  labs(y="Density", x=expression("-log10("*italic(q)*"-value)"))+
  coord_cartesian(x=c(0,10))

dev.off()

#can also check which of the principal components are most correlated with the outliers
#get outliers first
#outliers_joint <- pcadapt1_pq %>% dplyr::filter(log_q_val > 175)

#################################################################################
####################SNP LOADINGS BY PRNCIPAL COMPONENT - ########################
#################################################################################
#in order to get PC-wise test statistics you can use the method="componentwise" command which will calculate p-values not based on the Mahalonobis distance but on the distribution of loadings for each PC

pcadapt2 <- pcadapt("master_ss_GQ13_1miss_maf5p_bial.pcadapt", method="componentwise", K=3)
summary(pcadapt2)

##############SNP p-values and statistics (this time component-wise)#################

#######histogram of the p-values DO THIS FIRST!!!
#look at the shape of the histogram to check how your test went!
#extract pvalues
pcadapt2_pval <- cbind.data.frame(pos, c(1:NROW(pcadapt2$pvalues)), pcadapt2$pvalues[,1], pcadapt2$pvalues[,2], pcadapt2$pvalues[,3])
colnames(pcadapt2_pval) <- c("chrom","pos","snp_id","PC1", "PC2", "PC3")
#if we want to have PCs as factors we should convert the data to long format, use tidyr
pcadapt2_pval <- gather(pcadapt2_pval,"pc","p_val",4:6)
pcadapt2_pval$pc <- as.factor(pcadapt2_pval$pc)
#calc -log10 of p-vals
pcadapt2_pval <- mutate(pcadapt2_pval, log_p_val= -log(p_val))

myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/compo_pvaldens.png", res=300, width = 12, height = 6)

ggplot(pcadapt2_pval)+
  geom_density(aes(x=p_val, fill=pc), adjust=.25, colour="gray20", alpha=.8)+
  labs(x=expression(italic(p)*"-value"), y="Density")+
  scale_fill_discrete(guide=F)+
  facet_grid(.~pc)

dev.off()

#########manattan plot of p-values associated with snps. a p-value of zero means that R could not handle such a low value, it should be included though.
#ggplot

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/compo_pvalmanha.png", res=300, width = 12, height = 18)

ggplot(pcadapt2_pval)+
  geom_point(aes(x=snp_id, y=log_p_val, color=pc), alpha=.8)+
  geom_hline(yintercept = -log(0.001), colour="gray20", lty="twodash")+
  labs(x="SNP ID", y=expression("-log10("*italic(p)*"-value)"))+
  scale_color_discrete(guide=F)+
  annotate("text", x=2e5 ,y=-log(0.001)+17, label="p-value = 0.001", colour="gray20")+
  facet_grid(pc~.)

dev.off()

##########qqplot to check distribution of p-values. with the threshold command (vertical line) the top x percent of data are separated. Any point above the diagonal are smaller than the expcted p-values!
#here the pcadapt function is actually more useful

#myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/pvalqq_PC3_ss_GQ13_1miss_maf5p.png", res=300, width = 6, height = 6)

#plot(pcadapt2, option="qqplot", threshold=0.1, K=3)

#dev.off()

###########Plot the mahalonobis distance test statistic directly (should be chi-squared with K degrees of freedom)

#ggplot, $chi2.stat
pcadapt2_md2 <- cbind.data.frame(pos, c(1:NROW(pcadapt2$chi2.stat)), pcadapt2$chi2.stat[,1], pcadapt2$chi2.stat[,2], pcadapt2$chi2.stat[,3])
colnames(pcadapt2_md2) <- c("chrom","pos","snp_id","PC1", "PC2", "PC3")
#if we want to have PCs as factors we should convert the data to long format, use tidyr
pcadapt2_md2 <- gather(pcadapt2_md2,"pc","md2_gif_adj",4:6)
pcadapt2_md2$pc <- as.factor(pcadapt2_md2$pc)

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/compo_md2adjhist.png", res=300, width = 8, height = 12)

#CAREFUL: need to adjust the number of interpolation points (n) in the overlaid function to be larger than the number of bins used in the histogram!! Also use DF=1 since we are looking at single PCs
ggplot(pcadapt2_md2, aes(x=md2_gif_adj))+
  geom_histogram(aes(y=..density.., fill=pc), color="gray20", alpha=.8, binwidth=0.25)+
  stat_function(fun=dchisq, color="black", args=list(df=1), n = 100000, lty="twodash")+
  labs(x=expression("GIF-corrected Mahalanobis Distance,"~italic(D^2)), y="Density")+
  scale_fill_discrete(guide=F)+
  facet_grid(pc~.)+
  coord_cartesian(x=c(0,15))

dev.off()

##################Thresholding to define outliers################
#for this use the package "qvalue" which transforms p to q-values, we can add these qvalue calcs straight to the pvals
#Also calculate the -log10 of the qvals since we are interested in very small values
pcadapt2_pq <- pcadapt2_pval %>% mutate(q_val=qvalue(pcadapt2_pval$p_val)$qvalues) %>% mutate(log_q_val=-log(q_val))

#save this data
#save(pcadapt2_pq, file="/Users/Moritz/Documents/Academic/GEOMAR/Thesis/r_files/master_snps_ss_GQ13_1miss_maf5p_bial.pqvals_compo")

#plot the qvalues and see which values lie below alpha where alpha is the false discovery rate i.e. for alpha=0.1 10% of outlier SNPs are false positive
#have to plot -log(qval) in order get outliers to be large rather than small numbers, then use -log(alpha) as the threshold 

myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/compo_qvalmanha.png", res=300, width = 12, height = 18)

ggplot(pcadapt2_pq)+
  geom_point(aes(x=snp_id, y=log_q_val, color=pc), alpha=.6)+
  geom_hline(yintercept = -log(0.01), color="gray20", lty="twodash")+
  annotate("text", label="alpha = 0.01", x=1e4, y=-log(0.01)*3.5)+
  labs(x="SNP ID", y=expression("-log10("*italic(q)*"-value)"))+
  scale_color_discrete(guide=F)+
  facet_grid(pc~.)

dev.off()

myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/compo_qvalhist.png", res=300, width = 8, height = 12)

ggplot(pcadapt2_pq)+
  geom_histogram(aes(x=log_q_val, y=..density.., fill=pc), alpha=.8, binwidth=.5)+
  geom_vline(xintercept = -log(0.01), colour="gray20")+
  annotate("text", label="alpha = 0.01", x=-log(0.01)*1.3, y=0.5, angle=90)+
  labs(y="Density", x=expression("-log10("*italic(q)*"-value)"))+
  scale_fill_discrete(guide=F)+
  coord_cartesian(x=c(0,15))+
  facet_grid(pc~.)

dev.off()

############Getting outliers################
#can do this by setting a hard threshold, depends on data...
#outliers <- pcadapt2_pq[which(pcadapt2_pq$q_val<0.01),]

#or inspect graphs manually and select areas of largest outliers. for this its better to plot single graphs than facets
#PC1
myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/PC1_qvalmanha.png", res=300, width = 8, height = 4)

ggplot(subset(pcadapt2_pq, pc=="PC1"))+
  geom_point(aes(x=snp_id, y=log_q_val), color=threeggcols[1], alpha=.6)+
  geom_hline(yintercept = -log(0.01), color="gray20", lty="twodash")+
  annotate("text", label="alpha = 0.01", x=1e4, y=-log(0.01)*3.5)+
  labs(x="SNP ID", y=expression("-log10("*italic(q)*"-value)"))

dev.off()

#PC2
myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/PC2_qvalmanha.png", res=300, width = 8, height = 4)

ggplot(subset(pcadapt2_pq, pc=="PC2"))+
  geom_point(aes(x=snp_id, y=log_q_val), color=threeggcols[2], alpha=.6)+
  geom_hline(yintercept = -log(0.01), color="gray20", lty="twodash")+
  annotate("text", label="alpha = 0.01", x=1e4, y=-log(0.01)*3.5)+
  labs(x="SNP ID", y=expression("-log10("*italic(q)*"-value)"))

dev.off()

#PC3
myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/1miss/PC3_qvalmanha.png", res=300, width = 8, height = 4)

ggplot(subset(pcadapt2_pq, pc=="PC3"))+
  geom_point(aes(x=snp_id, y=log_q_val), color=threeggcols[3], alpha=.6)+
  geom_hline(yintercept = -log(0.01), color="gray20", lty="twodash")+
  annotate("text", label="alpha = 0.01", x=1e4, y=-log(0.01)*2)+
  labs(x="SNP ID", y=expression("-log10("*italic(q)*"-value)"))

dev.off()

####the MASS package from pcadapt masks the dplyr filter function, so have to call it with dplyr
outliers_pc3 <- pcadapt2_pq %>% dplyr::filter(pc == "PC3" & log_q_val > 30)
outliers_pc2 <- pcadapt2_pq %>% dplyr::filter(pc == "PC2" & log_q_val > 150)
outliers_pc1 <- pcadapt2_pq %>% dplyr::filter(pc == "PC1" & log_q_val > 600)

outliers_compo <- bind_rows(outliers_pc1,outliers_pc2,outliers_pc3) %>% arrange(snp_id)

write.table(outliers_compo, sep="\t", file="/Users/Moritz/Documents/Academic/GEOMAR/Thesis/r_files/pcadapt_master_snps_ss_GQ13_1miss_maf5p_bial.outliers_compo")

