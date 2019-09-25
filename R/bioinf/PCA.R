###Constructing and plotting a PCA based on genetic data using the SNPRelate package##
library(ggplot2)
library(SNPRelate)
library(dplyr)
setwd("/Users/Moritz/Documents/Academic/RSMAS/PhD/SF16_GBS/Plots/PCA")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#Convert VCF file to GDS format (more efficient computing)
#Once this file is generated you can skip these lines
#snpgdsVCF2GDS("/Users/Moritz/Fuse_vols/STORAGE/SF16_GBS/Results/r_Q30_DP5_bial_hwe_CR90_maf5p_thin300.vcf", "geno_thin300.gds", method="biallelic.only")
#snpgdsSummary("geno_thin300.gds")

#Open GDS file
genofile <- snpgdsOpen("geno_thin300.gds")

#Run the PCA
pca <- snpgdsPCA(genofile, autosome.only=F, num.thread=2)

#Plot the screeplot showing variance explained
#First extract the variance data from the PCA object
var_exp <- data.frame(PC=seq(1,length(pca$varprop),1), var_exp = pca$varprop)

png("var_expl_thin300.png", width = 1000, height = 1000, res = 300)
ggplot(var_exp, aes(x=PC, y=var_exp))+
  geom_bar(stat="identity")+
  labs(x="Principal Component", y="Variance Explained")+
  scale_y_continuous(limits = c(0,0.0085), labels = scales::percent_format(accuracy = 0.1))+
  theme_bw()
dev.off()

#Extract sample IDs
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

#Get population data
#pop_code <- read.delim("/Users/Moritz/Fuse_vols/STORAGE/SF16_GBS/Results/pop_post_filter.txt", header = F)
pop_code <- data.frame(ID = sample.id,
                       Habitat = c(rep("Basin", 51), rep("Pond 1", 46), rep("Pond 2", 49), rep("Pond 3", 47)),
                       Season = c(rep("Fall",26), rep("Spring",25), rep("Fall",27), rep("Spring",19), rep("Fall",30), rep("Spring",19), rep("Fall",24), rep("Spring",23)))

#Make a table of it all
pca.df <- data.frame(sample = sample.id,
                     Habitat = pop_code[,2],
                     Season = pop_code[,3],
                     PC1 = pca$eigenvect[,1],
                     PC2 = pca$eigenvect[,2],
                     PC3 = pca$eigenvect[,3],
                     PC4 = pca$eigenvect[,4],
                     PC5 = pca$eigenvect[,5],
                     stringsAsFactors = F)

#Biplot of PCs
png("PC2PC3_habseas_thin300.png", width = 1800, height = 1500, res = 300)
ggplot(pca.df, aes(x=PC2, y=PC3))+
  geom_point(aes(col = Habitat, shape = Season), alpha = 0.8)+
  stat_ellipse(aes(col = Habitat, lty = Season), alpha = 0.5)+
  labs(x="Principal Component 2", y = "Principal Component 3")+
  theme_bw()
dev.off()

#########################################################
###PCA of temporally significant SNPs only - what gives?
#Import temporal pvalue data and filter on pvalue. Generate a vector of SNP IDs that can be used to filter the GDS file
pvals <- read.delim("../pval_comparison/temporal_pvals.comp.txt")
temp_outliers_5p_ID <- which(pvals$Basin <= 0.05 | pvals$Pond1 <= 0.05 | pvals$Pond2 <= 0.05 | pvals$Pond3 <= 0.05)
temp_outliers_1p_ID <- which(pvals$Basin <= 0.01 | pvals$Pond1 <= 0.01 | pvals$Pond2 <= 0.01 | pvals$Pond3 <= 0.01)
temp_outliers_0.1p_ID <- which(pvals$Basin <= 0.001 | pvals$Pond1 <= 0.001 | pvals$Pond2 <= 0.001 | pvals$Pond3 <= 0.001)

#Run the PCA using only the temporal outlier SNPs
pca_5p <- snpgdsPCA(genofile, snp.id = temp_outliers_5p_ID, autosome.only=F, num.thread=2)
pca_1p <- snpgdsPCA(genofile, snp.id = temp_outliers_1p_ID, autosome.only=F, num.thread=2)
pca_0.1p <- snpgdsPCA(genofile, snp.id = temp_outliers_0.1p_ID, autosome.only=F, num.thread=2)

#Plot the screeplot showing variance explained
#First extract the variance data from the PCA object
var_exp_5p <- data.frame(PC=seq(1,length(pca_5p$varprop),1), var_exp = pca_5p$varprop)
var_exp_1p <- data.frame(PC=seq(1,length(pca_1p$varprop),1), var_exp = pca_1p$varprop)
var_exp_0.1p <- data.frame(PC=seq(1,length(pca_0.1p$varprop),1), var_exp = pca_0.1p$varprop)


png("var_expl_thin300_0.1p.png", width = 1000, height = 1000, res = 300)
ggplot(var_exp_0.1p, aes(x=PC, y=var_exp))+
  geom_bar(stat="identity")+
  labs(x="Principal Component", y="Variance Explained")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))+
  theme_bw()
dev.off()

#Make a table of it all
pca_5p.df <- data.frame(sample = sample.id,
                     Habitat = pop_code[,2],
                     Season = pop_code[,3],
                     PC1 = pca_5p$eigenvect[,1],
                     PC2 = pca_5p$eigenvect[,2],
                     PC3 = pca_5p$eigenvect[,3],
                     PC4 = pca_5p$eigenvect[,4],
                     PC5 = pca_5p$eigenvect[,5],
                     stringsAsFactors = F)

pca_1p.df <- data.frame(sample = sample.id,
                        Habitat = pop_code[,2],
                        Season = pop_code[,3],
                        PC1 = pca_1p$eigenvect[,1],
                        PC2 = pca_1p$eigenvect[,2],
                        PC3 = pca_1p$eigenvect[,3],
                        PC4 = pca_1p$eigenvect[,4],
                        PC5 = pca_1p$eigenvect[,5],
                        stringsAsFactors = F)

pca_0.1p.df <- data.frame(sample = sample.id,
                        Habitat = pop_code[,2],
                        Season = pop_code[,3],
                        PC1 = pca_0.1p$eigenvect[,1],
                        PC2 = pca_0.1p$eigenvect[,2],
                        PC3 = pca_0.1p$eigenvect[,3],
                        PC4 = pca_0.1p$eigenvect[,4],
                        PC5 = pca_0.1p$eigenvect[,5],
                        stringsAsFactors = F)

#Biplot of PCs
png("PC1PC2_habseas_thin300_5p.png", width = 1800, height = 1500, res = 300)
ggplot(pca_5p.df, aes(x=PC1, y=PC2))+
  geom_point(aes(col = Habitat, shape = Season), alpha = 0.8)+
  stat_ellipse(aes(col = Habitat, lty = Season), alpha = 0.5)+
  labs(x="Principal Component 1", y = "Principal Component 2")+
  theme_bw()
dev.off()

#Biplot of PCs
png("PC1PC2_habseas_thin300_1p.png", width = 1800, height = 1500, res = 300)
ggplot(pca_1p.df, aes(x=PC1, y=PC2))+
  geom_point(aes(col = Habitat, shape = Season), alpha = 0.8)+
  stat_ellipse(aes(col = Habitat, lty = Season), alpha = 0.5)+
  labs(x="Principal Component 1", y = "Principal Component 2")+
  theme_bw()
dev.off()

#Biplot of PCs
png("PC1PC2_habseas_thin300_0.1p.png", width = 1500, height = 1500, res = 300)
ggplot(pca_0.1p.df, aes(x=PC1, y=PC2))+
  geom_point(aes(col = Habitat, shape = Season), alpha = 0.8)+
  stat_ellipse(aes(col = Habitat, lty = Season), alpha = 0.5)+
  labs(x="Principal Component 1", y = "Principal Component 2")+
  theme_bw()
dev.off()

##############################################################################
##Just a test excluding the two weirdo samples
pca_test <- snpgdsPCA(genofile, sample.id = sample.id[which(sample.id != "P3SU08" & sample.id != "P2SU16")], autosome.only=F, num.thread=2)

var_exp_test <- data.frame(PC=seq(1,length(pca_test$varprop),1), var_exp = pca_test$varprop)

ggplot(var_exp_test, aes(x=PC, y=var_exp))+
  geom_bar(stat="identity")+
  labs(x="Principal Component", y="Variance Explained")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))+
  theme_bw()

pca_test.df <- data.frame(sample = sample.id[which(sample.id != "P3SU08" & sample.id != "P2SU16")],
                     Habitat = pop_code[,2][which(sample.id != "P3SU08" & sample.id != "P2SU16")],
                     Season = pop_code[,3][which(sample.id != "P3SU08" & sample.id != "P2SU16")],
                     PC1 = pca_test$eigenvect[,1],
                     PC2 = pca_test$eigenvect[,2],
                     PC3 = pca_test$eigenvect[,3],
                     PC4 = pca_test$eigenvect[,4],
                     PC5 = pca_test$eigenvect[,5],
                     stringsAsFactors = F)

#Biplot of PCs
#png("PC1PC3_habseas_thin300.png", width = 1500, height = 1500, res = 300)
ggplot(pca_test.df, aes(x=PC1, y=PC3))+
  geom_point(aes(col = Habitat, shape = Season), alpha = 0.8)+
  #stat_ellipse(aes(col = Habitat, lty = Season), alpha = 0.5)+
  labs(x="Principal Component 1", y = "Principal Component 3")+
  theme_bw()
#dev.off()
