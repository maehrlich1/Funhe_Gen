###Call libs###
setwd("/Users/Moritz/Documents/Academic/RSMAS/PhD/SF16_GBS/Plots/")
library(parallel)
library(dplyr)
library(adegenet)
library(pegas)
library(vcfR)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(Cairo)
cols <- brewer.pal(8,"Paired")
pops <- c(rep("Basin Fall",26), rep("Basin Spring",25), rep("Pond1 Fall",27), rep("Pond1 Spring",19), rep("Pond2 Fall",30), rep("Pond2 Spring",19), rep("Pond3 Fall",24), rep("Pond3 Spring",23))
bp_pops <- c(rep("Basin",51), rep("Pond",142))
sf_pops <- c(rep("Fall",26), rep("Spring",25), rep("Fall",27), rep("Spring",19), rep("Fall",30), rep("Spring",19), rep("Fall",24), rep("Spring",23))
myPng <- function(..., width=6, height=6, res=300) {png(..., width=width*res, height=height*res, res=res)}
threeggcols <- c("#00BA38","#F8766D","#619CFF")
twoggcols <- c("#00BFC4","#F8766D")

####PREPARATION OF GENLIGHT OBJECT#########

###Read VCF, convert to genlight and append pop data###

genlight <- vcfR2genlight(read.vcfR("/Users/Moritz/Fuse_vols/SCRATCH/SF16_GBS/Results/r_Q30_DP5_bial_hwe_CR90.vcf"))
pop(genlight) <- pops

save(genlight, file="r_Q30_DP5_bial_hwe_CR90.genlight")

#Alternatively if you have previously saved the genlight object

load(file="r_Q30_DP5_bial_hwe_CR90.genlight")

#if you need to subset your genlight object you can

fall_gl <- genlight[grep("*FT*",indNames(genlight))]
pop(fall_gl) <- c(rep("Basin",26), rep("Pond",81))
spring_gl <- genlight[grep("*SU*",indNames(genlight))]
pop(spring_gl) <- c(rep("Basin",25), rep("Pond",61))

#####Calculate PCA for use with DAPC#####

#nf determines the number of PCs to be retained. Ideally you want enough to explain most of the variation in as few PCs as possible but since you will not display >5-6 dimensions it is ok to go lower.
#BUT: if you want to use the results of this PCA later on in the DAPC, go for a high number of dims (>N) to not be limited later
spring_glpca <- glPca(spring_gl, nf=85, n.cores = 3, parallel = require("parallel"), useC=F)

#save(glpca, file="r_Q30_DP5_bial_hwe_CR90.glpca")

load(file="r_Q30_DP5_bial_hwe_CR90.glpca")

################################################################
#######DISCRIMINANT ANALYSIS OF PRINCIPAL COMPONENTS############
################################################################

####First identify the most likely number of clusters.

grp <- find.clusters(genlight, max.n.clust = 15)

#Limit the max number of clusters to be checked. At most there can be n-1 clusters. If populations are already known this point can be skipped though it is good practice to start with no assumptions.
#When choosing the no. of PCs to retain opt for a large number >N if computation is not an issue.
#The final choice for cluster number should be at the lowest point of the graph. A "true" value does however never exist, one should try several numbers of clusters.

#Check the independently derived cluster number and assignment against the 'real' populations

table(pop(genlight), grp$grp)

#png(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/DAPC/assignment_grid_K6.png", res=150, width = 800, height = 800)

table.value(table(pop(genlight), grp$grp), csize=3, clegend=0, col.lab=paste("Cluster", 1:6), row.lab=c("Miami","Sylt","Varna","Villefranche","Woods Hole"))

#dev.off()



#####
####Run the actual DAPC using the clusters identified previously. Here one can also use known populations!
######
#assign the pops you want to analyse
pop(genlight) <- bp_pops

dapc1 <- dapc(genlight, glPca=glpca, n.pca=36) #use grp$grp instead of the pops stored in genlight to use the inferred cluster assignment
#When choosing no. of PCs to retain look at the inflection point of the graph. Retaining too many PCs can lead to overfitting, too few to information loss. If unsure take a look at the a-score!
#When choosing the no. of discriminant functions, which should have a maximum no.of.clusters-1, retain as many as you realistically would like to display later e.g. 3 dimensions = 3 DF.
#Can also use the output of a previous PCA to save time!!!

#ASIDE: The a-score for determining the amount of PCs to be retained.

#By analysing the a-score the optimal number of PCs to be retined can be determined, similar to cluster number.

#Create a new dapc with the amount of PCs you want to test. No. of DF can be high.
temp <- dapc(spring_gl, glPca=spring_glpca, n.da=85, n.pca=85)

#Run a-score optimization

asc <- optim.a.score(temp, smart=T, n.pca=36, n.sim=100, n=36)

###
###Plot the DAPC
###

myPng("dapc_pops_DF2+3.png", res=300, width = 6, height = 6)

scatter(dapc1, xax=1, yax=4, col=cols, posi.da = "topleft", ratio.da = 0.2, 
        mstree=F, cstar = T, solid=.6, cex=3, clab=0,
        leg=T, posi.leg = "topleft", cleg=1, txt.leg = unique(genlight$pop))

dev.off()

###1 DIMENSIONAL PLOT
#If you are using only 1 dimension i.e. 1 DF then you can display indv. densities along the DF.
#get the difference in group means too
del <- abs(signif(dapc1$grp.coord[1]-dapc1$grp.coord[2],3))

#myPng("dapc_fallBP_DF1.png", res=300, width = 8, height = 6)
png("dapc_habitat_DF1.png", res=300, width = 1800, height = 1800)

scatter(dapc1,1,1, col=twoggcols, bg="white", scree.da=F, leg=TRUE, posi.leg = "topright", solid=.6, cleg = 1.2)
abline(v=dapc1$grp.coord[1], col=twoggcols[1], lty="twodash")
abline(v=dapc1$grp.coord[2], col=twoggcols[2], lty="twodash")
text(x=3, y=0.3, expression(paste(Delta~Group~Means~"=")))
text(x=3, y=0.27, del)
#title("Spring")

dev.off()

###Loading plots to determine discriminant SNPs###

png(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/DAPC/loadings_preBQSR_realpops_DF3.png", res=300, width = 2400, height = 1600)

loadingplot(dapc1$var.contr, axis=3, thr=quantile(dapc1$var.contr[,"LD3"],0.9999, type=8))
#axis determines which DF should be analyzed. thr is used to determine whether labelling should take place, here 99.9% percentile.

dev.off()

###STRUCTURE-like composition plots to determine group assignments###

compo1 <- compoplot(dapc1, col=cols)
