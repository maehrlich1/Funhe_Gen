###Call libs###

library(dplyr)
library(adegenet)
library(pegas)
library(vcfR)
library(ape)
library(ggplot2)
cols <- c("#E41A1C","#984EA3","#FF7F00","#377EB8","#4DAF4A")
pops <- c(rep("Miami",15), rep("Sylt",16), rep("Varna",16), rep("Villefranche",9), rep("Woods Hole",16))

####PREPARATION OF GENLIGHT OBJECT#########

###Read VCF, convert to genlight and append pop data###

genlight <- vcfR2genlight(read.vcfR("/Users/Moritz/NEC_WORK/variants/merged/master_snps_ss_GQ13_90pCall_maf5p.vcf.gz"))
pop(genlight) <- pops

##############PRINCIPAL COMPONENT ANALYSIS################

###Calculate PCA###

#nf determines the number of PCs to be retained. Ideally you want enough to explain most of the variation in as few PCs as possible but since you will not display >5-6 dimensions it is ok to go lower.
#BUT: if you want to use the results of this PCA later on in the DAPC, go for a high number of dims (>N) to not be limited later
pca1 <- glPca(genlight, nf=80)

###Plotting PCA###
#Unlike the dapc-object, the pca-object does not have a nice inherent plotting function. we must therefore extract the data from the object beforehand.

eigen <- cbind.data.frame(seq(1,NROW(pca1$eig),1),pca1$eig)
colnames(eigen) <- c("Principal Component","Eigenvalue")

pca1_simple <- cbind.data.frame(rownames(pca1$scores), pops, pca1$scores[,1], pca1$scores[,2], pca1$scores[,3], pca1$scores[,4], pca1$scores[,5])
colnames(pca1_simple) <- c("ID","Population","PC1","PC2","PC3","PC4","PC5")

#now we can plot using ggplot

png(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/eigen_ss_GQ13_90pCall_maf5p.png", res=300, width = 1600, height = 1600)

ggplot(eigen)+
    geom_bar(aes(x=`Principal Component`, y=Eigenvalue), stat="identity", colour=c(rep("red",3),rep("gray20",13)))+
    xlim(0,16)

dev.off()

png(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/pca_ss_GQ13_90pCall_maf5p_PC2+3.png", res=300, width = 1600, height = 1600)

ggplot(pca1_simple, aes(x=PC2, y=PC3))+
  geom_point(aes(col=Population), size=3, alpha=.8)+
  geom_hline(yintercept=0, alpha=.6)+
  geom_vline(xintercept=0, alpha=.6)+
  coord_fixed(x=c(-220,220), y=c(-220,220))+
  scale_color_manual(values=cols)

dev.off()

###Loading plots of influential SNPs

png(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/PCA/loadings_ss_GQ13_90pCall_maf5p_PC1.png", res=300, width = 2400, height = 1600)

loadingplot(pca1$loadings, axis=2)

plot(-log(abs(pca1$loadings[,"Axis1"]))) #thr=quantile(pca1$loadings[,"Axis1"],0.9999, type=8))
#axis determines which DF should be analyzed. thr is used to determine whether labelling should take place, here 99.9% percentile.

#######DISCRIMINANT ANALYSIS OF PRINCIPAL COMPONENTS############

#First identify the most likely number of clusters.

grp <- find.clusters(genlight, max.n.clust = 10)

#Limit the max number of clusters to be checked. At most there can be n-1 clusters. If populations are already known this point can be skipped though it is good practice to start with no assumptions.
#When choosing the no. of PCs to retain opt for a large number >N if computation is not an issue.
#The final choice for cluster number should be at the lowest point of the graph. A "true" value does however never exist, one should try several numbers of clusters.

#Check the independently derived cluster number and assignment against the 'real' populations

table(pop(genlight), grp$grp)

png(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/DAPC/assignment_grid_K6.png", res=150, width = 800, height = 800)

table.value(table(pop(genlight), grp$grp), csize=3, clegend=0, col.lab=paste("Cluster", 1:6), row.lab=c("Miami","Sylt","Varna","Villefranche","Woods Hole"))

dev.off()

#Run the actual DAPC using the clusters identified previously. Here one can also use known populations!

dapc1 <- dapc(genlight, glPca = pca1) #use grp$grp instead of the inherent genlight pops to use other cluster assignments
#When choosing no. of PCs to retain look at the inflection point of the graph. Retaining too many PCs can lead to overfitting, too few to information loss. If unsure take a look at the a-score!
#When choosing the no. of discriminant functions, which should be equal to no.of.clusters-1, retain as many as you realistically would like to display later e.g. 3 dimensions = 3 DF.
#Can also use the output of a previous PCA to save time!!!

#ASIDE: The a-score for determining the amount of PCs to be retained.

#By analysing the a-score the optimal number of PCs to be retined can be determined, similar to cluster number.

#Create a new dapc with the amount of PCs you want to test. No. of DF can be high.
temp <- dapc(genlight, glPca=pca1)

#Run a-score optimization

asc <- optim.a.score(temp, smart=F, n.sim=1000)

#Plot the results

png(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/DAPC/dapc_ss_GQ13_90pCall_maf5p_realpops_DF2+3.png", res=300, width = 1600, height = 1600)

scatter(dapc1, xax=2, yax=3, col=cols, posi.da = "bottomleft", 
        mstree=F, cstar = T, solid=.6, cex=3, clab=0,
        leg=T, posi.leg = "bottomright", txt.leg = c("Miami","Sylt","Varna","Villefranche","Woods Hole")) #txt.leg = paste("Cluster", 1:6))

dev.off()

###ASIDE
#If you are using only 1 dimension i.e. 1 DF then you can display indv. desities along the DF.

png(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/DAPC/DAPC_K2.png", res=150, width = 1200, height = 800)

scatter(dapc1,1,1, col=cols, bg="white", scree.da=FALSE, leg=TRUE, txt.leg = paste("Cluster",1:2), solid=.6)

dev.off()

###Loading plots to determine discriminant SNPs###

png(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/DAPC/loadings_preBQSR_realpops_DF3.png", res=300, width = 2400, height = 1600)

loadingplot(dapc1$var.contr, axis=3, thr=quantile(dapc1$var.contr[,"LD3"],0.9999, type=8))
#axis determines which DF should be analyzed. thr is used to determine whether labelling should take place, here 99.9% percentile.

dev.off()

###STRUCTURE-like composition plots to determine group assignments###

compo1 <- compoplot(dapc1, col=cols)
