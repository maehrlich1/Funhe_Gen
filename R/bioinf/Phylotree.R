###Call libs###
library(poppr)
library(adegenet)
library(vcfR)
library(ape)
library(parallel)
myPng <- function(..., width=6, height=6, res=300) {png(..., width=width*res, height=height*res, res=res)}
pops <- c(rep("Miami",15), rep("Sylt",16), rep("Varna",16), rep("Villefranche",9), rep("Woods Hole",16))
popcols <- c(rep("#E41A1C",15), rep("#984EA3",16), rep("#FF7F00",16), rep("#377EB8",9), rep("#4DAF4A",16))
sample_size <- cbind.data.frame(Population=c("Sylt","Woods Hole","Miami","Varna","Villefranche"), sample_size=c(16,16,15,16,9))
cols <- c("#E41A1C","#984EA3","#FF7F00","#377EB8","#4DAF4A")
threeggcols <- c("#F8766D","#00BA38","#619CFF")
twoggcols <- c("#F8766D","#00BFC4")
set.seed(1)

####PREPARATION OF GENLIGHT OBJECT#########
#take SNP set with minimal missing data
###Read VCF and convert to genlight
genlight <- vcfR2genlight(read.vcfR("/Users/Moritz/NEC_WORK/variants/merged/no_mito/master_snps_ss_GQ13_1miss_maf5p_bial.vcf"))
##Append pop data
pop(genlight) <- pops

##########PLOTTING A PHYLOGENETIC TREE########
#the poppr function bitwise.dist() can operate directly on a genlight object to calculate a distance matrix
D <- bitwise.dist(genlight)

#OPTIONAL: Calculate distance matrix by hand in R
#we don't really use the entire genlight object but rather just a genotype matrix so extracting this information is easy
#geno <- as.matrix(genlight)

#then create the distance matrix. its just a normal R function
#D <- dist(geno)

#OPTIONAL 2: split huge genlight object into smaller blocks to make calculations more manageable. Each block can have approx. 10,000 SNPs.
#blocks <- seploc(genlight, n.block=117)
#names(blocks) #just to check

#compute pairwise dist for each and sum up across all blocks
#dist_list <- lapply(blocks, function(e) dist(as.matrix(e)))
#D <- Reduce("+", dist_list)

#CLUSTERING ALGORITHM
#Neighbour-joining algorithm. bionjs() is an improved NJ algorithm which
#ladderizing just orders the branches
tree <- bionjs(D) %>% ladderize()

###PLOTTING. Here we erase the tiplabels and use the custom color instead

myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/dendro/NJ_unrooted_ss_GQ13_1miss_maf5p_bial.png", res=300, width = 8, height = 8)

plot(tree, typ="unrooted", show.tip=F)
tiplabels(pch=19, col=popcols)
legend(x=-.05, y=.18, c("Miami","Sylt","Varna","Villefranche","Woods Hole"), pch=19, col=cols)

dev.off()

#Now for a rooted tree, with default root

myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/dendro/NJ_rooted_ss_GQ13_90pCall_maf5p.png", res=300, width = 6, height = 6)

plot(tree, show.tip=F)
tiplabels(pch=19, col=popcols)
legend("bottomleft", c("Miami","Sylt","Varna","Villefranche","Woods Hole"), pch=19, pt.cex=1.2, col=cols)
add.scale.bar(0.15,70)

dev.off()

#can also specify the root by clicking on the plot!

roottree <- root(tree,interactive = T)

myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/dendro/NJ_rooted_ss_GQ13_1miss_maf5p_bial.png", res=300, width = 6, height = 6)

plot(roottree, show.tip=F)
tiplabels(pch=20, col=popcols, cex=1.2)
legend(0.15, 70, c("Miami","Sylt","Varna","Villefranche","Woods Hole"), pch=19, col=cols)
add.scale.bar(0,20)

dev.off()

#Hierarchical clustering algorithm, the NJ method is the preferred distance-based clustering algorithm but this can be used in some cases too
#hclust is not a genetics-specific function so we have to convert to a phylo object afterwards (nj() is part of the ape package and does this automatically)
P <- as.phylo(hclust(D, method = "ward.D2"))

myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/dendro/HC_unrooted_ss_GQ13_1miss_maf5p_bial.png", res=300, width = 6, height = 8)

plot(P, type="unrooted", label.offset=50, show.tip=F)
tiplabels(pch=19, col=popcols)
legend(.4,.8, c("Miami","Sylt","Varna","Villefranche","Woods Hole"), pch=19, col=cols)

dev.off()

############BOOTSTRAPPING###############
#to check branch support we can subsample SNPs and run the analysis again many hundreds or thousands of times.
#if the branches still resolve equally this gives us high confidence in our tree

#this function outputs bootstrap statistics per node.
boottree <- aboot(genlight, sample=1000, distance="bitwise.dist", tree="bionjs", root=F) %>% ladderize()

#extract the bootstrap values and only use those above 95 in this case
boots <- ifelse(round(boottree$node.label)>=95, round(boottree$node.label), NA)
#create a  vector where the symbol type is given for bootstrap values above 95 (here 25=downward triangle)
boots_symbol <- ifelse(round(boottree$node.label)>=95, 25, NA)

#plot
myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/dendro/NJ_unrooted_boot95_ss_GQ13_1miss_maf5p_bial.png", res=300, width = 6, height = 8)

plot(boottree, typ="unrooted", show.tip=F)
tiplabels(pch=19, col=popcols)
nodelabels(pch=boots_symbol, cex=1, bg="yellow", frame="none")
legend(x=-.02, y=.18, c("Miami","Sylt","Varna","Villefranche","Woods Hole"), pch=19, col=cols)

dev.off()

#OPTIONAL: Use the ape package bootstrapping algorithm instead.As input it takes the tree, the genotype matrix used and the function that was used on the genotype matrix to get the tree.
#Here we do 10 bootstrap runs
#boots <- boot.phylo(tree, geno, FUN= function(x) bionjs(dist(x)), B=10)

#makes sense to run the above command on the cluster. can either write the output to a table with write.table() or save the actual boots object using save()
#load the cluster output here

#load("/Users/Moritz/NEC_WORK/r_output/boots_1000_NJ_ss_GQ13_90pCall_maf5p.boots")
#boots <- round(boots/10,0) #to get a percentage bootstrap support

#now replot the tree and add bootstrap values to the nodes

#myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/dendro/NJ_unrooted_boot_ss_GQ13_90pCall_maf5p.png", res=300, width = 6, height = 6)

#plot(tree, typ="unrooted", show.tip=F)
#tiplabels(pch=19, col=popcols)
#legend("bottomright", c("Miami","Sylt","Varna","Villefranche","Woods Hole"), pch=19, cex=0.7, col=cols)
#nodelabels(boots, cex=1, bg="transparent", frame="none")

#dev.off()
