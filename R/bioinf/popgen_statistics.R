###########POPGEN STATISTICS#############
library(pegas)
library(ggplot2)
library(dplyr)
library(tidyr)
myPng <- function(..., width=6, height=6, res=300) {png(..., width=width*res, height=height*res, res=res)}
pops <- c(rep("Miami",15), rep("Sylt",16), rep("Varna",16), rep("Villefranche",9), rep("Woods Hole",16))
sample_size <- cbind.data.frame(Population=c("Sylt","Woods Hole","Miami","Varna","Villefranche"), sample_size=c(16,16,15,16,9))
cols <- c("#E41A1C","#984EA3","#FF7F00","#377EB8","#4DAF4A")
threeggcols <- c("#F8766D","#00BA38","#619CFF")
twoggcols <- c("#F8766D","#00BFC4")
###IMPORT METADATA FOR COORDINATE PLOTS#########ON A HUGE FAKE GENOME WITH MULTIPLE SCAFFOLDS THIS DOESNT MAKE SENSE
#since position values in the data is relative to the start of the chromosome we had to calculate global position coordinates which are imported below
glob_coord <- read.table("~/NEC_WORK/reference/glob_coord.txt", col.names = c("CHROM", "length", "glob_coord"), colClasses = c("factor", "integer", "integer"))

#HERE ACTUAL STATISTICS ARE CALCULATED RATHER THAN PARAMETERS CALCULATED#

####PREPARATION OF GENLIGHT OBJECT#########
###Read VCF and convert to genlight
genlight <- vcfR2genlight(read.vcfR("/Users/Moritz/NEC_WORK/variants/merged/no_mito/master_snps_ss_GQ13_1miss_maf5p_bial.vcf"))
##Append pop data
pop(genlight) <- pops

#create distance matrix
D <- bitwise.dist(genlight)

#create a minimum spanning network
msn <- poppr.msn(genlight, D)

#plot it with more options
plot_poppr_msn(genlight, msn, inds="somethingrandom", palette = cols)

###play

D <- dist(as.matrix(genlight))
