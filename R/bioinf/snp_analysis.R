########BASIC SNP STATS##########
library(ggplot2)
library(dplyr)
library(tidyr)
myPng <- function(..., width=6, height=6, res=300) {png(..., width=width*res, height=height*res, res=res)}
pops <- c(rep("Miami",15), rep("Sylt",16), rep("Varna",16), rep("Villefranche",9), rep("Woods Hole",16))
cols <- c("#E41A1C","#984EA3","#FF7F00","#377EB8","#4DAF4A")
threeggcols <- c("#F8766D","#00BA38","#619CFF")
twoggcols <- c("#F8766D","#00BFC4")

###IMPORT METADATA FOR COORDINATE PLOTS#########ON A HUGE FAKE GENOME WITH MULTIPLE SCAFFOLDS THIS DOESNT MAKE SENSE
#since position values in the data is relative to the start of the chromosome we had to calculate global position coordinates which are imported below

glob_coord <- read.table("~/NEC_WORK/reference/glob_coord.txt", col.names = c("CHROM", "length", "glob_coord"), colClasses = c("factor", "integer", "integer"))
gen_size <- tail(glob_coord$length,1)+tail(glob_coord$glob_coord,1)-1

#####CUMSUM OF REFERENCE GENOME########
#calculated the cumsum from the .fai index file of the reference genome. use awk.
glob_coord_full_raw <- read.table("~/NEC_WORK/reference/glob_coord_full.txt", col.names = c("CHROM", "length", "glob_coord"), colClasses = c("factor", "integer", "integer"))
glob_coord_full <- glob_coord_full_raw %>% arrange(desc(length)) %>% mutate(cumsum=cumsum(length))

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/snp_stats/reference_cumsum.png", res=300, width = 6, height = 6)

ggplot(glob_coord_full)+
  geom_line(aes(x=c(1:nrow(glob_coord_full)), y=cumsum))+
  geom_area(data=subset(glob_coord_full, length>=10000), aes(x=c(1:nrow(subset(glob_coord_full, length>=10000))), y=cumsum), fill="springgreen3", alpha=.6)+
  geom_vline(xintercept=nrow(subset(glob_coord_full, length>=10000)), lty="twodash")+
  labs(x="Summed Scaffolds", y="Cumulative Length (Mb)")+
  annotate("text", label=paste(round(max(subset(glob_coord_full, length>=10000)$cumsum)/max(glob_coord_full$cumsum)*100,2),"% of Genome", sep=""), x=750, y=50000000, angle=90)+
  scale_y_continuous(breaks=seq(0,150000000,50000000), labels = seq(0,150,50))+
  scale_x_continuous(breaks=c(0,500,nrow(subset(glob_coord_full, length>=10000)),seq(2000,5000,1000)))

dev.off()

##########GQ distribution############
gq_raw <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/gq_ss.GQ.FORMAT", header = T, colClasses = c("factor","integer", rep("integer",72)), na.strings=".")
#change NAs to zero since a missing genotype is the same as a crap genotype. alo change inv. names to population names
gq_raw[is.na(gq_raw)] <- 0

#calculate mean GQs per snp per pop
gq_means <- gq_raw %>% transmute("Miami"=rowMeans(gq_raw[,3:17]), "Sylt"=rowMeans(gq_raw[,18:33]), "Varna"=rowMeans(gq_raw[,34:49]), "Villefranche"=rowMeans(gq_raw[,50:58]), "Woods Hole"=rowMeans(gq_raw[,59:74]))
#reappend position data
gq_means <- bind_cols(gq_raw[,1:2], gq_means)
#transform data to long format
gq_means <- gather(gq_means, "pop", "mean_gq", 3:7)

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/snp_stats/gq_ss.png", res=300, width = 8, height = 10)

ggplot(gq_means)+
  geom_density(aes(x=mean_gq, fill=pop), bw=0.5, alpha=.8)+
  geom_vline(xintercept = 13, lty=2, color="gray20", alpha=.8)+
  labs(x="Mean Genotype Quality", y="Density")+
  scale_x_continuous(breaks=c(0,10,13,20,30,50,75,100))+
  scale_fill_manual(values=cols, guide=F)+
  facet_grid(pop~.)

dev.off()

#overlaid
#shading density plot...need to do some tricks. Shading an area requires a ymax and ymin value as the area to shade between. Since we calculated density we don't have a yvalue/variable, but we can extract the y-coordinates from the plot itself and use those!!
prelimplot <- ggplot(gq_means)+
                  geom_density(aes(x=mean_gq, color=pop), bw=0.5)+
                  scale_color_manual(values=cols, name="Population")
prelimdata <- ggplot_build(prelimplot)$data[[1]]

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/snp_stats/pop_gq_ss.png", res=300, width = 8, height = 6)

ggplot()+
  geom_density(data=gq_means, aes(x=mean_gq, color=pop), bw=0.5, size=1.2)+
  geom_area(data=subset(prelimdata, colour==cols[1] & x<13), aes(x=x, y=y), color="gray20", alpha=.5)+
  geom_area(data=subset(prelimdata, colour==cols[2] & x<13), aes(x=x, y=y), color="gray20", alpha=.5)+
  geom_area(data=subset(prelimdata, colour==cols[3] & x<13), aes(x=x, y=y), color="gray20", alpha=.5)+
  geom_area(data=subset(prelimdata, colour==cols[4] & x<13), aes(x=x, y=y), color="gray20", alpha=.5)+
  geom_area(data=subset(prelimdata, colour==cols[5] & x<13), aes(x=x, y=y), color="gray20", alpha=.5)+
  labs(x="Mean Genotype Quality", y="Density")+
  scale_x_continuous(breaks=c(0,10,13,20,30,50,75,100))+
  scale_color_manual(values=cols, name="Population")+
  theme(legend.position = c(0.85,0.8))

dev.off()

#could also convert data to long format before taking means (Takes ages...do on cluster?)
gq <- gather(gq_raw,"Pops","GQ",3:74)

##########SNP Density#############
#best to calculate with vcftools as adegenet cannot do a window-based approach. However windows screw up at chromosome/scaffold ends so you should delete those faulty windows, see below

snpden_raw <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/global_snpden_ss_GQ13.snpden", header = T)

snpden <- left_join(snpden_raw, glob_coord, by="CHROM") %>% filter(length-BIN_START>=10000)

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/snp_stats/global_snpden.png", res=300, width = 6, height = 6)

ggplot(snpden)+
  geom_density(aes(x=VARIANTS.KB), bw = 1, fill="coral", alpha=.8)+
  #geom_vline(xintercept = mean(snpden$VARIANTS.KB), lty="twodash", col="gray20")+
  labs(x="SNP Density (per Kb)", y="Proportion of Genome")

dev.off()

######Missing Data#######
##fraction of missing genotypes per snp

miss_sites_raw <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/global_misssite_ss_GQ13.lmiss", header = T)

prelimplot <- ggplot(miss_sites_raw, aes(x=(1-F_MISS)*100))+
          geom_density(adjust=3)

#Here we created the object prelimplot in order to extract information from it. Since we want to shade the area under a curve we need to give ymin and ymax values to the geom_area command.
#However as this is a density plot we do not have y-values!! But we can extract these from the ggplot graphing background data using ggplot_build and looking in the data section

prelimdata <- ggplot_build(prelimplot)$data[[1]]

#Now we can use this information to shade appropriately.

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/snp_stats/global_callrate_ss_GQ13.png", res=300, width = 6, height = 6)

ggplot()+
  geom_density(data=miss_sites_raw, aes(x=(1-F_MISS)*100), adjust=3, alpha=0)+
  labs(x="Call Rate (%)", y="Proportion of SNPs")+
  geom_area(data=subset(prelimdata, x>=70 & x<=98), aes(x=x, y=y), fill=twoggcols[2], alpha=.8)+
  geom_area(data=subset(prelimdata, x>=98), aes(x=x, y=y), fill=twoggcols[1], alpha=.8)
  #annotate("text", label=c("70% Call Rate"), x=80, y=0.1, angle=90)

dev.off()

##fraction of missing genotypes split per population
SY_miss_sites <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/SY_misssites_ss_GQ13.lmiss", header = T )%>% mutate(pop="Sylt")
WH_miss_sites <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/WH_misssites_ss_GQ13.lmiss", header = T )%>% mutate(pop="Woods Hole")
MI_miss_sites <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/MI_misssites_ss_GQ13.lmiss", header = T )%>% mutate(pop="Miami")
VA_miss_sites <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VA_misssites_ss_GQ13.lmiss", header = T )%>% mutate(pop="Varna")
VF_miss_sites <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VF_misssites_ss_GQ13.lmiss", header = T )%>% mutate(pop="Villefranche")

pop_miss_sites <- bind_rows(SY_miss_sites, WH_miss_sites, MI_miss_sites, VA_miss_sites, VF_miss_sites)

prelimplot <- ggplot(pop_miss_sites)+
                geom_density(aes(x=(1-F_MISS)*100, color=pop), adjust=18)+
                scale_color_manual(values=cols)
prelimdata <- ggplot_build(prelimplot)$data[[1]]

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/snp_stats/pop_miss_site_ss_GQ13.png", res=300, width = 6, height = 6)

ggplot()+
  geom_density(data=pop_miss_sites, aes(x=(1-F_MISS)*100, color=pop), adjust=18, size=1.2)+
  geom_area(data=subset(prelimdata, colour==cols[1] & x<70), aes(x=x, y=y), color="gray20", alpha=.5)+
  geom_area(data=subset(prelimdata, colour==cols[2] & x<70), aes(x=x, y=y), color="gray20", alpha=.5)+
  geom_area(data=subset(prelimdata, colour==cols[3] & x<70), aes(x=x, y=y), color="gray20", alpha=.5)+
  geom_area(data=subset(prelimdata, colour==cols[4] & x<70), aes(x=x, y=y), color="gray20", alpha=.5)+
  geom_area(data=subset(prelimdata, colour==cols[5] & x<70), aes(x=x, y=y), color="gray20", alpha=.5)+
  labs(x="Call Rate (%)", y="Density")+
  scale_x_continuous(breaks=c(0,50,60,70,80,90,100))+
  scale_color_manual(values=cols, name="Population")+
  theme(legend.position = c(0.2,0.8))

dev.off()

##fraction of missing data per pop before and after filtering
miss_full <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/global_missindv_ss_GQ13.imiss", header = T) %>% mutate(fltr="Base Set", pop=pops)
SY_miss_70p <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/SY_missindv_ss_GQ13_70pCall.imiss", header = T )%>% mutate(fltr="70% Call Rate", pop="Sylt")
WH_miss_70p <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/WH_missindv_ss_GQ13_70pCall.imiss", header = T )%>% mutate(fltr="70% Call Rate", pop="Woods Hole")
MI_miss_70p <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/MI_missindv_ss_GQ13_70pCall.imiss", header = T )%>% mutate(fltr="70% Call Rate", pop="Miami")
VA_miss_70p <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VA_missindv_ss_GQ13_70pCall.imiss", header = T )%>% mutate(fltr="70% Call Rate", pop="Varna")
VF_miss_70p <- read.table("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VF_missindv_ss_GQ13_70pCall.imiss", header = T )%>% mutate(fltr="70% Call Rate", pop="Villefranche")

#join tables and reorder the fltr column to make sense
miss <- bind_rows(miss_full, SY_miss_70p, WH_miss_70p, MI_miss_70p, VA_miss_70p, VF_miss_70p)
miss$fltr <- factor(miss$fltr, levels = c("Base Set", "70% Call Rate"))


myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/snp_stats/pop_miss_indv.png", res=300, width = 8, height = 6)

ggplot(miss)+
  geom_boxplot(aes(y=F_MISS*100, x=pop, fill=fltr), alpha=.9)+
  labs(x="Population", y="Missing Data per Individual (%)")+
  scale_fill_discrete(name="SNP Set")+
  theme(legend.position = c(0.2,0.8))

dev.off()

#SNP numbers before and after filtering by population
#calculate and add number snps manually. Just copy from vcftools output
Population <- c("Sylt", "Woods Hole", "Miami", "Varna", "Villefranche")
seventyCall <- c(7127219, 6727347, 6257046, 5998825, 3254160)
oneMiss <- c(1356839,748712, 1378873, 968322, 575194)
total <- rep(7310753,5)

filter_counts <- cbind.data.frame(Population, seventyCall, oneMiss, total)

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/snp_stats/pop_snp_filters.png", res=300, width = 6, height = 6)

ggplot(filter_counts)+
  geom_bar(aes(x=Population, y=total), stat="identity", fill="grey20", alpha=.6)+
  geom_bar(aes(x=Population, y=seventyCall), stat="identity", fill=twoggcols[2])+
  geom_bar(aes(x=Population, y=oneMiss), stat= "identity", fill=twoggcols[1])+
  scale_y_continuous(breaks=seq(0,7000000,1000000), labels = c(0:7))+
  labs(x="Population", y="Number of SNPs (Million)")
  
dev.off()

#can use this function to create darker versions of our colors
#darken <- function(color, factor){
#  col <- col2rgb(color)
#  col <- col/factor
#  col <- rgb(t(col), maxColorValue=255)
#  col
#}