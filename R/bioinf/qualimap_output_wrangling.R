#######Qualimap output wrangling#########
library(dplyr)
library(ggplot2)
myPng <- function(..., width=6, height=6, res=300) {png(..., width=width*res, height=height*res, res=res)}
pops <- c(rep("Miami",15), rep("Sylt",16), rep("Varna",16), rep("Villefranche",9), rep("Woods Hole",16))
sample_size <- cbind.data.frame(Population=c("Sylt","Woods Hole","Miami","Varna","Villefranche"), sample_size=c(16,16,15,16,9))
cols <- c("#E41A1C","#984EA3","#FF7F00","#377EB8","#4DAF4A")
threeggcols <- c("#F8766D","#00BA38","#619CFF")
twoggcols <- c("#F8766D","#00BFC4")
#Iport metadata for coordinate plots. CAREFUL: ON A HUGE FAKE GENOME WITH MULTIPLE SCAFFOLDS THIS DOESNT MAKE SENSE
#since position values in the data are relative to the start of the chromosome we had to calculate global position coordinates which are imported below
glob_coord <- read.table("~/NEC_WORK/reference/glob_coord.txt", col.names = c("CHROM", "length", "glob_coord"), colClasses = c("factor", "integer", "integer"))
gen_size <- tail(glob_coord$length,1)+tail(glob_coord$glob_coord,1)-1

######Mapping quality across reference######

mq_across_ref <- read.delim("/Users/Moritz/NEC_WORK/mapped/merged/qualimap/mq_across_ref.txt")

myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/mapping_stats/mapping_quality_across_reference.png", res=300, width = 8, height = 6)

ggplot(mq_across_ref, aes(x=Pos,y=MQ))+
  geom_point(col="gray20", alpha=.005)+
  geom_smooth(aes(group=Population, col=Population, fill=Population))+
  labs(x="Genomic Position (Mb)", y="Mapping Quality")+
  scale_x_continuous(breaks = seq(0,140000000,20000000), labels = seq(0,140,20))+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)+
  theme(legend.position = c(0.75,0.2))

dev.off()

#####Coverage across reference#####

cov_across_ref <- read.delim("/Users/Moritz/NEC_WORK/mapped/merged/qualimap/cov_across_ref.txt")

myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/mapping_stats/coverage_across_reference.png", res=300, width = 8, height = 6)

ggplot(cov_across_ref, aes(x=Pos,y=Cov))+
  geom_point(col="gray20", alpha=.005)+
  geom_smooth(aes(group=Population, col=Population, fill=Population))+
  coord_cartesian(y=c(0,140))+
  labs(x="Genomic Position (Mb)", y="Sequencing Depth")+
  scale_x_continuous(breaks = seq(0,140000000,20000000), labels = seq(0,140,20))+
  scale_y_continuous(breaks = c(0,10,20,30,40,50,75,100,125))+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)+
  theme(legend.position = c(0.75,0.8))

dev.off()

####Mapping Quality Histogram####

mq_hist_raw <- read.delim("/Users/Moritz/NEC_WORK/mapped/merged/qualimap/mq_hist.txt")
 
mq_hist <- mq_hist_raw %>% group_by(MQ,Population) %>% summarise(sum_reads=sum(Num_Reads))

myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/mapping_stats/mapping_quality_histogram.png", res=300, width = 6, height = 12)

ggplot(mq_hist)+
  geom_histogram(aes(x=MQ, y=sum_reads, fill=Population), stat = "identity", binwidth = 1)+
  labs(x="Mapping Quality", y="Number of Reads")+
  scale_fill_manual(values=cols, guide=F)+
  facet_grid(Population~.)

dev.off()

########Coverage Histogram########

cov_hist_raw <- read.delim("/Users/Moritz/NEC_WORK/mapped/merged/qualimap/cov_hist.txt")

cov_hist <- cov_hist_raw %>% group_by(Cov,Population) %>% summarise(sum_loci=sum(Num_Reads)) %>%
                              full_join(sample_size, by="Population") %>%
                              mutate(prop_loci = sum_loci/(gen_size*sample_size)*100)

myPng(filename = "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/mapping_stats/coverage_histogram.png", res=300, width = 6, height = 12)

ggplot(cov_hist)+
  geom_histogram(aes(x=Cov, y=prop_loci, fill=Population), stat="identity", binwidth = 1)+
  coord_cartesian(x=c(0,60))+
  labs(x="Sequencing Depth", y="Proportion of Genome (%)")+
  scale_x_continuous(breaks = seq(0,60,10))+
  scale_fill_manual(values=cols, guide=F)+
  facet_grid(Population~.)

dev.off()

####Population means of coverage
pop_covmean <- cov_hist %>% group_by(Population) %>% summarise(meancov = sum(Cov*sum_loci/sum(sum_loci)))

#Mean percentage of non-covered loci
sum(subset(cov_hist, Cov==0)$sample_size*subset(cov_hist, Cov==0)$prop_loci)/72
