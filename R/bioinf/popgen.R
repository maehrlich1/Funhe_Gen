###########POPGEN ANALYSIS#############
library(ggplot2)
library(dplyr)
library(tidyr)
library(qvalue)
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

##########AF distribution###########
####Plot the number of SNPs that are bi-,tri- or n-allelic###
freq_raw <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/global_freq_ss_GQ13_70pCall.frq", col.names = c("CHROM", "POS", "N_ALLELES", "N_CHR", paste("freq_",1:5, sep="")), header=F, skip=1)

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/allele_num_ss_GQ13_70pCall.png", res=300, width = 6, height = 6)

ggplot(freq_raw)+
  geom_bar(aes(x=N_ALLELES), fill="coral", colour="gray20", alpha=.8)+
  labs(x="Number of Alleles", y="Number of SNPs (millions)")+
  scale_y_continuous(breaks = seq(0,6000000,1000000), labels = c(0:6))

dev.off()

###Global AF distribution (Multi-allelic)###
freq <- freq_raw %>% gather("allele","freq", 5:9, na.rm=T)

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/global_af_ss_GQ13_70pCall.png", res=300, width = 6, height = 6)

ggplot(freq)+
  geom_histogram(aes(x=freq), binwidth=0.01, fill="coral", color="gray20", alpha=.8)+
  labs(x="Allele Frequency", y="Number of SNPs (millions)")+
  scale_y_continuous(breaks = seq(0,1600000,200000), labels = seq(0,1.6,0.2))

dev.off()

#Folded global MAF-distribution (Biallelic)
#don't include MAF=0.5 as here the minor allele does not exist

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/global_maf_ss_GQ13_70pCall.png", res=300, width = 6, height = 6)

ggplot(subset(freq, freq<0.5))+
  geom_histogram(aes(x=freq), binwidth=0.01, fill="coral", color="gray20", alpha=.8)+
  labs(x="Minor Allele Frequency", y="Number of SNPs (millions)")+
  scale_y_continuous(breaks = seq(0,1600000,200000), labels = seq(0,1.6,0.2))

dev.off()

#POP-WISE MAF SITE FREQUENCY SPECTRUM (Folded SFS) (biallelic)

#CAREFUL!! The MAF is not a continuous variable if we have a finite sample size! 
#E.g. the lowest possible MAF is 1/2n. As a result our binwidth should be 1/2n, however we have not always sampled 2n chromosomes beacause of missing data! 
#The discretisation is therefore specific to each SNP and its sampling depth!!!
#best to work with populations of equal sample size (downsample if necessary) and with 100% call rate! See reduced_n.R script...

########HETEROZYGOSITY AND INBREEDING COEFFICIENT##########
#We know that LD between closely spaced SNPs can distort reality.
#SNPs in LD will be weighted as much as those that are not although they do not contain additional information.
#This can be avoided by setting up a minimum distance between SNPs used in the analysis. Ideally this distance should be based on previously measured LD but we can set an arbitrary (large) distance too.
#Given that we have to reduce the SNP set by thinning we might aswell pick those SNPs that are A: polymorphic and B: have a high call rate

##Individual Multilocus Heterozygosity (MLH = Proportion of heterozygous markers). This evaluates the total degree of hetero/homozygosity and confounds the effect of identity by descent with "normal" HW processes
##However, MLH is dependent on sample size. In order to get an unbiased result we need to downsample to a common sample size.
##See red_n script

##Individual Multilocus Inbreeding Coefficient, F (MoM Calculation Yang, 2010). This estimates the degree of homozygosity that is due to two alleles being identical by descent
#In summary this estimator compares the NUMBER OF SITES per Individual that are homozygous to the number of sites that are expected given no inbreeding. A per individual statistic!
SY_het_thin <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/SY_het_thin_ss_GQ13_1miss_bial.het") %>% mutate(Population="Sylt")
WH_het_thin <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/WH_het_thin_ss_GQ13_1miss_bial.het") %>% mutate(Population="Woods Hole")
MI_het_thin <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/MI_het_thin_ss_GQ13_1miss_bial.het") %>% mutate(Population="Miami")
VA_het_thin <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VA_het_thin_ss_GQ13_1miss_bial.het") %>% mutate(Population="Varna")
VF_het_thin <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VF_het_thin_ss_GQ13_1miss_bial.het") %>% mutate(Population="Villefranche")

pop_indv_F <- bind_rows(SY_het_thin, WH_het_thin, MI_het_thin, VA_het_thin, VF_het_thin)

#plot for fun but we can add the population based F too
ggplot(pop_indv_F)+
  geom_hline(yintercept=0, color="gray20", alpha=.8, lty="twodash")+
  geom_boxplot(aes(x=Population, y=F, fill=Population))+
  scale_fill_manual(values=cols, guide=F)+
  labs(x="Population", y=expression(Inbreeding~Coefficient~"("*italic(F)*")"))

##Population wide Inbreeding Coefficient after Weir and Cockerham (1984). A per population statistic!
#Here we calculate the inbreeding coefficient for every locus by comparing the observed number of heterozygotes to the expected number of hets in that population.
#Also we should correct the expected number of hets for sample size. (Expect more hets at smaller sample sizes) this is incorporated in the WC84 estimator

SY_hardy_thin <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/SY_hardy_thin_ss_GQ13_1miss_bial.hwe") %>% mutate(Population="Sylt")
WH_hardy_thin <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/WH_hardy_thin_ss_GQ13_1miss_bial.hwe") %>% mutate(Population="Woods Hole")
MI_hardy_thin <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/MI_hardy_thin_ss_GQ13_1miss_bial.hwe") %>% mutate(Population="Miami")
VA_hardy_thin <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VA_hardy_thin_ss_GQ13_1miss_bial.hwe") %>% mutate(Population="Varna")
VF_hardy_thin <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VF_hardy_thin_ss_GQ13_1miss_bial.hwe") %>% mutate(Population="Villefranche")

pop_Fis <- bind_rows(SY_hardy_thin, WH_hardy_thin, MI_hardy_thin, VA_hardy_thin, VF_hardy_thin) %>%
                  separate(OBS.HOM1.HET.HOM2., c("obs_hom1", "obs_het", "obs_hom2"), sep="/", convert=T) %>%
                  separate(E.HOM1.HET.HOM2., c("exp_hom1", "exp_het", "exp_hom2"), sep="/", convert=T) %>%
                  inner_join(sample_size, by="Population") %>%
                  group_by(Population) %>%
                  summarise(mean_obs_het = mean(obs_het), mean_denom = mean((((sample_size/(sample_size-1))*exp_het)-((1/(2*sample_size-2))*obs_het)))) %>%
                  mutate(pop_Fis = 1-(mean_obs_het/mean_denom))

#as you can see Fis is estimated for multilocus data by taking the ratio of means. this performs better than calculating Fis for every SNP and then taking the average. The same goes for Fst, hence weighted Fst!!
#add our calculated Fis to the data table from before to include it in our plot

combo_Fis <- inner_join(pop_indv_F, pop_Fis, by="Population")

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/pop_F_thin_ss_GQ13_1miss_bial.png", res=300, width = 6, height = 6)

ggplot(combo_Fis)+
  geom_hline(yintercept=0, color="gray20", alpha=.8, lty="twodash")+
  geom_boxplot(aes(x=Population, y=F, fill=Population))+
  geom_point(aes(x=Population, y=pop_Fis), fill="yellow", pch=23, size=2)+
  scale_fill_manual(values=cols, guide=F)+
  labs(x="Population", y=expression(Inbreeding~Coefficient~"("*italic(F)*")"))

dev.off()

########SITE-SPECIFIC HWE TESTING#########
#Here we use the HWE exact test (Wigginton, Cutler and Abecasis (2005)).
#GLOBAL TEST (likely to fail as we have strong pop structure i.e. lack of hets)
hwe_raw <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/global_hardy_ss_GQ13_70pCall.hwe")

#Next use the package qvalue to convert pvalues to qvalues
hwe <- hwe_raw %>% mutate(qval_hwe = qvalue(hwe_raw$P_HWE)$qvalues) %>%
                    mutate(qval_het_def = qvalue(hwe_raw$P_HET_DEFICIT)$qvalues) %>%
                    mutate(qval_het_ex = qvalue(hwe_raw$P_HET_EXCESS)$qvalues)

prelimplot <- ggplot(hwe)+
                  geom_density(aes(x=-log(qval_hwe)), bw=0.5)

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/global_hardy_ss_GQ13_70pCall.png", res=300, width = 6, height = 6)

ggplot()+
  geom_density(data=hwe, aes(x=-log(qval_hwe)), bw=0.5)+
  geom_vline(xintercept = -log(0.05), lty="twodash")+
  geom_area(data=subset(ggplot_build(prelimplot)$data[[1]], x<(-log(0.05))), aes(x=x, y=y), fill="gray20", alpha=.8)+
  geom_area(data=subset(ggplot_build(prelimplot)$data[[1]], x>(-log(0.05))), aes(x=x, y=y), fill="coral", alpha=.8)+
  labs(x=expression(-log[10](italic(q)*"-"*"value")), y="Density")

dev.off()

#Compare Het Excess vs Het Deficit - in which direction are the points that fail HWE biased?
hwe_fail <- hwe %>% filter(qval_hwe < 0.05) %>%
                    mutate(bias_dir=ifelse(qval_het_def<qval_het_ex,"Het. Deficit","Het. Excess")) %>%
                    mutate(bias_mag=ifelse(bias_dir=="Het. Deficit", qval_het_def, qval_het_ex))

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/global_hardy_biasdir_ss_GQ13_70pCall.png", res=300, width = 8, height = 6)

ggplot(hwe_fail)+
  geom_histogram(aes(x=-log(bias_mag), fill=bias_dir), binwidth=1)+
  labs(x=expression(-log[10](italic(q)*"-"*"value")), y="Count (100,000)")+
  scale_y_continuous(breaks=seq(0,250000,50000), labels=seq(0,2.5,0.5))+
  scale_fill_discrete(name="Bias Direction")+
  theme(legend.position = c(0.8,0.8))

dev.off()

##Population specific HW testing
#SY_hwe_raw <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/SY_hardy_ss_GQ13_70pCall.hwe")
SY_hwe <- SY_hwe_raw %>% filter(!is.na(ChiSq_HWE)) %>% #need to make sure qvalues are not based on artifical pvalues. if these nan are not removed we get a lot of high pvalues that are not true
                          mutate(qval_hwe = qvalue(P_HWE)$qvalues) %>%
                          mutate(qval_het_def = qvalue(P_HET_DEFICIT)$qvalues) %>%
                          mutate(qval_het_ex = qvalue(P_HET_EXCESS)$qvalues) %>%
                          mutate(Population = "Sylt")

#WH_hwe_raw <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/WH_hardy_ss_GQ13_70pCall.hwe")
WH_hwe <- WH_hwe_raw %>% filter(!is.na(ChiSq_HWE)) %>%
                          mutate(qval_hwe = qvalue(P_HWE)$qvalues) %>%
                          mutate(qval_het_def = qvalue(P_HET_DEFICIT)$qvalues) %>%
                          mutate(qval_het_ex = qvalue(P_HET_EXCESS)$qvalues) %>%
                          mutate(Population = "Woods Hole")

#MI_hwe_raw <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/MI_hardy_ss_GQ13_70pCall.hwe")
MI_hwe <- MI_hwe_raw %>% filter(!is.na(ChiSq_HWE)) %>%
                          mutate(qval_hwe = qvalue(P_HWE)$qvalues) %>%
                          mutate(qval_het_def = qvalue(P_HET_DEFICIT)$qvalues) %>%
                          mutate(qval_het_ex = qvalue(P_HET_EXCESS)$qvalues) %>%
                          mutate(Population = "Miami")

#VA_hwe_raw <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VA_hardy_ss_GQ13_70pCall.hwe")
VA_hwe <- VA_hwe_raw %>% filter(!is.na(ChiSq_HWE)) %>%
                          mutate(qval_hwe = qvalue(P_HWE)$qvalues) %>%
                          mutate(qval_het_def = qvalue(P_HET_DEFICIT)$qvalues) %>%
                          mutate(qval_het_ex = qvalue(P_HET_EXCESS)$qvalues) %>%
                          mutate(Population = "Varna")

#VF_hwe_raw <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VF_hardy_ss_GQ13_70pCall.hwe")
VF_hwe <- VF_hwe_raw %>% filter(!is.na(ChiSq_HWE)) %>%
                          mutate(qval_hwe = qvalue(P_HWE)$qvalues) %>%
                          mutate(qval_het_def = qvalue(P_HET_DEFICIT)$qvalues) %>%
                          mutate(qval_het_ex = qvalue(P_HET_EXCESS)$qvalues) %>%
                          mutate(Population = "Villefranche")

pop_hwe <- bind_rows(SY_hwe, WH_hwe, MI_hwe, VA_hwe, VF_hwe)

#All loci

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/pop_hardy_ss_GQ13_70pCall.png", res=300, width = 6, height = 10)

ggplot(pop_hwe)+
  geom_density(aes(x=-log(qval_hwe), fill=Population), adjust=8, alpha=.8)+
  labs(x=expression(-log[10](italic(q)*"-"*"value")), y="Density")+
  geom_vline(xintercept = -log(0.05), lty="twodash", col="gray20")+
  facet_grid(Population~.)+
  scale_fill_manual(values=cols, guide=F)

dev.off()

#Loci Failing HWE
pop_hwe_fail <- pop_hwe %>% filter(qval_hwe < 0.05) %>%
                            mutate(bias_dir=ifelse(qval_het_def<qval_het_ex,"Het. Deficit","Het. Excess")) %>%
                            mutate(bias_mag=ifelse(bias_dir=="Het. Deficit", qval_het_def, qval_het_ex))

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/pop_hardy_biasdir_ss_GQ13_70pCall.png", res=300, width = 8, height = 10)

ggplot(pop_hwe_fail)+
  geom_histogram(aes(x=-log(bias_mag), fill=bias_dir), breaks=seq(-log(0.05),-log(0.05)+6.5, by=.5))+
  labs(x=expression(-log[10](italic(q)*"-"*"value")), y="Count (10,000)")+
  scale_y_continuous(breaks=seq(0,60000,10000), labels=seq(0,6,1))+
  scale_fill_discrete(name="Bias Direction")+
  coord_cartesian(x=c(-log(0.05),-log(0.05)+6.5))+
  facet_grid(Population~.)

dev.off()

##########Watterson's Theta##############
#This is the number of segregating sites per population scaled by the sample size of individuals.
#We can get the number of seg sites straight from the vcf file or (if we want windowed data) from the snp density calculations
SY_dens_nomono <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/SY_dens_nomono_ss_GQ13_70pCall.snpden")%>% mutate(Population= "Sylt")
WH_dens_nomono <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/WH_dens_nomono_ss_GQ13_70pCall.snpden")%>% mutate(Population= "Woods Hole")
MI_dens_nomono <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/MI_dens_nomono_ss_GQ13_70pCall.snpden")%>% mutate(Population= "Miami")
VA_dens_nomono <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VA_dens_nomono_ss_GQ13_70pCall.snpden")%>% mutate(Population= "Varna")
VF_dens_nomono <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VF_dens_nomono_ss_GQ13_70pCall.snpden")%>% mutate(Population= "Villefranche")

#join everything and remove non-sensical bins
pop_thetaw <- bind_rows(SY_dens_nomono,WH_dens_nomono,MI_dens_nomono,VA_dens_nomono,VF_dens_nomono) %>%
  left_join(glob_coord, by="CHROM") %>%
  filter(length-BIN_START>=10000) %>%
  full_join(sample_size, by="Population")

#We have to write our own function for Watterson's Theta. The problem is that R interprets seq(VECTOR) as seq(1:length(VECTOR)) rather than using every element of the vector as an input to create a seq. Have to use sapply to return a vector.
theta.w <- function(sites,sample_size) {
  N <- 2*sample_size
  correctors <- sapply(N, function(x) sum(1/seq(x-1)))
  sites/correctors
}

#now calculate theta_w using our own function. we can even scale the result to the length of the sequence under consideration.
pop_thetaw <- pop_thetaw %>% mutate(theta_w.KB=theta.w(VARIANTS.KB, sample_size)) %>% mutate(theta_w=theta_w.KB/1000)

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/pop_thetaW_ss_GQ13_70pCall.png", res=300, width = 6, height = 10)

ggplot(pop_thetaw)+
  geom_density(aes(x=theta_w, fill=Population), alpha=.8)+
  labs(y="Density", x=expression("Scaled Watterson's Theta,"~italic(theta[W])))+
  facet_grid(Population~.)+
  scale_fill_manual(values=cols, guide=F)+
  coord_cartesian(x=c(0,0.015))

dev.off()

################NUCLEOTIDE DIVERSITY################
#global nucleotide diversity
pi_raw <- read.delim("~/NEC_WORK/variants/merged/no_mito/global_nucdiv_ss_GQ13_70pCall.windowed.pi")

#the vcftools command calculates pi past the ends of the scaffolds. need to delete those bins that dont lie within the true scaffold
pi <- pi_raw %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0)

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/global_nucdiv_ss_GQ13_70pCall.png", res=300, width = 6, height = 6)

ggplot(pi)+
  geom_density(aes(x=PI, y=..density..), fill="coral", alpha=.8)+
  labs(x=expression("Nucleotide Diversity, "*pi), y="Density")+
  geom_vline(xintercept=median(pi$PI), lty=2, color="gray20")+
  stat_function(fun=dnorm, args=list(mean=mean(pi$PI), sd=sd(pi$PI)), col="blue", lty=2)+
  annotate("text", label=paste("Median =", signif(median(pi$PI), 3)), x=1.1*median(pi$PI), y=50, angle=90, color="gray20")

dev.off()

#population-specific nucleotide diversity, also remove nonsensical bins
SY_pi_raw <- read.table("~/NEC_WORK/variants/merged/no_mito/SY_nucdiv_ss_GQ13_70pCall.windowed.pi", header=T)
SY_pi <- SY_pi_raw %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(pop="Sylt")

WH_pi_raw <- read.table("~/NEC_WORK/variants/merged/no_mito/WH_nucdiv_ss_GQ13_70pCall.windowed.pi", header=T)
WH_pi <- WH_pi_raw %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(pop="Woods Hole")

MI_pi_raw <- read.table("~/NEC_WORK/variants/merged/no_mito/MI_nucdiv_ss_GQ13_70pCall.windowed.pi", header=T)
MI_pi <- MI_pi_raw %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(pop="Miami")

VA_pi_raw <- read.table("~/NEC_WORK/variants/merged/no_mito/VA_nucdiv_ss_GQ13_70pCall.windowed.pi", header=T)
VA_pi <- VA_pi_raw %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(pop="Varna")

VF_pi_raw <- read.table("~/NEC_WORK/variants/merged/no_mito/VF_nucdiv_ss_GQ13_70pCall.windowed.pi", header=T)
VF_pi <- VF_pi_raw %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(pop="Villefranche")

#now join all of the individual pop data tables
pop_pi <- bind_rows(SY_pi,WH_pi,MI_pi,VA_pi,VF_pi)

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/pop_nucdiv_ss_GQ13_70pCall.png", res=300, width = 6, height = 10)

ggplot(pop_pi)+
  geom_density(aes(x=PI, fill=pop), alpha=.8)+
  facet_grid(pop~.)+
  scale_fill_manual(values=cols, guide=F)+
  coord_cartesian(x=c(0,0.015))+
  labs(x=expression("Nucleotide Diversity, "*pi), y="Density")

dev.off()

#site-specific nucleotide diversity...prefer windowed

#all_sitepi <- read.table("~/NEC_WORK/variants/merged/all_preBQSR.sites.pi", header=T)
#all_sitepi <- mutate(all_sitepi, ID=seq(1,nrow(all_sitepi),1))

#all_sitepi_plot <- ggplot(all_sitepi)+
#geom_point(aes(x=ID, y=PI), alpha=.2)+
#geom_hline(yintercept=mean(all_sitepi$PI), col="red")+
#coord_cartesian(ylim=c(0,1))

###########TAJIMAS D#####################
#Tajimas D is the ratio of Wattersons theta and the nucleotide diversity. It is quite robust against sample size and can be an indicator of selection/demography.
SY_tajD_raw <- read.table("~/NEC_WORK/variants/merged/no_mito/SY_tajD_ss_GQ13_70pCall.Tajima.D", header=T)
SY_tajD <- SY_tajD_raw %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-(BIN_START+10000)>=0) %>% mutate(pop="Sylt", pos=BIN_START+glob_coord)

WH_tajD_raw <- read.table("~/NEC_WORK/variants/merged/no_mito/WH_tajD_ss_GQ13_70pCall.Tajima.D", header=T)
WH_tajD <- WH_tajD_raw %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-(BIN_START+10000)>=0) %>% mutate(pop="Woods Hole", pos=BIN_START+glob_coord)

MI_tajD_raw <- read.table("~/NEC_WORK/variants/merged/no_mito/MI_tajD_ss_GQ13_70pCall.Tajima.D", header=T)
MI_tajD <- MI_tajD_raw %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-(BIN_START+10000)>=0) %>% mutate(pop="Miami", pos=BIN_START+glob_coord)

VA_tajD_raw <- read.table("~/NEC_WORK/variants/merged/no_mito/VA_tajD_ss_GQ13_70pCall.Tajima.D", header=T)
VA_tajD <- VA_tajD_raw %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-(BIN_START+10000)>=0) %>% mutate(pop="Varna", pos=BIN_START+glob_coord)

VF_tajD_raw <- read.table("~/NEC_WORK/variants/merged/no_mito/VF_tajD_ss_GQ13_70pCall.Tajima.D", header=T)
VF_tajD <- VF_tajD_raw %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-(BIN_START+10000)>=0) %>% mutate(pop="Villefranche", pos=BIN_START+glob_coord)

#now join all of the individual pop data tables
pop_tajD <- bind_rows(SY_tajD,WH_tajD,MI_tajD,VA_tajD,VF_tajD)

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/pop_tajD_ss_GQ13_70pCall.png", res=300, width = 6, height = 10)

ggplot(pop_tajD)+
  geom_density(aes(x=TajimaD, fill=pop), alpha=.8)+
  geom_vline(xintercept = 0, color="gray20", alpha=.6)+
  labs(x=expression("Tajima's"~italic(D)~"(10 Kb Windows)"), y="Density")+
  scale_fill_manual(values=cols, guide=F)+
  coord_cartesian(x=c(-4,4))+
  facet_grid(pop~.)

dev.off()

#can also plot Tajimas D along the genome. SHOULDNT DO THIS SINCE ORDER OF CONTIGS NOT KNOWN!!!

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/pop_tajD_manha_nomono_ss_GQ13_70pCall.png", res=300, width = 12, height = 10)

ggplot(pop_tajD)+
  geom_point(aes(x=pos, y=TajimaD, color=pop), alpha=.4)+
  geom_hline(yintercept = 0, color="gray20", alpha=.8)+
  labs(y=expression("Tajima's"~italic(D)~"(10 Kb Windows)"), x="Arbitrary Position (Mb)")+
  scale_x_continuous(breaks=seq(0,150000000,50000000), labels = seq(0,150,50))+
  scale_color_manual(values=cols, guide=F)+
  facet_grid(pop~.)

dev.off()


#############F STATISTICS######################
#Here the output of vcftools fst calculations (Weir and Cockerham 1984) are processed.
#in general the weighted fst is preferred for multilocus calculations as the mean fst underestimates the true value. (monomorphic loci have fst=0 and should not form part of the average since they are not SNPs)
#GLOBAL STATISTICS
fst_raw <- read.table("~/NEC_WORK/variants/merged/no_mito/global_fst_ss_GQ13_70pCall.windowed.weir.fst", header=T)
fst <- fst_raw %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% gather("method","fst",WEIGHTED_FST:MEAN_FST)

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/global_fst_ss_GQ13_70pCall.png", res=300, width = 6, height = 6)

ggplot(fst)+
  geom_density(aes(x=fst, fill=method), alpha=.5)+
  geom_vline(xintercept = mean(subset(fst, method=="MEAN_FST")$fst), colour=twoggcols[1], lty=2)+
  geom_vline(xintercept = mean(subset(fst, method=="WEIGHTED_FST")$fst), colour=twoggcols[2], lty=2)+
  labs(x=expression(italic(F[ST])~"(10 Kb windows)"), y="Density")+
  scale_x_continuous(breaks = c(0,signif(mean(subset(fst, method=="MEAN_FST")$fst),3),0.25,signif(mean(subset(fst, method=="WEIGHTED_FST")$fst),3), 0.5, 0.75))+
  scale_fill_discrete(label=c("Unweighted Fst","Weighted Fst"))+
  theme(legend.title = element_blank(), legend.position = c(0.75, 0.75))

dev.off()

#PAIRWISE COMPARISONS
SYvsWH_fst <- read.table("~/NEC_WORK/variants/merged/no_mito/pwfst_SYvsWH_ss_GQ13_70pCall.windowed.weir.fst", header=T) %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(comp="Sylt vs. Woods Hole")
SYvsMI_fst <- read.table("~/NEC_WORK/variants/merged/no_mito/pwfst_SYvsMI_ss_GQ13_70pCall.windowed.weir.fst", header=T) %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(comp="Sylt vs. Miami")
SYvsVA_fst <- read.table("~/NEC_WORK/variants/merged/no_mito/pwfst_SYvsVA_ss_GQ13_70pCall.windowed.weir.fst", header=T) %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(comp="Sylt vs. Varna")
SYvsVF_fst <- read.table("~/NEC_WORK/variants/merged/no_mito/pwfst_SYvsVF_ss_GQ13_70pCall.windowed.weir.fst", header=T) %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(comp="Sylt vs. Villefranche")
WHvsMI_fst <- read.table("~/NEC_WORK/variants/merged/no_mito/pwfst_WHvsMI_ss_GQ13_70pCall.windowed.weir.fst", header=T) %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(comp="Woods Hole vs. Miami")
WHvsVA_fst <- read.table("~/NEC_WORK/variants/merged/no_mito/pwfst_WHvsVA_ss_GQ13_70pCall.windowed.weir.fst", header=T) %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(comp="Woods Hole vs. Varna")
WHvsVF_fst <- read.table("~/NEC_WORK/variants/merged/no_mito/pwfst_WHvsVF_ss_GQ13_70pCall.windowed.weir.fst", header=T) %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(comp="Woods Hole vs. Villefranche")
MIvsVA_fst <- read.table("~/NEC_WORK/variants/merged/no_mito/pwfst_MIvsVA_ss_GQ13_70pCall.windowed.weir.fst", header=T) %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(comp="Miami vs. Varna")
MIvsVF_fst <- read.table("~/NEC_WORK/variants/merged/no_mito/pwfst_MIvsVF_ss_GQ13_70pCall.windowed.weir.fst", header=T) %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(comp="Miami vs. Villefranche")
VAvsVF_fst <- read.table("~/NEC_WORK/variants/merged/no_mito/pwfst_VAvsVF_ss_GQ13_70pCall.windowed.weir.fst", header=T) %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(comp="Varna vs. Villefranche")

#Here we can calculate the absolute pairwise Fst values
pw_fst <- bind_rows(SYvsWH_fst, SYvsMI_fst, SYvsVA_fst, SYvsVF_fst, WHvsMI_fst, WHvsVA_fst, WHvsVF_fst, MIvsVA_fst, MIvsVF_fst, VAvsVF_fst)
summary_pw_fst <- pw_fst %>% group_by(comp) %>% summarise(mean_weighted= mean(WEIGHTED_FST), mean_unweighted=mean(MEAN_FST))
summary_pw_fst[,2:3] <- round(summary_pw_fst[,2:3], 3)

write.table(summary_pw_fst, "/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/pwfst_summary_ss_GQ13_70pCall.txt", sep="\t", row.names = F)

###Fst MANHATTAN PLOTS
#We can also plot the individual Fst windows by giving them an ID
#The equivalent of the PCAdapt 1. PC would be the comparison between north and south
#THIS IS PROBLEMATIC SINCE WE HAVE STRONG SUBSTRUCTURE HERE! DOESNT REALLY WORK FOR FST!
northvssouth_fst <- read.table("~/NEC_WORK/variants/merged/no_mito/pwfst_northvssouth_ss_GQ13_70pCall.windowed.weir.fst", header=T) %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(comp="North vs. South")
northvssouth_fst <- mutate(northvssouth_fst, BIN_ID=seq(1:nrow(northvssouth_fst)))

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/northvssouth_fst_manha_ss_GQ13_70pCall.png", res=300, width = 12, height = 6)

ggplot(northvssouth_fst)+
  geom_point(aes(y=WEIGHTED_FST, x=BIN_ID), color=threeggcols[1] ,alpha=.4)+
  labs(x="Window ID", y=expression("Weighted"~italic(F[ST])~"(10 Kb Windows)"))

dev.off()

#The equivalent of the PCAdapt 2. PC would be the comparison between MI and the southern invasives
MIvssouthinv_fst <- read.table("~/NEC_WORK/variants/merged/no_mito/pwfst_MIvssouthinv_ss_GQ13_70pCall.windowed.weir.fst", header=T) %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(comp="Miami vs. Southern Invasives")
MIvssouthinv_fst <- mutate(MIvssouthinv_fst, BIN_ID=seq(1:nrow(MIvssouthinv_fst)))

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/MIvssouthinv_fst_manha_ss_GQ13_70pCall.png", res=300, width = 12, height = 6)

ggplot(MIvssouthinv_fst)+
  geom_point(aes(y=WEIGHTED_FST, x=BIN_ID), color=threeggcols[2] ,alpha=.4)+
  labs(x="Window ID", y=expression("Weighted"~italic(F[ST])~"(10 Kb Windows)"))

dev.off()

#The equivalent of the PCAdapt 3. PC would be the comparison of Varna and Villefranche'
VAvsVF_fst <- mutate(VAvsVF_fst, BIN_ID=seq(1:nrow(VAvsVF_fst)))

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/VAvsVF_fst_manha_ss_GQ13_70pCall.png", res=300, width = 12, height = 6)

ggplot(VAvsVF_fst)+
  geom_point(aes(y=WEIGHTED_FST, x=BIN_ID), color=threeggcols[3] ,alpha=.4)+
  labs(x="Window ID", y=expression("Weighted"~italic(F[ST])~"(10 Kb Windows)"))

dev.off()

#Can also try natives vs invasives
nativeVSinvasive_fst <- read.table("~/NEC_WORK/variants/merged/no_mito/pwfst_nativeVSinvasive_ss_GQ13_70pCall.windowed.weir.fst", header=T) %>% left_join(glob_coord, by="CHROM") %>% dplyr::filter(length-BIN_END>=0) %>% mutate(comp="Native vs. Invasive")
nativeVSinvasive_fst <- mutate(nativeVSinvasive_fst, BIN_ID=seq(1:nrow(nativeVSinvasive_fst)))

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/nativeVSinvasive_fst_manha_ss_GQ13_70pCall.png", res=300, width = 12, height = 6)

ggplot(nativeVSinvasive_fst)+
  geom_point(aes(y=WEIGHTED_FST, x=BIN_ID), color="black" ,alpha=.4)+
  labs(x="Window ID", y=expression("Weighted"~italic(F[ST])~"(10 Kb Windows)"))

dev.off()

########FIXED RECIPROCAL SITES##########
#output allele frequencies for each population by using vcftools --freq and import this here
north_freq_raw <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/north_freq_ss_GQ13_70pCallinPop.frq", col.names = c("CHROM", "POS", "N_ALLELES", "N_CHR", paste("freq_",1:5, sep="")), header=F, skip=1)
south_freq_raw <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/south_freq_ss_GQ13_70pCallinPop.frq", col.names = c("CHROM", "POS", "N_ALLELES", "N_CHR", paste("freq_",1:5, sep="")), header=F, skip=1)

#filter for fixed alleles in each set
north_freq_fixed <- north_freq_raw %>%
                        gather("allele_no_north","freq_north",freq_1:freq_5) %>%
                        dplyr::filter(grepl(':1',freq_north))

south_freq_fixed <- south_freq_raw %>%
                        gather("allele_no_south","freq_south",freq_1:freq_5) %>%
                        dplyr::filter(grepl(':1',freq_south))

#join sites present in both sets and filter sites where different alleles are fixed
reciprocal_freq_fixed <- inner_join(north_freq_fixed, south_freq_fixed, by=c("CHROM","POS")) %>%
                            filter(freq_north != freq_south)