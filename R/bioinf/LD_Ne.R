#########LINKAGE DISEQUILIBRIUM AND ESTIMATION OF EFFECTIVE POPULATION SIZE##########
library(ggplot2)
library(dplyr)
myPng <- function(..., width=6, height=6, res=300) {png(..., width=width*res, height=height*res, res=res)}
pops <- c(rep("Miami",15), rep("Sylt",16), rep("Varna",16), rep("Villefranche",9), rep("Woods Hole",16))
sample_size <- cbind.data.frame(Population=c("Sylt","Woods Hole","Miami","Varna","Villefranche"), sample_size=c(16,16,15,16,9))
cols <- c("#E41A1C","#984EA3","#FF7F00","#377EB8","#4DAF4A")

#Import raw data. The ZOOM files were a bit large and hat to be downsampled by a factor 10.
SY_ld_wide <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/SY_ld_WIDE_ss_GQ13_1miss_maf5p_bial.geno.ld") %>% mutate(Population="Sylt")
#SY_ld_zoom <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/SY_ld_ZOOM_ss_GQ13_1miss_maf5p_bial.geno.ld") %>% mutate(Population="Sylt")

WH_ld_wide <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/WH_ld_WIDE_ss_GQ13_1miss_maf5p_bial.geno.ld") %>% mutate(Population="Woods Hole")
#WH_ld_zoom <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/WH_ld_ZOOM_ss_GQ13_1miss_maf5p_bial.geno.ld") %>% mutate(Population="Woods Hole")

MI_ld_wide <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/MI_ld_WIDE_ss_GQ13_1miss_maf5p_bial.geno.ld") %>% mutate(Population="Miami")
#MI_ld_zoom <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/MI_ld_ZOOM_ss_GQ13_1miss_maf5p_bial.geno.ld") %>% mutate(Population="Miami")

VA_ld_wide <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VA_ld_WIDE_ss_GQ13_1miss_maf5p_bial.geno.ld") %>% mutate(Population="Varna")
#VA_ld_zoom <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VA_ld_ZOOM_ss_GQ13_1miss_maf5p_bial.geno.ld") %>% mutate(Population="Varna")

VF_ld_wide <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VF_ld_WIDE_ss_GQ13_1miss_maf5p_bial.geno.ld") %>% mutate(Population="Villefranche")
#VF_ld_zoom <- read.delim("/Users/Moritz/NEC_WORK/variants/merged/no_mito/VF_ld_ZOOM_ss_GQ13_1miss_maf5p_bial.geno.ld") %>% mutate(Population="Villefranche")

#merge data and calculate distances
pop_ld_wide <- bind_rows(SY_ld_wide, WH_ld_wide, MI_ld_wide, VA_ld_wide, VF_ld_wide) %>% mutate(dist=POS2-POS1)
#pop_ld_zoom <- bind_rows(SY_ld_zoom, WH_ld_zoom, MI_ld_zoom, VA_ld_zoom, VF_ld_zoom) %>% mutate(dist=POS2-POS1)

###Corbin (2012) method for Ne estimation using unphased markers
#Parameter definition
bw_wide=1000      #this determines over which distance width r2 values will be averaged (should be adjusted depending on resolution of SNPs tested and the scale)
res_wide=500      #this is the snp resolution in the dataset. it is used to shift the bins upwards as there can be no data below this value.
bw_zoom=100
res_zoom=50
alpha=2           #this takes values between 1-2.2 and is a correction for mutations (Ohta and Kimura 1971). Essentially a lower number means a lower mutation rate and 2.2 a high rate. A low value fits short-term data better 20-200 generations (less mut.), higher values are more accurate in the long term. Tenesa 2007 suggest a=2. can calculate this empirically by considering very short distances only. there the r2 must be purely due to mutation and hence a = 1/r2
beta=1            #this takes either the value of 1 or 2 depending on whether the data is unphased or phased respectively
recrate=3e-08     #the recombination rate, a value of 1e-08 is common. Here I took the value for C. elegans (3 cM/Mb, Prachumwat 2003)
mat_time=21/365   #the time to maturity in years (used to convert generation times to years) (for Mleidyi Baker & Reeve 1974)

#binning is done by taking the dist value and substracting the residual of (dist-resolution) divided by the binwidth then centering the bin by adding half of the binwidth e.g. 5347-297+500=5550 this bin will contain values between 5050 and 6050. see where the res parameter comes in?
#next we group by population and bin and take the weighted average of r2 values. we need to weight them since we have missing data and some r2 values are therefore worth more.
#instead of using the sample size we also have to calculate a harmonic mean of sample size for the same reason as above
#linkage distance is calculated and corrected according to Sved & Feldman 1973 (different methods! check SNeP software manual/paper)
#generations in the past are calculated according to Hayes 2003

pop_ld_wide_bin <- pop_ld_wide %>%
                        mutate(dist_bin = dist - ((dist-res_wide) %% bw_wide) + bw_wide/2) %>%   #binning
                        group_by(Population, dist_bin) %>%
                        summarise(mean_r2 = sum(N_INDV*R.2)/sum(N_INDV), harm_N = 1/mean(1/N_INDV)) %>%  #calculate r2 weighted mean and harmonic mean of N
                        mutate(r_adj = mean_r2 - (1/(beta*harm_N)), d = dist_bin*recrate) %>%   #correct r2 for sample size and calculate linkage distance, d
                        mutate(c = 0.5*(1-exp(-2*d))) %>% #Haldanes mapping function which converts linkage distance into recombination frequency, c
                        mutate(corr_c = c*(1-(c/2))) %>%   #correct recombination frequency according to (Sved & Feldman 1973)
                        mutate(gen = 1/(2*corr_c)) %>%   #calculate generations from linkage distance. ASSUMES CONSTANT OR LINEAR GROWTH! (Hayes 2003)
                        mutate(Ne = (gen/2)*((1/r_adj)-alpha)+0.5) %>%    #the formula for Ne from Corbin 2012
                        mutate(yr = 2016-(gen*mat_time))   #converting generations to years using time to maturity

pop_ld_zoom_bin <- pop_ld_zoom %>%
                        mutate(dist_bin = dist - ((dist-res_zoom) %% bw_zoom) + bw_zoom/2) %>%   #binning
                        group_by(Population, dist_bin) %>%
                        summarise(mean_r2 = sum(N_INDV*R.2)/sum(N_INDV), harm_N = 1/mean(1/N_INDV)) %>%  #calculate r2 weighted mean and harmonic mean of N
                        mutate(r_adj = mean_r2 - (1/(beta*harm_N)), d = dist_bin*recrate) %>%   #correct r2 for sample size and calculate linkage distance, d
                        mutate(c = 0.5*(1-exp(-2*d))) %>% #Haldanes mapping function which converts linkage distance into recombination frequency, c
                        mutate(corr_c = c*(1-(c/2))) %>%   #correct recombination frequency according to (Sved & Feldman 1973)
                        mutate(gen = 1/(2*corr_c)) %>%   #calculate generations from rec. frequency. ASSUMES CONSTANT OR LINEAR GROWTH! (Hayes 2003)
                        mutate(Ne = (gen/2)*((1/r_adj)-alpha)+0.5) %>%    #the formula for Ne from Corbin 2012
                        mutate(yr = 2016-(gen*mat_time))   #converting generations to years using time to maturity

#plot the LD decay plots
#first the WIDE long-range LD

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/pop_ld_WIDE_ss_GQ13_1miss_maf5p_bial.png", res=300, width = 8, height = 6)

ggplot(pop_ld_wide_bin)+
  geom_line(aes(x=dist_bin, y=r_adj, color=Population), size=1.1, alpha=.8)+
  labs(x="Physical Distance (Kb)", y=expression("Weighted Mean Correlation Coefficient,"~italic(r^2)))+
  scale_x_continuous(breaks=seq(0,400000,50000), labels=seq(0,400,50))+
  scale_color_manual(values=cols)+
  theme(legend.position = c(0.8,0.8))

dev.off()

#now ZOOM into the short-range LD

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/pop_ld_ZOOM_ss_GQ13_1miss_maf5p_bial.png", res=300, width = 8, height = 6)

ggplot(pop_ld_zoom_bin)+
  geom_line(aes(x=dist_bin, y=r_adj, color=Population), size=1.1, alpha=.8)+
  labs(x="Physical Distance (Kb)", y=expression("Weighted Mean Correlation Coefficient,"~italic(r^2)))+
  scale_x_continuous(breaks=seq(0,20000,5000), labels=seq(0,20,5))+
  scale_color_manual(values=cols)+
  theme(legend.position = c(0.8,0.8))

dev.off()

#Can also calculate half-decay distances by getting the half-max r_adj

half_decay <- pop_ld_wide_bin %>% group_by(Population) %>% summarise(half_decay_r = max(r_adj)/2)

####Finally plot Ne
#WIDE data is better for recent times (LRLD) but has lower resolution in deep time (SRLD). Dont really care that much about deep time anyway...

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/pop_Ne_GEN150_ss_GQ13_1miss_maf5p_bial.png", res=300, width = 8, height = 6)

ggplot(subset(pop_ld_wide_bin, gen<=150))+
  geom_line(aes(x=gen, y=Ne, color=Population), alpha=.5)+
  geom_smooth(aes(x=gen, y=Ne, color=Population), alpha=.8, se=F, span=.2)+
  labs(x="Generations Ago", y=expression("Estimated Effective Population Size,"~italic(N[e])))+
  scale_color_manual(values=cols)+
  theme(legend.position = c(0.15,0.85))+
  scale_x_continuous(breaks = seq(50,150,20))

dev.off()

#Can also plot Ne against Years by converting the generations using time to maturity

myPng("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Graphs/popgen_stats/pop_Ne_YEAR1980_ss_GQ13_1miss_maf5p_bial.png", res=300, width = 8, height = 6)

ggplot(subset(pop_ld_wide_bin, yr>=1980))+
  geom_line(aes(x=yr, y=Ne, color=Population), alpha=.5)+
  geom_smooth(aes(x=yr, y=Ne, color=Population), alpha=.8, se=F, span=.1)+
  labs(x="Year", y=expression("Estimated Effective Population Size,"~italic(N[e])))+
  scale_color_manual(values=cols)+
  theme(legend.position = c(0.9,0.85))+
  scale_x_continuous(breaks = seq(1980,2015,5))

dev.off()

#ZOOM data has no for recent past (LRLD) but a higher resolution in deep time (SRLD). Not really useful for us right now.
