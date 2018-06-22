####### IN SILICO DIGEST #################
library(SimRAD)
library(ggplot2)
library(dplyr)
library(stringr)
myPng <- function(..., width=6, height=6, res=300) {png(..., width=width*res, height=height*res, res=res)}

##get FASTA
#genome reduction factor
ratio <- 1
seq <- ref.DNAseq("~/Documents/Academic/RSMAS/PhD/Seq/reference/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fna", subselect.contigs = T, prop.contigs = ratio)

###########################
###Restriction enzymes#####
###########################
###TA Overhang#########
###4-Base Cutters######

#MseI family   5'  T'TA A  3'
#              3'  A AT'T  5'
MseI_5p <- "T"
MseI_3p <- "TAA"

#MaeI family  5'  C'TA G  3'
#             3'  G AT'C  5'
MaeI_5p <- "C"
MaeI_3p <- "TAG"

#CviQI family   5'  G'TA C  3'
#               3'  C AT'G  5'
CviQI_5p <- "G"
CviQI_3p <- "TAC"

###6-Base Cutters####
#AseI family  5'  AT'TA AT  3'
#             3'  TA AT'TA  5'
AseI_5p <- "AT"
AseI_3p <- "TAAT"

#NdeI family  5'  CA'TA TG  3'
#             3'  GT AT'AC  5'
NdeI_5p <- "CA"
NdeI_3p <- "TATG"

#####CG Overhang######
###4-base cutters###

#HinP1I family   5'  G'CG C  3'
#                3'  C GC'G  5'
HinP1I_5p <- "G"
HinP1I_3p <- "CGC"

#HpaII family   5'  C'CG G  3'
#               3'  G GC'C  5'
HpaII_5p <- "C"
HpaII_3p <- "CGG"

#MaeII family   5'  A'CG T  3'
#               3'  T GC'A  5'
MaeII_5p <- "A"
MaeII_3p <- "CGT"

#TaqI family  5'  T'CG A  3'
#             3'  A GC'T  5'
TaqI_5p <- "T"
TaqI_3p <- "CGA"

###6-base cutters######

#AclI family  5'  AA'CG TT  3'
#             3'  TT GC'AA  5'
AclI_5p <- "AA"
AclI_3p <- "CGTT"

#AcyI family  5'  GR'CG YC  3'
#             3'  CY GC'RG  5'
AcyI_5p <- "GR"
AcyI_3p <- "CGYC"

#AsuII family   5'  TT'CG AA  3'
#               3'  AA GC'TT  5'
AsuII_5p <- "TT"
AsuII_3p <- "CGAA"

#ClaI family  5'  AT'CG AT  3'
#             3'  TA GC'TA  5'
ClaI_5p <- "AT"
ClaI_3p <- "CGAT"

#NarI family  5'  GG'CG CC  3'
#             3'  CC GC'GG  5'
NarI_5p <- "GG"
NarI_3p <- "CGCC"

####Plug-in enzyme to run
cs_5p <- AseI_5p
cs_3p <- AseI_3p
deg <- 6 #plug-in the degeneracy of the enzyme i.e. 4-base cutter
min_size <- 250 #minimum fragment size
max_size <- 480

###In-silico digest

digest <- insilico.digest(seq, cs_5p, cs_3p, verbose = T)
#can call the length of all fragments in the output by using nchar() function

#Get only fragments that have flanking restriction sites (otherwise no overhang)
#digest <- adapt.select(digest,type = "AA", MseI_5p, MseI_3p)
#not really necessary when using only 1 enzyme

#This part allows us to exclude fragments that have many N bases.
#We assume that an x-base cutter cuts every x^4 bases. To approximate, we exclude any fragments that have more than x^4 N's in their sequence.
digest_noN <- subset(digest, str_count(digest, "N") < 4^deg)

### select a certain insert size for library prep
size_sel <- size.select(digest_noN,  min.size = min_size, max.size = max_size, verbose = T, graph = F)
#Illumina systems work best with insert sizes of 300-400 bp. Check your specific system!

###Graph
#make dataframe for ggplot
fl_all <- cbind.data.frame(frag_length=nchar(digest_noN))
fl_ss <- cbind.data.frame(frag_length=nchar(size_sel))

myPng("/Users/Moritz/Documents/Academic/RSMAS/PhD/GBS/insilico_digest/TA_overhang/6_base_cutters/AseI/AseI_250-480.png", res=300, width = 8, height = 6)

ggplot()+
  geom_histogram(data=fl_all, aes(x=frag_length), stat = 'count')+
  geom_histogram(data=fl_ss, aes(x=frag_length), stat= 'count', fill='coral')+
  labs(x="Fragment Length", y="Count")+
  annotate("text", x=3000, y=500, label=paste(nrow(fl_ss), "fragments between", min_size, "and", max_size, "bp."), col="coral")+
  coord_cartesian(x=c(0,quantile(fl_all$frag_length, 0.95)), y=c(0,1000))

dev.off()


#####################################
#####PCR with Extended Primers#########
######################################
#Can further subsample by using primers that only amplify fragments that begin/end with a certain base sequence

###1-BP Extensions#######
#Define base extensions to try
bext <- c("A","T","C","G")

#initialise output list to write frags to
ext_sel <- vector("list",length(bext))
names(ext_sel) <- bext

for(i in 1:length(bext))
{
ext_sel[[i]] <- grep(paste(bext[i],cs_5p,"$", sep = ""), size_sel, value=T)
}

#get the number of fragments by running the length function over the entire list
lapply(ext_sel,length)

#Plot the size distribution of your selected fragments
#first get the frag sizes for your specific bp extension
d <- cbind.data.frame(frag_length=width(ext_sel[["C"]]))

myPng("/Users/Moritz/Documents/Academic/RSMAS/PhD/GBS/insilico_digest/TA_overhang/6_base_cutters/AseI/AseI_250-480_C_ext.png", res=300, width = 8, height = 6)

ggplot()+
  geom_histogram(data=d, aes(x=frag_length), binwidth = 1)+
  labs(x="Fragment Length", y="Count")+
  annotate("text", x=mean(c(max_size,min_size)), y=110, label=paste(nrow(d), "fragments"), col="coral")

dev.off()

###2-BP Extension###
####For 2-bp extensions you can double-loop through the bext vector.
bext <- c("A","T","C","G")

#first we need a list in the form of a matrix, use the dim function for this
ext_sel <- vector("list", length(bext)^2)
dim(ext_sel) <- c(rep(length(bext),2))
#this creates a two-dimensional list. now just rename the row and column headers.
rownames(ext_sel) <- bext
colnames(ext_sel) <- bext

#now we can start adding data to it

for(i in 1:length(bext))
{
  for(j in 1:length(bext))
  {
    ext_sel[[i,j]] <- grep(paste(bext[j],bext[i],cs_5p,"$", sep = ""), size_sel, value=T)
  }
}

print(ext_sel)

#plot the size distribution of your selected fragments
#first get the frag sizes for your specific bp extension
d <- cbind.data.frame(frag_length=width(ext_sel[["T","T"]]))

myPng("/Users/Moritz/Documents/Academic/RSMAS/PhD/GBS/insilico_digest/TA_overhang/4_base_cutters/MaeI/MaeI_250-500_TT_ext.png", res=300, width = 8, height = 6)

ggplot()+
  geom_histogram(data=d, aes(x=frag_length), binwidth = 1)+
  labs(x="Fragment Length", y="Count")+
  annotate("text", x=375, y=900, label=paste(nrow(d), "fragments"), col="coral")

dev.off()
