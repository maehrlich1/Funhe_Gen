setwd("/Users/Moritz/Documents/Academic/GEOMAR/Thesis/Seq/miseq/mito_var/cons_seq")

library(pegas)
library(RColorBrewer)

input <- "all_mt_rotated.maf"
d <- ape::read.dna(input, format='fasta')
e <- dist.dna(d)
h <- pegas::haplotype(d)
h <- sort(h, what = "label")
net <- pegas::haploNet(h)
ind.hap <- with(stack(setNames(attr(h, "index"), rownames(h))),
                table(hap=ind, pop=rownames(d)[values]))

png(filename="~/Desktop/HapNet_mt.png", width=800, height=600)

plot(net, size=attr(net, "freq"), scale.ratio=0.2, pie=ind.hap)
legend(6, -3, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=2)

dev.off()