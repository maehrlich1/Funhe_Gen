# Assumes you've already run bedtools -genomecov and grepped "genome". E.g. something like:
# find *.bam | parallel 'bedtools -genomecov -ibam {} -g genome.fa | grep 'genome' > {}_cov_genome.txt'

# Get a list of the bedtools output files you'd like to read in
file_names <- list.files(pattern="*genome.txt")
file_names_TS <- grep(file_names, pattern="*TS*", value=T)
file_names_NX <- grep(file_names, pattern="*TS*", value=T, inv=T)

# Optional, create short sample names from the filenames. 
# For example, in this experiment, my sample filenames might look like this:
# prefixToTrash-01.pe.on.pos.dedup.realigned.recalibrated.bam
# prefixToTrash-02.pe.on.pos.dedup.realigned.recalibrated.bam
# prefixToTrash-03.pe.on.pos.dedup.realigned.recalibrated.bam
# This regular expression leaves me with "samp01", "samp02", and "samp03" in the legend.
#print(labs <- paste("samp", gsub("prefixToTrash-0|\\.pe\\.on\\.pos\\.dedup\\.realigned\\.recalibrated\\.bam\\.cov\\.hist\\.txt\\.all\\.txt", "", file_names, perl=TRUE), sep=""))

# Create lists to hold coverage and cumulative coverage for each alignment,
# and read the data into these lists. Creating a vecor containing column classes can help structure the data frames.
CC <- c("factor","numeric","numeric","numeric","numeric")
cov_TS <- list()
cov_NX <- list()
#cov_cumul <- list()
for (i in 1:length(file_names_TS)) {
  cov_TS[[i]] <- read.table(file_names_TS[i], header=F, colClasses=CC)
  #cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
}
for (i in 1:length(file_names_NX)) {
  cov_NX[[i]] <- read.table(file_names_NX[i], header=F, colClasses=CC)
  #cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
}

# Pick some colors
# Ugly:
#cols <- 1:length(cov)
# Prettier:
# ?colorRampPalette
# display.brewer.all()
#library(RColorBrewer)
#cols <- brewer.pal(length(cov), "Dark2")
#If a lot of (>8) differentiable colours are needed:
cols <- c("dodgerblue2","#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2","#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")

# Save the graph to a file.
#png("coverage_hist_TS_vs_NX.png", h=1000, w=1500, pointsize=20)

# Create plot area, but do not plot anything. Add gridlines and axis labels. Can also plot the cumulative stats by adding cov_cumul[[1]][1:400].
plot(cov_TS[[1]][,5]~cov_TS[[1]][,2], type='n', xlab="Depth", xlim=c(0,4), ylab="Fraction of bases", ylim=c(0,1.0), main="Genome Coverage")
abline(v = 20, col = "gray60")
abline(v = 50, col = "gray60")
abline(h = 0.50, col = "gray60")
abline(h = 0.90, col = "gray60")
axis(1, at=c(20,50,80), labels=c(20,50,80))
axis(2, at=c(0.50,0.90), labels=c(0.50,0.90))

# Actually plot the data for each of the alignments (stored in the lists).
for (i in 1:length(cov_TS)) points(cov_TS[[i]][,5]~cov_TS[[i]][,2], type='l', lwd=3, col="dodgerblue2")
for (i in 1:length(cov_NX)) points(cov_NX[[i]][,5]~cov_NX[[i]][,2], type='l', lwd=3, col="orange")

# Add a legend.
legend("topright", legend=c("TruSeq","Nextera"), col=c("dodgerblue2","orange"), cex=1.5, lty=1, lwd=4)

#dev.off()