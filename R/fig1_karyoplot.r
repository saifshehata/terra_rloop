# set working directory.
setwd("/cluster/khiom/sshehata001/proj/terra_rloop/results/figures/Rplots")

# load packages.
library("data.table")
library("ggplot2")
library("dplyr")
# if not installed, install karyoploteR with: 
#BiocManager::install("karyoploteR", update=FALSE)
library("karyoploteR")

# read data
terra_intersect <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/peak_lengths/terra_intersect_peak_lengths.bed", header = F,sep="\t")


# plot karyo.
pdf("fig1_karyoplot.pdf", width = 5.6, height = 3.15)

mm10karyo <- plotKaryotype(genome = "mm10")
kpPlotRegions(mm10karyo,terra_intersect, col = "red")
#56B4E9 #E69F00

dev.off()
