# set working directory.
setwd("/proj/nb_storage/private/terra_rloop_project/results/figures/karyoploteR/")

# load packages.
library("data.table")
library("ggplot2")
library("dplyr")
#install karyoploteR with: BiocManager::install("karyoploteR", update=FALSE)
library("karyoploteR")

# read data
terra_all <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_peaks_peakLength.bed", header = F,sep="\t")
rloop_all <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/rloop_peaks_peakLength.bed", header = F,sep="\t")
terra_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_overlap_peakLength.bed", header = F,sep="\t")
rloop_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/rloop_overlap_peakLength.bed", header = F,sep="\t")
terra_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_noOverlap_peakLength.bed", header = F,sep="\t")
rloop_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/rloop_noOverlap_peakLength.bed", header = F,sep="\t")

# plot karyo.
pdf("Fig1C_peakDistr_karyoplot.pdf", width = 8, height = 4.5)

mm10karyo<-plotKaryotype(genome = "mm10")
kpPlotRegions(mm10karyo,terra_overlap,col = "red")
#D55E00, #0072B2
dev.off()

# mm10karyo<-plotKaryotype(genome = "mm10")
# kpPlotRegions(mm10karyo,terra_overlap,col = "blue")
# kpPlotRegions(mm10karyo,rloop_all,col = "yellow")

