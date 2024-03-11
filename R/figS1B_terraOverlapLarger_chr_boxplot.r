# set working directory.
setwd("/proj/nb_storage/private/terra_rloop_project/results/figures/")

# load packages.
library("data.table")
library("ggplot2")
library("dplyr")

# read bed files with peak coordinates.
terra_all <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_peaks_peakLength.bed", header = F,sep="\t")
rloop_all <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/rloop_peaks_peakLength.bed", header = F,sep="\t")
terra_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_overlap_peakLength.bed", header = F,sep="\t")
rloop_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/rloop_overlap_peakLength.bed", header = F,sep="\t")
terra_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_noOverlap_peakLength.bed", header = F,sep="\t")
rloop_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/rloop_noOverlap_peakLength.bed", header = F,sep="\t")

# add new column to existing data frames to be able to create boxplots.
terra_all$Sample="Terra"
terra_all$Peak_Group="All"

terra_overlap$Sample="Terra"
terra_overlap$Peak_Group="Overlapping"

terra_noOverlap$Sample="Terra"
terra_noOverlap$Peak_Group="Non-Overlapping"

rloop_all$Sample="R-loop"
rloop_all$Peak_Group="All"

rloop_overlap$Sample="R-loop"
rloop_overlap$Peak_Group="Overlapping"

rloop_noOverlap$Sample="R-loop"
rloop_noOverlap$Peak_Group="Non-Overlapping"

# combine all data frames.
combined_all <- rbind(terra_all,terra_overlap,terra_noOverlap,rloop_all,rloop_overlap,rloop_noOverlap)

# define color-blind palettes.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



##################################################################################################################
# boxplot of peak length distribution showing TERRA overlapping peaks are generally smaller than non-overlapping
##################################################################################################################

pdf("FigS1B_chr_terraOverlapLarger_boxplot.pdf", width = 8, height = 4.5)
combined_all %>% filter(Peak_Group != "All" & Sample == "Terra") %>% 
  filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Peak_Group)) + 
  geom_boxplot(outlier.size = 0) + 
  xlab("") + ylab("Peak Length (bp)") + 
  scale_y_continuous(limits = c(100,3000), breaks=seq(0,3000, by=500)) + 
  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.title = element_blank(), legend.position = "top") + 
  facet_wrap(~Sample) +
  scale_fill_manual(values=cbPalette)
dev.off()

