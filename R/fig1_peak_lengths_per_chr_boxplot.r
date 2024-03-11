# set working directory.
setwd("/cluster/khiom/sshehata001/proj/terra_rloop/results/figures/Rplots")

# load packages.
library("data.table")
library("ggplot2")
library("dplyr")

# read bed files with peak coordinates.
terra_all <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/peak_lengths/terra_peaks_peak_lengths.bed", header = F,sep="\t")
rloop_all <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/peak_lengths/rloop_peaks_peak_lengths.bed", header = F,sep="\t")
terra_no_intersect <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/peak_lengths/terra_no_intersect_peak_lengths.bed", header = F,sep="\t")
rloop_no_intersect <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/peak_lengths/rloop_no_intersect_peak_lengths.bed", header = F,sep="\t")
terra_intersect <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/peak_lengths/terra_intersect_peak_lengths.bed", header = F,sep="\t")
rloop_intersect <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/peak_lengths/rloop_intersect_peak_lengths.bed", header = F,sep="\t")

# add new column to existing data frames to be able to create plots.
terra_all$Sample="TERRA"
terra_all$Peak_Group="All"

terra_no_intersect$Sample="TERRA"
terra_no_intersect$Peak_Group="Non-Intersecting"

terra_intersect$Sample="TERRA"
terra_intersect$Peak_Group="Intersecting"

rloop_all$Sample="R-loop"
rloop_all$Peak_Group="All"

rloop_no_intersect$Sample="R-loop"
rloop_no_intersect$Peak_Group="Non-Intersecting"

rloop_intersect$Sample="R-loop"
rloop_intersect$Peak_Group="Intersecting"

# combine all data frames.
combined_all <- rbind(terra_all, terra_no_intersect, terra_intersect, rloop_all, rloop_no_intersect, rloop_intersect)

# define color-blind palettes.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##################################################################################################################
# boxplots of peak length distribution among all_peaks, overlapping_peaks and non-overlapping_peaks 
# expanding upon the above boxplot to compare peak length distr. among individual chromosomes between 
# terra and rloop.
##################################################################################################################

# for Overlapping Peaks.
pdf("fig1_peak_lengths_per_chr_boxplot.pdf", width = 8, height = 4.5)
combined_all %>% filter(Peak_Group == "Intersecting") %>% 
  #filter(!grepl("(random|chrUn|chrM)",V1, perl = T)) %>% 
  #ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = factor(Sample, levels = unique(Sample)))) + 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Sample)) + 
  geom_boxplot(outlier.size = 0) + 
  xlab("") + ylab("Peak Length (bp)") + 
  scale_y_continuous(limits = c(100,3000), breaks=seq(0,3000, by=500)) + 
  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.title = element_blank(), legend.position = "top", text = element_text(size = 22)) + 
  facet_wrap(~Peak_Group) + 
  #guides(fill = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values=cbPalette)
dev.off()
