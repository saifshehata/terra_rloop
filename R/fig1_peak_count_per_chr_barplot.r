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



###########################################################################################################################
# bar chart of peak distribution per chromosome (nr of peaks per chromosome)
###########################################################################################################################

pdf("fig1_peak_count_per_chr_barplot.pdf", width = 8, height = 4.5)
combined_all %>% filter(Peak_Group != "All") %>% 
  #filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(y = factor(V1, levels = unique(terra_all$V1)), fill = Sample)) + 
  #ggplot(.,aes(y = factor(V1, levels = unique(terra_all$V1)), fill = factor(Sample, levels = unique(Sample)))) + 
  geom_bar(position=position_dodge()) +
  facet_grid(factor(Peak_Group, levels = unique(Peak_Group)) ~., scales = "free") + 
  #facet_grid(Peak_Group ~., scales = "free") + 
  xlab("Nr. of Peaks") + ylab("") + 
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.title = element_blank(), legend.position = "top", text = element_text(size = 22)) + 
  scale_fill_manual(values=cbPalette)
dev.off()

