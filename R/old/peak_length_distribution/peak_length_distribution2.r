# set working directory.
setwd("/proj/nb_storage/private/terra_rloop_project/results/figures/peak_length_distribution/")

# load packages.
library("data.table")
library("ggplot2")
library("dplyr")

# read data
terra_all <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_peaks_peakLength.bed", header = F,sep="\t")
rloop_all <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/rloop_peaks_peakLength.bed", header = F,sep="\t")
terra_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_overlap_peakLength.bed", header = F,sep="\t")
rloop_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/rloop_overlap_peakLength.bed", header = F,sep="\t")
terra_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_noOverlap_peakLength.bed", header = F,sep="\t")
rloop_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/rloop_noOverlap_peakLength.bed", header = F,sep="\t")

# add new column to existing data frames to be able to create boxplots.
terra_all$Sample="Terra"
terra_all$Peaks="All"

terra_overlap$Sample="Terra"
terra_overlap$Peaks="Overlapping"

terra_noOverlap$Sample="Terra"
terra_noOverlap$Peaks="Non-Overlapping"

rloop_all$Sample="R-loop"
rloop_all$Peaks="All"

rloop_overlap$Sample="R-loop"
rloop_overlap$Peaks="Overlapping"

rloop_noOverlap$Sample="R-loop"
rloop_noOverlap$Peaks="Non-Overlapping"

# combine all data frames.
combined_all <- rbind(terra_all,terra_overlap,terra_noOverlap,rloop_all,rloop_overlap,rloop_noOverlap)

# define color-blind palettes.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


##################################################################################################################
# boxplot: peak length distribution per peak group (divided by all_peaks, overlap and noOverlap), 
# and separated by sample (terra vs rloop)
##################################################################################################################

# wrap by peaks.
pdf("peakLengths_byPeaks_boxplot.pdf", width = 8.74, height = 5.03)
combined_all %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = Sample, y = V5, fill = Peaks)) + geom_boxplot() + 
  facet_wrap(~Sample, scales="free_x") + 
  #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x", strip.position = "bottom") + 
  xlab("Peaks") + ylab("Peak Length (bp)") + ggtitle("Peak Length Distribution") + ylim(0,3000) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), legend.title = element_blank()) + 
  scale_fill_manual(values=cbPalette)
dev.off()

# wrap by sample.
pdf("peakLengths_bySample_boxplot.pdf", width = 8.74, height = 5.03)
combined_all %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = Peaks, y = V5, fill = Sample)) + geom_boxplot() + 
  facet_wrap(~Peaks, scales="free_x") + 
  
  xlab("Peaks") + ylab("Peak Length (bp)") + ggtitle("Peaks Length Distribution") + ylim(0,3000) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), legend.title = element_blank()) + 
  scale_fill_manual(values=cbPalette)
dev.off()


##################################################################################################################
# boxplots of peak length distribution among all_peaks, overlapping_peaks and non-overlapping_peaks 
# expanding upon the above boxplot to compare peak length distr. among individual chromosomes between 
# terra and rloop.
##################################################################################################################

# for All Peaks.
pdf("peakLengths_chr_allPeaks_boxplot.pdf", width = 8.74, height = 5.03)
combined_all %>% filter(Peaks == "All") %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Sample)) + 
  geom_boxplot(position = "dodge") + 
  xlab("Peaks by Chromosome") + ylab("Peak Length (bp)") + ggtitle("Peak Length Distribution") + ylim(0,3000) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "right", legend.title = element_blank()) + 
  facet_wrap(~Peaks) + 
  scale_fill_manual(values=cbPalette)
dev.off()

# for Non-Overlapping Peaks.
pdf("peakLengths_chr_noOverlap_boxplot.pdf", width = 8.74, height = 5.03)
combined_all%>% filter(Peaks =="Non-Overlapping") %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Sample)) + 
  geom_boxplot(position = "dodge") + 
  xlab("Peaks by Chromosome") + ylab("Peak Length (bp)") + ggtitle("Peak Length Distribution") + ylim(0,3000) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "right", legend.title = element_blank()) + 
  facet_wrap(~Peaks) + 
  scale_fill_manual(values=cbPalette)
dev.off()

# for Overlapping Peaks.
pdf("peakLengths_chr_overlap_boxplot.pdf", width = 8.74, height = 5.03)
combined_all %>% filter(Peaks == "Overlapping") %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Sample)) + 
  geom_boxplot(position = "dodge") + 
  xlab("Peaks by Chromosome") + ylab("Peak Length (bp)") + ggtitle("Peak Length Distribution") + ylim(0,3000) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.title = element_blank()) + 
  facet_wrap(~Peaks) + 
  #guides(fill = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values=cbPalette)
dev.off()


##################################################################################################################
# boxplot of peak length distribution showing R-LOOP overlapping peaks are generally smaller than non-overlapping
##################################################################################################################

pdf("peakLengths_chr_rloopOverlapSmall_boxplot.pdf", width = 8.74, height = 5.03)
combined_all %>% filter(Peaks != "All" & Sample == "R-loop") %>% 
  filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Peaks)) + 
  geom_boxplot() + 
  xlab("Peaks by Chromosome") + ylab("Peak Length (bp)") + ggtitle("Peak Length Distribution") + ylim(0,3000) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.title = element_blank()) + 
  facet_wrap(~Sample) +
  scale_fill_manual(values=cbPalette)
dev.off()


##################################################################################################################
# boxplot of peak length distribution showing TERRA overlapping peaks are generally smaller than non-overlapping
##################################################################################################################

pdf("peakLengths_chr_terraOverlapLarge_boxplot.pdf", width = 8.74, height = 5.03)
combined_all %>% filter(Peaks != "All" & Sample == "Terra") %>% 
  filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Peaks)) + 
  geom_boxplot() + 
  xlab("Peaks by Chromosome") + ylab("Peak Length (bp)") + ggtitle("Peak Length Distribution") + ylim(0,3000) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.title = element_blank()) + 
  facet_wrap(~Sample) +
  scale_fill_manual(values=cbPalette)
dev.off()


#########################################################################################################################################
# peak distribution histogram showing terra overlapping peaks are longer than rloop ones
# i.e. there is a larger amount of long terra peaks, and a larger amount of short rloop peaks
#########################################################################################################################################

pdf("peakLengths_histogram.pdf", width = 8.74, height = 5.03)
combined_all %>% filter(Peaks == "Overlapping") %>% ggplot(.,aes(x=V5, fill=Sample)) + 
  geom_histogram(binwidth=100, position = "dodge") + 
  scale_x_continuous(limits = c(0,1500), breaks=seq(0,1500, by=100)) + 
  xlab("Peak Length (bp)") + ylab("Peak Count") + ggtitle("R-loop Overlapping Peaks") + 
  facet_wrap(~Peaks) + 
  theme_bw() + 
  theme(legend.title = element_blank()) + 
  scale_fill_manual(values=cbPalette)
dev.off()


#########################################################################################################################################
# bar chart of peak distribution per chromosome (nr of peaks per chromosome)
#########################################################################################################################################

pdf("peakCount_chr.pdf", width = 8.74, height = 5.03)
combined_all %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(y = factor(V1, levels = unique(terra_all$V1)), fill = Sample)) + 
  geom_bar(position=position_dodge()) +
  facet_grid(Peaks ~., scales = "free") + 
  xlab("Peak Count") + ylab("Peaks per Chromosome") + ggtitle("Nr. Peaks per Chromosome") + 
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.title = element_blank()) + 
  scale_fill_manual(values=cbPalette)
dev.off()
