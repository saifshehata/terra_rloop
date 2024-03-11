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
combined_everything <- rbind(terra_all,terra_overlap,terra_noOverlap,rloop_all,rloop_overlap,rloop_noOverlap)
combined_allPeaks_only <- rbind(terra_all,rloop_all)
combined_noOverlap_only <- rbind(terra_noOverlap,rloop_noOverlap)
combined_overlap_only <- rbind(terra_overlap,rloop_overlap)
combined_terra <- rbind(terra_all,terra_noOverlap,terra_overlap)
combined_rloop <- rbind(rloop_all,rloop_noOverlap,rloop_overlap)
combined_overlap_noOverlap <- rbind(terra_overlap,terra_noOverlap,rloop_overlap,rloop_noOverlap)
combined_overlap_noOverlap_terra <- rbind(terra_overlap,terra_noOverlap)
combined_overlap_noOverlap_rloop <- rbind(rloop_overlap,rloop_noOverlap)

#########################################################################################################################################
# define color-blind palettes.
#########################################################################################################################################

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
#scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)

##################################################################################################################
# boxplot: peak length distribution per peak group (divided by all_peaks, overlap and noOverlap), 
# and separated by sample (terra vs rloop)
##################################################################################################################

# draft
combined_everything %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = Sample, y = V5, fill = Peaks)) + geom_boxplot() + 
  xlab("") + ylab("Peak Length") + ggtitle("Peaks Length Distribution") + ylim(0,3000) + 
  scale_fill_manual(values=cbPalette)

# final
pdf("peakLengths_byPeaks.pdf", width = 8.74, height = 5.03)
combined_everything %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = Sample, y = V5, fill = Peaks)) + geom_boxplot() + 
  facet_wrap(~Sample, scales="free_x") + 
  #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x", strip.position = "bottom") + 
  xlab("") + ylab("Peak Length") + ggtitle("Peaks Length Distribution") + ylim(0,3000) + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=cbPalette)
dev.off()


# draft
combined_everything %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = Peaks, y = V5, fill = Sample)) + geom_boxplot() + 
  xlab("") + ylab("Peak Length") + ggtitle("Peaks Length Distribution") + ylim(0,3000) + 
  scale_fill_manual(values=cbPalette)

# final
pdf("peakLengths_bySample.pdf", width = 8.74, height = 5.03)
combined_everything %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = Peaks, y = V5, fill = Sample)) + geom_boxplot() + 
  facet_wrap(~Peaks, scales="free_x", strip.position = "bottom") + 
  
  xlab("") + ylab("Peak Length") + ggtitle("Peaks Length Distribution") + ylim(0,3000) + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=cbPalette)
dev.off()


##################################################################################################################
# boxplots of peak length distribution among all_peaks, overlapping_peaks and non-overlapping_peaks 
# expanding upon the above boxplot to compare peak length distr. among individual chromosomes between 
# terra and rloop.
##################################################################################################################

# for All Peaks.
pdf("peakLengths_chr_allPeaks.pdf", width = 8.74, height = 5.03)
combined_allPeaks_only %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Sample)) + 
  geom_boxplot(position = "dodge") + 
  xlab("") + ylab("Peak Length") + ggtitle("Peak Length Distribution") + ylim(0,3000) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "right") + 
  facet_wrap(~Peaks) + 
  scale_fill_manual(values=cbPalette)
dev.off()

# for Non-Overlapping Peaks.
pdf("peakLengths_chr_noOverlap.pdf", width = 8.74, height = 5.03)
combined_noOverlap_only %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Sample)) + 
  geom_boxplot(position = "dodge") + 
  xlab("") + ylab("Peak Length") + ggtitle("Peak Length Distribution") + ylim(0,3000) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "right") + 
  facet_wrap(~Peaks) + 
  scale_fill_manual(values=cbPalette)
dev.off()

# for Overlapping Peaks.
pdf("peakLengths_chr_overlap.pdf", width = 8.74, height = 5.03)
combined_overlap_only %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Sample)) + 
  geom_boxplot(position = "dodge") + 
  xlab("") + ylab("Peak Length") + ggtitle("Peak Length Distribution") + ylim(0,3000) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "right") + 
  facet_wrap(~Peaks) + 
  #guides(fill = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values=cbPalette)
dev.off()


##################################################################################################################
# boxplot of peak length distribution showing R-LOOP overlapping peaks are generally smaller than non-overlapping
##################################################################################################################

pdf("peakLengths_chr_rloopSmallOverlap.pdf", width = 8.74, height = 5.03)
combined_overlap_noOverlap_rloop %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Peaks)) + 
  geom_boxplot() + 
  xlab("") + ylab("Peak Length") + ggtitle("Peak Length Distribution") + ylim(0,3000) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "right") + 
  facet_wrap(~Sample) +
  scale_fill_manual(values=cbPalette)
dev.off()

##################################################################################################################
# boxplot of peak length distribution showing TERRA overlapping peaks are generally smaller than non-overlapping
##################################################################################################################

pdf("peakLengths_chr_terraLargeOverlap.pdf", width = 8.74, height = 5.03)
combined_overlap_noOverlap_terra %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Peaks)) + 
  geom_boxplot() + 
  xlab("") + ylab("Peak Length") + ggtitle("Peak Length Distribution") + ylim(0,3000) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "right") + 
  facet_wrap(~Sample) +
  scale_fill_manual(values=cbPalette)
dev.off()

#########################################################################################################################################
# density of peak length distribution per chromosome.
#########################################################################################################################################

pdf("peakLengths_density.pdf", width = 8.74, height = 5.03)
combined_overlap_only %>% ggplot(.,aes(x=V5, fill=Sample)) + 
  geom_density(alpha=0.5) + 
  xlab("Peak Length (bp)") + ylab("Density") + ggtitle("Overlapping Peak Length Density") + 
  xlim(0,2500) + 
  theme_bw() + 
  scale_fill_manual(values=cbPalette)
dev.off()


#########################################################################################################################################
# peak distribution histogram showing terra overlapping peaks are longer than rloop ones
# i.e. there is a larger amount of long terra peaks, and a larger amount of short rloop peaks
#########################################################################################################################################

pdf("peakLengths_histogram.pdf", width = 8.74, height = 5.03)
combined_overlap_only %>% ggplot(.,aes(x=V5, fill=Sample)) + 
  geom_histogram(binwidth=100, position = "dodge") + 
  scale_x_continuous(limits = c(0,1500), breaks=seq(0,1500, by=100)) + 
  xlab("Peak Length Distribution (bp)") + ylab("Peak Count") + ggtitle("R-loop Overlapping Peaks") + 
  facet_wrap(~Peaks) + 
  theme_bw() + 
  scale_fill_manual(values=cbPalette)
dev.off()


# For summarizing the mean and median peak lengths
peaks %>% group_by(V1) %>% summarise(mean_peak_size=mean(V5),median_peak_size=median(V5)) %>% as.data.frame() %>% arrange(V1)

#########################################################################################################################################
# bar chart of peak distribution per chromosome (nr of peaks per chromosome)
#########################################################################################################################################

pdf("peakCount_chr.pdf", width = 8.74, height = 5.03)
combined_everything %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(y = factor(V1, levels = unique(terra_all$V1)), fill = Sample)) + 
  geom_bar(position=position_dodge()) +
  facet_grid(Peaks ~., scales = "free") + 
  xlab("Peak Count") + ylab("") + ggtitle("Nr. Peaks per Chromosome") + 
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "top", legend.title = element_blank()) + 
  scale_fill_manual(values=cbPalette)
dev.off()




#########################################################################################################################################
# END
#########################################################################################################################################

combined_terra %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = V1)) + 
  geom_boxplot() + 
  xlab("") + ylab("Peak Length") + ggtitle("Peak Length Distribution") + ylim(0,3000) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "none") + 
  facet_grid(Peaks ~factor(V1, levels = unique(terra_all$V1)), scales = "free") 

###

rloop_overlap %>% ggplot(.,aes(V5)) + 
  geom_histogram(aes(fill=..count..),alpha=0.5, binwidth=100) + 
  scale_fill_gradient("Count", low="green", high="red") + 
  geom_density(alpha=0, fill="red") + 
  scale_x_continuous(limits = c(0,1500), breaks=seq(0,1500, by=200)) + 
  xlab("Peak Length Distribution (bp)") + ylab("Peak Count") + ggtitle("TERRA Overlapping Peaks")

terra_overlap %>% ggplot(.,aes(V5)) + 
  geom_histogram(aes(y=..density.., fill=..count..),alpha=0.5, binwidth=100) + 
  scale_fill_gradient("Count", low="green", high="red") + geom_density(alpha=0, fill="red") + 
  scale_x_continuous(limits = c(0,1500), breaks=seq(0,1500, by=200)) + 
  xlab("Peak Length Distribution (bp)") + ylab("Peak Count") + ggtitle("TERRA Overlapping Peaks")

rloop_overlap %>% ggplot(.,aes(V5)) + 
  geom_histogram(aes(fill=..count..),alpha=0.5, binwidth=100) + 
  scale_fill_gradient("Count", low="green", high="red") + 
  geom_density(alpha=0, fill="red") + 
  scale_x_continuous(limits = c(0,1500), breaks=seq(0,1500, by=200)) + 
  xlab("Peak Length Distribution (bp)") + ylab("Peak Count") + ggtitle("R-loop Overlapping Peaks")


###

combined_everything %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(y = factor(V1, levels = unique(terra_all$V1)), fill = Sample)) + 
  geom_bar(position=position_dodge()) +
  facet_wrap(~Peaks, scales = "free") + 
  xlab("Peak Count") + ylab("") + ggtitle("Nr. Peaks per Chromosome") + 
  scale_fill_manual(values=cbPalette)

combined_everything %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(y = factor(V1, levels = unique(terra_all$V1)), fill = Sample)) + 
  geom_bar(position=position_dodge()) +
  facet_grid(Peaks ~factor(V1, levels = unique(terra_all$V1)), scales = "free") + 
  xlab("Peak Count") + ylab("") + ggtitle("Nr. Peaks per Chromosome") + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "top") + 
  scale_fill_manual(values=cbPalette)

