#########################################################################################################################################
# IMPORTANT NOTE: files being read below only contain peaks with at least 2 repeats. all peaks with only one repeat were excluded, as they
#                 were too many and shifted the mean in the figures significantly such that important data were not visible. Also, we did
#                 not consider peaks with only one repeat to be of importance/significance in R-loop formation, as more consecutive
#                 repeats are needed for that.
#########################################################################################################################################

# set working directory.
setwd("/proj/nb_storage/private/terra_rloop_project/results/figures/repeat_count_distribution/")

# load packages.
library("data.table")
library("ggplot2")
library("dplyr")

# read files.
terra_all <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/repeat_count_distribution/terra_peaks_repCount.txt", header = F,sep="\t")
terra_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/repeat_count_distribution/terra_overlap_repCount.txt", header = F,sep="\t")
terra_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/repeat_count_distribution/terra_noOverlap_repCount.txt", header = F,sep="\t")
rloop_all <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/repeat_count_distribution/rloop_peaks_repCount.txt", header = F,sep="\t")
rloop_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/repeat_count_distribution/rloop_overlap_repCount.txt", header = F,sep="\t")
rloop_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/repeat_count_distribution/rloop_noOverlap_repCount.txt", header = F,sep="\t")

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

############################################################################################################################
# boxplot of repeat count distribution per peak file
############################################################################################################################

pdf("Fig2_repeatCountDistribution_boxplot.pdf", width = 4, height = 4.5)

combined_all %>% filter(Sample == "Terra" & Peak_Group != "All") %>% 
  ggplot(.,aes(x=Peak_Group,y=V2,fill=Peak_Group)) + geom_boxplot() + 
  xlab("") + ylab("Nr. of Repeats") + 
  labs(caption = "Only peaks with at least 1 repeat were included") + 
  
  facet_wrap(~Sample, scales = "free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), legend.title = element_blank(), legend.position = "top", plot.caption = element_text(hjust = 0.5, vjust = 141)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_y_continuous(limits = c(0,120), breaks=seq(0,120, by=20))

dev.off()

# to lock aspect ratio while resizing during savind plots.
# theme(aspect.ratio = 4.5/8)
#############################################################################################################################
# density plot: repeat distribution among peak groups
#############################################################################################################################


# pdf("repeatCountDistribution_densityplot.pdf", width = 8, height = 4.5)
# 
# combined_all %>% filter(Peak_Group != "All") %>% ggplot(.,aes(x=V2,fill=Peak_Group)) + geom_density(alpha=0.5) + 
#   xlab("Repeats") + ylab("Density") + ggtitle("Repeat Count Distribution") + 
#   ylim(0,0.2) + 
#   facet_wrap(~Sample) + 
#   theme_bw() + 
#   scale_fill_manual(values=cbPalette)
# 
# dev.off()
