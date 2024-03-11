#########################################################################################################################################
# IMPORTANT NOTE: files being read below only contain peaks with at least 2 repeats. all peaks with only one repeat were excluded, as they
#                 were too many and shifted the mean in the figures significantly such that important data were not visible. Also, we did
#                 not consider peaks with only one repeat to be of importance/significance in R-loop formation, as more consecutive
#                 repeats are needed for that.
#########################################################################################################################################

# set working directory.
setwd("/proj/nb_storage/private/terra_rloop_project/results/figures/repeat_count_distribution/")

library("data.table")
library("ggplot2")
# tibble needed for add_column() function.
#library("tibble", lib.loc="/domus/h1/shehata/R/x86_64-redhat-linux-gnu-library/3.6")

# read files.
terra_all <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/repeat_count_distribution/terra_peaks_repCount.txt", header = F,sep="\t")
terra_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/repeat_count_distribution/terra_overlap_repCount.txt", header = F,sep="\t")
terra_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/repeat_count_distribution/terra_noOverlap_repCount.txt", header = F,sep="\t")
rloop_all <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/repeat_count_distribution/rloop_peaks_repCount.txt", header = F,sep="\t")
rloop_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/repeat_count_distribution/rloop_overlap_repCount.txt", header = F,sep="\t")
rloop_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/repeat_count_distribution/rloop_noOverlap_repCount.txt", header = F,sep="\t")

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
combined_all_noOverlap <- rbind(terra_all,terra_noOverlap,rloop_all,rloop_noOverlap)
combined_overlap <- rbind(terra_overlap,rloop_overlap)
combined_terra <- rbind(terra_noOverlap,terra_overlap)
combined_rloop <- rbind(rloop_noOverlap,rloop_overlap)
combined_overlap_noOverlap <- rbind(terra_overlap,terra_noOverlap,rloop_overlap,rloop_noOverlap)

#########################################################################################################################################
# define color-blind palettes.
#########################################################################################################################################

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#########################################################################################################################################
# boxplot of repeat count distribution per peak file
#########################################################################################################################################

# + scale_y_continuous(limits = c(0,30), breaks=seq(0,30, by=5))

pdf("repeatCountDistribution_boxplot.pdf", width = 8.74, height = 5.03)

combined_all %>% ggplot(.,aes(x=Sample,y=V2,fill=Peaks)) + geom_boxplot() + 
  xlab("") + ylab("Repeat Count") + ggtitle("Repeat Count Distribution") + 
  facet_wrap(~Sample, scales = "free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=cbPalette)

dev.off()

combined_all_noOverlap %>% ggplot(.,aes(x=Sample,y=V2,fill=Peaks)) + geom_boxplot() + 
  xlab("Peaks") + ylab("Repeat Count") + ggtitle("Repeat Count Distribution")

combined_overlap_noOverlap %>% ggplot(.,aes(x=Sample,y=V2,fill=Peaks)) + geom_boxplot() + 
  xlab("Peaks") + ylab("Repeat Count") + ggtitle("Repeat Count Distribution")
#+ theme(axis.text.x = element_text(angle = 60, hjust = 1))


#########################################################################################################################################
# peak distribution density plot
#########################################################################################################################################

pdf("", width = 8.74, height = 5.03)

combined_overlap_noOverlap %>% ggplot(.,aes(x=V2,fill=Peaks)) + geom_density(alpha=0.5) + 
  xlab("Repeat Count") + ylab("Density") + ggtitle("Repeat Count Distribution") + 
  ylim(0,0.2) + 
  facet_wrap(~Sample) + 
  theme_bw() + 
  scale_fill_manual(values=cbPalette)

dev.off()

combined_terra %>% ggplot(.,aes(x=V2,fill=Peaks)) + geom_density(alpha=0.5) + 
  xlab("Repeat Count") + ylab("Density") + ggtitle("Terra Repeat Count Distribution")  + 
  xlim(0,40) + 
  theme_minimal() 

combined_rloop %>% ggplot(.,aes(x=V2,fill=Peaks)) + geom_density(alpha=0.5) + 
  xlab("Repeat Count") + ylab("Density") + ggtitle("R-Loop Repeat Count Distribution")  + 
  ylim(0,0.2) + xlim(0,40) + 
  theme_minimal() 

combined_all %>% ggplot(.,aes(x=V2,fill=Peaks)) + geom_density(alpha=0.5) + 
  xlab("Repeat Count") + ylab("Density") + ggtitle("Repeat Count Distribution") + 
  ylim(0,0.2) + xlim(0,40) + 
  facet_wrap(~Sample)

#########################################################################################################################################
# END
#########################################################################################################################################


terra_overlap %>% ggplot(.,aes(x=V2,fill=Peaks)) + geom_density() + xlab("Repeat Count") + ggtitle("Repeat Count Distribution")


#########################################################################################################################################

# peak distribution histogram

#########################################################################################################################################

# + geom_histogram(binwidth=10)+ scale_x_continuous(limits = c(0,1500), breaks=seq(0,1500, by=200))
combined_terra %>% ggplot(.,aes(x=V2,fill=Peaks)) + 
  geom_histogram(binwidth = 1) + 
  xlab("Repeat Count") + ylab("Peak Count") + ggtitle("Repeat Count Distribution")

combined_all %>% ggplot(.,aes(x=V2,fill=Peaks)) + 
  geom_histogram(binwidth = 1,alpha=0.5) + 
  xlab("Repeat Count") + ylab("Peak Count") + ggtitle("Repeat Count Distribution")

