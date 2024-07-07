#########################################################################################################################################
# IMPORTANT NOTE: files being read below only contain peaks with at least 1 repeats. all peaks with only one repeat were excluded, as they
#                 were too many and shifted the mean in the figures significantly such that important data were not visible. Also, we did
#                 not consider peaks with only one repeat to be of importance/significance in R-loop formation, as more consecutive
#                 repeats are needed for that.
#########################################################################################################################################

# set working directory.
setwd("/cluster/khiom/sshehata001/proj/terra_rloop/results/figures/Rplots")

# load packages.
library("data.table")
library("ggplot2")
library("dplyr")

# read files.
peak_statistics<-fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/repeat_stats/peak_statistics.tsv",header = T,sep="\t")

terra_all <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/repeat_counts/terra_peaks_motif_counts.txt", header = F,sep="\t")
terra_no_intersect <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/repeat_counts/terra_no_intersect_motif_counts.txt", header = F,sep="\t")
terra_intersect <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/repeat_counts/terra_intersect_motif_counts.txt", header = F,sep="\t")
rloop_all <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/repeat_counts/rloop_peaks_motif_counts.txt", header = F,sep="\t")
rloop_no_intersect <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/repeat_counts/rloop_no_intersect_motif_counts.txt", header = F,sep="\t")
rloop_intersect <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/repeat_counts/rloop_intersect_motif_counts.txt", header = F,sep="\t")

# add new column to existing data frames to be able to create boxplots.
terra_all$Sample="TERRA"
terra_all$Peak_Group="All"

terra_intersect$Sample="TERRA"
terra_intersect$Peak_Group="Intersecting"

terra_no_intersect$Sample="TERRA"
terra_no_intersect$Peak_Group="Non-Intersecting"

rloop_all$Sample="R-loop"
rloop_all$Peak_Group="All"

rloop_intersect$Sample="R-loop"
rloop_intersect$Peak_Group="Intersecting"

rloop_no_intersect$Sample="R-loop"
rloop_no_intersect$Peak_Group="Non-Intersecting"

# combine all data frames.
combined_all <- rbind(terra_all,terra_no_intersect,terra_intersect,rloop_all,rloop_no_intersect,rloop_intersect)

# extract peak numbers with at least 1 repeat for non-intersecting and intersecting terra peaks
terra_no_intersect_1r <- as.numeric((peak_statistics %>% filter(Sample == "TERRA" & Peak_Group == "Non-Intersecting")) %>% select(Min_1R))
terra_intersect_1r <- as.numeric((peak_statistics %>% filter(Sample == "TERRA" & Peak_Group == "Intersecting")) %>% select(Min_1R))

# define color-blind palettes.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

############################################################################################################################
# boxplot of repeat count distribution per peak file
############################################################################################################################

# create plot
pdf("fig2_repeat_counts_violinplot.pdf", width = 4, height = 4.5)

combined_all %>% filter(Sample == "TERRA" & Peak_Group != "All") %>% 
  ggplot(.,aes(x=factor(Peak_Group, levels = unique(Peak_Group)), y=V2, fill=factor(Peak_Group, levels = unique(Peak_Group)))) + 
  geom_violin() + 
  xlab("") + ylab("Nr. of repeats per peak") + 
  # labs(caption = "Only peaks with repeats >= 1 were included") + 
  
  facet_wrap(~Sample, scales = "free_x") + 
  theme_bw() + 
  theme(legend.title = element_blank(), 
        legend.position = "top", 
        # plot.caption = element_text(hjust = 0.5, vjust = 70), 
        text = element_text(size = 22)) + 
  scale_fill_manual(values=cbPalette) + 
  scale_y_continuous(limits = c(0,120), breaks=seq(0,120, by=20)) +
  scale_x_discrete(labels=c("Non-Intersecting" = terra_no_intersect_1r, "Intersecting" = terra_intersect_1r)) +
  stat_compare_means(label = "p.format", method = "t.test")
  # stat_compare_means(label = "p.signif", method = "t.test")


dev.off()

#  geom_dotplot(binaxis = 'y', stackdir = "center", binwidth = 2, stackratio = 0.1)


### Extra calculating median.
combined_all %>% filter(Sample == "TERRA" & Peak_Group == "Intersecting") %>% pull(V2) %>% median()# [1] 11

combined_all %>% filter(Sample == "TERRA" & Peak_Group == "Non-Intersecting") %>% pull(V2) %>% median()# [1] 1


combined_all %>%
  filter(V2 > 1) %>% 
  group_by(Sample, Peak_Group) %>%
  summarise(nrs = paste(V2, collapse = ",")) %>% 
  mutate(
    nrs = str_split(nrs, ","),
    median = median(as.integer(unlist(nrs))))

### Extra testing distribution.
combined_all %>%
  filter(Sample == "TERRA" & Peak_Group != "All" & V2 > 1 & V2 < 30) %>%
  ggplot((aes(V2))) +
  geom_bar()

