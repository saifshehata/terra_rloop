# set working directory.
setwd("/cluster/khiom/sshehata001/proj/terra_rloop/results/figures/Rplots")

# load packages.
library("data.table")
library("ggplot2")
library("dplyr")


# read files.
peak_statistics<-fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/repeat_stats/peak_statistics.tsv",header = T,sep="\t")

terra_all <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/tandem_repeats/terra_peaks_tandem_motif_counts.txt", header = T,sep="\t")
terra_no_intersect <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/tandem_repeats/terra_no_intersect_tandem_motif_counts.txt", header = T,sep="\t")
terra_intersect <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/tandem_repeats/terra_intersect_tandem_motif_counts.txt", header = T,sep="\t")
rloop_all <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/tandem_repeats/rloop_peaks_tandem_motif_counts.txt", header = T,sep="\t")
rloop_no_intersect <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/tandem_repeats/rloop_no_intersect_tandem_motif_counts.txt", header = T,sep="\t")
rloop_intersect <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/tandem_repeats/rloop_intersect_tandem_motif_counts.txt", header = T,sep="\t")

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
terra_no_intersect_4r <- as.numeric((peak_statistics %>% filter(Sample == "TERRA" & Peak_Group == "Non-Intersecting")) %>% select(Min_4tR))
terra_intersect_4r <- as.numeric((peak_statistics %>% filter(Sample == "TERRA" & Peak_Group == "Intersecting")) %>% select(Min_4tR))

# define color-blind palettes.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##################################################################################################################
# violinplot: maximum tandem repeat count distribution per peak group
##################################################################################################################

pdf("Fig2_longest_tandem_repeats_violinplot.pdf", width = 4, height = 4.5)

combined_all %>% filter(Max_Tandem_Repeats >= 4) %>%  # filter to remove peaks with only 1 repeat. Max_Tandem_Repeats > 1
  filter(Sample == "TERRA" & Peak_Group != "All") %>%  
  ggplot(.,aes(x=factor(Peak_Group, levels = unique(Peak_Group)),y=Max_Tandem_Repeats,fill=factor(Peak_Group, levels = unique(Peak_Group)))) + 
  geom_violin() + 
  xlab("") + 
  ylab("Longest stretch of tandem repeats per peak") +
  # ylab("Tandem repeats per peak") + 
  # labs(caption = "Only peaks with repeats >= 4 tandem were included") + 
  facet_wrap(~Sample, scales = "free_x") + 
  theme_bw() + 
  theme(legend.title = element_blank(),
        legend.position = "top",
        # plot.caption = element_text(hjust = 0.5, vjust = 70),
        text = element_text(size = 22) ) +           
  scale_fill_manual(values=cbPalette) +             # colour by cbPalette
  scale_y_continuous(limits = c(0,90), breaks=seq(0,90, by=20)) +
  scale_x_discrete(labels=c("Non-Intersecting" = terra_no_intersect_4r, "Intersecting" = terra_intersect_4r))

dev.off()


### Extra calculating median.
combined_all %>% filter(Sample == "TERRA" & Peak_Group == "Intersecting") %>% pull(Max_Tandem_Repeats) %>% median()# [1] 8

combined_all %>% filter(Sample == "TERRA" & Peak_Group == "Non-Intersecting") %>% pull(Max_Tandem_Repeats) %>% median()# [1] 1


# combined_all %>%
#   filter(Max_Tandem_Repeats > 1) %>%
#   group_by(Sample, Peak_Group) %>%
#   summarise(nrs = paste(Max_Tandem_Repeats, collapse = ",")) %>% 
#   mutate(
#     nrs = str_split(nrs, ","),
#     median = median(as.integer(unlist(nrs))))

### Extra testing distribution.
combined_all %>%
  filter(Sample == "TERRA" & Peak_Group != "All" & Max_Tandem_Repeats > 1 & Max_Tandem_Repeats < 30) %>%
  ggplot((aes(Max_Tandem_Repeats))) +
  geom_bar()
