# set working directory.
setwd("/proj/nb_storage/private/terra_rloop_project/results/figures/consecutive_repeats_distribution/")

# load packages.
library("data.table")
library("ggplot2")
library("dplyr")
#library(tibble) # to add columns in between existing ones.


# read files.
#max_fields <- max(count.fields("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/consecutive_repeats_distribution/terra_peaks_consRepCount.txt", sep="\t"))
#, col.names = paste0("V",seq_len(max_fields)
terra_all<- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/consecutive_repeats_distribution/terra_peaks_consRepCount.txt", header = T, sep="\t", fill = T)
terra_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/consecutive_repeats_distribution/terra_overlap_consRepCount.txt", header = T, sep="\t", fill = T)
terra_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/consecutive_repeats_distribution/terra_noOverlap_consRepCount.txt", header = T, sep="\t", fill = T)

rloop_all <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/consecutive_repeats_distribution/rloop_peaks_consRepCount.txt", header = T, sep="\t", fill = T)
rloop_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/consecutive_repeats_distribution/rloop_overlap_consRepCount.txt", header = T, sep="\t", fill = T)
rloop_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/consecutive_repeats_distribution/rloop_noOverlap_consRepCount.txt", header = T, sep="\t", fill = T)

# to combine/bind data frames, they must have the same nr of columns. Since this is not the case, select only columns 1 to 3 form each.
#terra_all <- select(terra_all,Peak_ID,Cons_Repeat_Groups,Max_Cons_Repeats)
#terra_overlap <- select(terra_overlap,Peak_ID,Cons_Repeat_Groups,Max_Cons_Repeats)
#terra_noOverlap <- select(terra_noOverlap,Peak_ID,Cons_Repeat_Groups,Max_Cons_Repeats)

#rloop_all <- select(rloop_all,Peak_ID,Cons_Repeat_Groups,Max_Cons_Repeats)
#rloop_overlap <- select(rloop_overlap,Peak_ID,Cons_Repeat_Groups,Max_Cons_Repeats)
#rloop_noOverlap <- select(rloop_noOverlap,Peak_ID,Cons_Repeat_Groups,Max_Cons_Repeats)


# add new columns to existing data frames to be able to create boxplots.
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

# filter to remove rows with ZERO groups of consecutive (4+) repeats. Cons_Repeat_Groups != 0
# filter to remove rows with only 1 repeat. Max_Cons_Repeats > 1
# terra_all <- filter(terra_all, Cons_Repeat_Groups != 0)
# terra_overlap <- filter(terra_overlap, Cons_Repeat_Groups != 0)
# terra_noOverlap <- filter(terra_noOverlap, Cons_Repeat_Groups != 0)
# 
# rloop_all <- filter(rloop_all, Cons_Repeat_Groups != 0)
# rloop_overlap <- filter(rloop_overlap, Cons_Repeat_Groups != 0)
# rloop_noOverlap <- filter(rloop_noOverlap, Cons_Repeat_Groups != 0)

# combine all data frames.
combined_all <- rbind(terra_all,terra_overlap,terra_noOverlap,rloop_all,rloop_overlap,rloop_noOverlap)
combined_all_noOverlap <- rbind(terra_all,terra_noOverlap,rloop_all,rloop_noOverlap)
combined_overlap <- rbind(terra_overlap,rloop_overlap)
combined_terra_overlap_noOverlap <- rbind(terra_noOverlap,terra_overlap)
combined_rloop_overlap_noOverlap <- rbind(rloop_noOverlap,rloop_overlap)
combined_overlap_noOverlap <- rbind(terra_overlap,terra_noOverlap,rloop_overlap,rloop_noOverlap)
combined_terra <- rbind(terra_all,terra_overlap,terra_noOverlap)

#########################################################################################################################################
# define color-blind palettes.
#########################################################################################################################################

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#########################################################################################################################################
# boxplot of max consecutive repeat count distribution per peak file
#########################################################################################################################################



combined_terra_overlap_noOverlap %>% filter(Max_Cons_Repeats > 1) %>%    # filter to remove rows with only 1 repeat. Max_Cons_Repeats > 1
  ggplot(.,aes(x=Sample,y=Max_Cons_Repeats,fill=Peaks)) + geom_boxplot() + 
  xlab("Peaks") + ylab("Maximum Consecutive Repeats Count") + ggtitle("Consecutive (4+) Repeat Distribution") + 
  facet_wrap(~Sample) + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=cbPalette)

# final
pdf("maxConsRepeatCount.pdf", width = 8.74, height = 5.03)

combined_all %>% filter(Max_Cons_Repeats > 1) %>%  # filter to remove peaks with only 1 repeat. Max_Cons_Repeats > 1
  ggplot(.,aes(x=Sample,y=Max_Cons_Repeats,fill=Peaks)) + 
  geom_boxplot() + 
  xlab("Peaks") + ylab("Maximum Nr. of Consecutive Repeats") + 
  ggtitle("Consecutive Repeat Distribution", 
          subtitle = "Only considers peaks with at least 2 repeats. Consecutive = at least 4 in a row") + 
  facet_wrap(~Sample, scales = "free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),              # remove text from x-axis ticks
        legend.title = element_blank()) +           # remove legend (fill) title
  scale_fill_manual(values=cbPalette) +             # colour by cbPalette
  scale_y_continuous(limits = c(0,90), breaks=seq(0,90, by=10))


dev.off()

pdf("ConsRepeatGroupCount_overlapVSnoOverlap.pdf", width = 8.74, height = 5.03)

combined_all %>% filter(Max_Cons_Repeats > 1) %>%  # filter to remove peaks with only 1 repeat. Max_Cons_Repeats > 1
  filter(Peaks != "All") %>%         # filter to remove 'All' peaks and only keep 'Overlapping' and 'Non-Overlapping'
  
  ggplot(.,aes(x=Sample,y=Cons_Repeat_Groups,fill=Peaks)) + 
  geom_boxplot() + 
  xlab("Peaks") + ylab("Nr. of Groups of Consecutive Repeats ") + 
  ggtitle("Consecutive Repeat Groups Distribution", 
          subtitle = "Only considers peaks with at least 2 repeats. Consecutive = at least 4 in a row") + 
  facet_wrap(~Sample, scales = "free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),              # remove text from x-axis ticks
        legend.title = element_blank()) +           # remove legend (fill) title
  scale_fill_manual(values=cbPalette) +             # colour by cbPalette
  scale_y_continuous(limits = c(0,10), breaks=seq(0,10, by=1))


dev.off()

pdf("ConsRepeatGroupCount_allVSnoOverlap.pdf", width = 8.74, height = 5.03)

combined_all %>% filter(Max_Cons_Repeats > 1) %>%  # filter to remove peaks with only 1 repeat. Max_Cons_Repeats > 1
  filter(Peaks != "Overlapping") %>%         # filter to remove 'All' peaks and only keep 'Overlapping' and 'Non-Overlapping'
  
  ggplot(.,aes(x=Sample,y=Cons_Repeat_Groups,fill=Peaks)) + 
  geom_boxplot() + 
  xlab("Peaks") + ylab("Nr. of Groups of Consecutive Repeats ") + 
  ggtitle("Consecutive Repeat Groups Distribution", 
          subtitle = "Only considers peaks with at least 2 repeats. Consecutive = at least 4 in a row") + 
  facet_wrap(~Sample, scales = "free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),              # remove text from x-axis ticks
        legend.title = element_blank()) +           # remove legend (fill) title
  scale_fill_manual(values=cbPalette) +             # colour by cbPalette
  scale_y_continuous(limits = c(0,10), breaks=seq(0,10, by=1))


dev.off()


############################################### START TEST #############################################################################

combined_terra_overlap_noOverlap %>% ggplot(.,aes(x=Sample,y=Cons_Repeat_Groups,fill=Peaks)) + geom_boxplot() + 
  xlab("Peaks") + ylab("Nr. Groups of 4+ Consecutive Repeats") + ggtitle("Consecutive (4+) Repeat Distribution")

combined_all_noOverlap %>% ggplot(.,aes(x=Sample,y=Cons_Repeat_Groups,fill=Peaks)) + geom_boxplot() + 
  xlab("Peaks") + ylab("Repeat Count") + ggtitle("Repeat Count Distribution")

combined_overlap_noOverlap %>% ggplot(.,aes(x=Sample,y=Cons_Repeat_Groups,fill=Peaks)) + geom_boxplot() + 
  xlab("Peaks") + ylab("Repeat Count") + ggtitle("Repeat Count Distribution")
#+ theme(axis.text.x = element_text(angle = 60, hjust = 1))

############################################### End TEST #############################################################################


#########################################################################################################################################
# peak distribution density plot
#########################################################################################################################################

combined_overlap_noOverlap %>% filter(Cons_Repeat_Groups > 0) %>% 
  ggplot(.,aes(x=Max_Cons_Repeats,fill=Peaks)) + geom_density(alpha=0.5) + 
  xlab("Consecutive Repeat Count") + ylab("Density") + ggtitle("Maximum Consecutive Repeats Count Distribution") +
  facet_wrap(~Sample)

combined_overlap_noOverlap %>% filter(Cons_Repeat_Groups > 0) %>% 
  ggplot(.,aes(x=Max_Cons_Repeats,fill=Peaks)) + geom_bar(alpha=0.5, color="grey") + 
  xlab("Consecutive Repeat Count") + ylab("Nr. Peaks") + ggtitle("Maximum Consecutive Repeats Count Distribution") +
  facet_wrap(~Sample)


combined_overlap_noOverlap %>% filter(Cons_Repeat_Groups > 0) %>% 
  ggplot(.,aes(x=Max_Cons_Repeats,fill=Peaks)) + geom_density(alpha=0.5) + 
  xlab("Consecutive Repeat Count") + ylab("Density") + ggtitle("Maximum Consecutive Repeats Count Distribution") +
  facet_wrap(~Sample) + 
  xlim(0,25)

combined_overlap_noOverlap %>% filter(Cons_Repeat_Groups > 0) %>% 
  ggplot(.,aes(x=Cons_Repeat_Groups,fill=Peaks)) + geom_density(alpha=0.5) + 
  xlab("Nr. Groups of 4+ Consecutive Repeats") + ylab("Density") + ggtitle("Groups of Consecutive (4+) Repeats") +
  facet_wrap(~Sample) + ylim(0,1)
  #xlim(0,5)

# draft density plot
combined_terra %>% filter(Peaks != "All")%>% filter(Cons_Repeat_Groups > 0) %>% 
  ggplot(.,aes(x=Cons_Repeat_Groups,fill=Peaks)) + geom_density(alpha=0.5) + 
  xlab("Nr. Groups of 4+ Consecutive Repeats") + ylab("Density") + ggtitle("Groups of Consecutive (4+) Repeats") +
  facet_wrap(~Peaks, scales = "free_y") + 
  scale_x_continuous(limits = c(0,10), breaks=seq(0,10, by=2))
  
#########################################################################################################################################
# peak distribution bar plot
#########################################################################################################################################

# bar plot of TERRA overlapping vs non-overlapping peaks showing how many peaks contain >=1 group of consecutive repeats.

pdf("ConsRepeatGroupCount_overlapVSnoOverlap_barplot.pdf", width = 8.74, height = 5.03)

combined_terra %>% filter(Peaks != "All")%>% filter(Cons_Repeat_Groups >= 1) %>% 
  ggplot(.,aes(x=Cons_Repeat_Groups, fill=Peaks)) + geom_bar() + 
  xlab("Nr. of Groups") + ylab("Nr. of Peaks") + 
  ggtitle("Terra Peaks", subtitle =  "Groups of Consecutive (4+) Repeats") +
  facet_wrap(~Peaks) + 
  theme_bw() + 
  theme(strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,135), breaks=seq(0,135, by=10)) +
  scale_x_continuous(limits = c(0,11), breaks=seq(0,11, by=2)) + 
  scale_fill_manual(name = "Terra Peaks", labels = c("Non-Overlapping", "Overlapping"), values=cbPalette)

dev.off()

# bar plot of TERRA overlapping vs non-overlapping peaks showing how many peaks contain >=2 group of consecutive repeats.

pdf("ConsRepeatGroupCount_overlapVSnoOverlap_min2groups_barplot.pdf", width = 8.74, height = 5.03)

combined_terra %>% filter(Peaks != "All")%>% filter(Cons_Repeat_Groups >= 2) %>% 
  ggplot(.,aes(x=Cons_Repeat_Groups,fill=Peaks)) + geom_bar() + 
  xlab("Nr. of Groups") + ylab("Nr. of Peaks") + 
  ggtitle("Terra Peaks", subtitle =  "Groups of Consecutive (4+) Repeats") +
  facet_wrap(~Peaks) + 
  theme_bw() + 
  theme(strip.text.x = element_blank()) +
  #theme(strip.background = element_blank(), strip.text.x = element_blank()  ) +
  scale_y_continuous(limits = c(0,24), breaks=seq(0,24, by=2)) +
  scale_x_continuous(limits = c(1,11), breaks=seq(0,11, by=2)) + 
  #scale_fill_discrete(name = "Terra Peaks", labels = c("Non-Overlapping", "Overlapping")) + 
  scale_fill_manual(name = "Terra Peaks", labels = c("Non-Overlapping", "Overlapping"), values=cbPalette)

dev.off()

# draft
combined_terra %>% filter(Peaks != "All")%>% filter(Cons_Repeat_Groups >= 1) %>% 
  ggplot(.,aes(x=Cons_Repeat_Groups,fill=Peaks)) + geom_bar() + 
  xlab("Nr. of Groups") + ylab("Nr. of Peaks") + ggtitle("Groups of Consecutive (4+) Repeats") +
  facet_wrap(~Peaks, scales = "free") + 
  theme_bw() 


#########################################################################################################################################

# END

#########################################################################################################################################

