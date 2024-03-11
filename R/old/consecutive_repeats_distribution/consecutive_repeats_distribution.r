# set working directory.
setwd("/proj/nb_storage/private/terra_rloop_project/results/figures/consecutive_repeats_distribution/")

# load packages.
library("data.table")
library("ggplot2")
library("dplyr")


# read files.
terra_all<- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/consecutive_repeats_distribution/terra_peaks_consRepCount.txt", header = T, sep="\t", fill = T)
terra_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/consecutive_repeats_distribution/terra_overlap_consRepCount.txt", header = T, sep="\t", fill = T)
terra_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/consecutive_repeats_distribution/terra_noOverlap_consRepCount.txt", header = T, sep="\t", fill = T)

rloop_all <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/consecutive_repeats_distribution/rloop_peaks_consRepCount.txt", header = T, sep="\t", fill = T)
rloop_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/consecutive_repeats_distribution/rloop_overlap_consRepCount.txt", header = T, sep="\t", fill = T)
rloop_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/consecutive_repeats_distribution/rloop_noOverlap_consRepCount.txt", header = T, sep="\t", fill = T)

# add new columns to existing data frames to be able to create boxplots.
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
# boxplot: maximum consecutive repeat count distribution per peak type
##################################################################################################################

pdf("Fig3_maxConsRepeatCount_boxplot.pdf", width = 4, height = 4.5)

combined_all %>% filter(Max_Cons_Repeats >= 4) %>%  # filter to remove peaks with only 1 repeat. Max_Cons_Repeats > 1
  filter(Sample == "Terra" & Peak_Group != "All") %>%  
  ggplot(.,aes(x=Peak_Group,y=Max_Cons_Repeats,fill=Peak_Group)) + 
  geom_boxplot() + 
  xlab("") + ylab("Longest Stretch of Consecutive Repeats") + 
  labs(caption = "Only peaks with at least 4 consecutive repeats were included") + 
  facet_wrap(~Sample, scales = "free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),              # remove text from x-axis ticks
        legend.title = element_blank(),             # remove legend (fill) title
        legend.position = "top",
        plot.caption = element_text(hjust = 0.5, vjust = 141) ) +           
  scale_fill_manual(values=cbPalette) +             # colour by cbPalette
  scale_y_continuous(limits = c(0,90), breaks=seq(0,90, by=20))

dev.off()

#################################################################################################################
# bar plot: how many peaks contain several groups of consecutive repeats. 
#################################################################################################################


# TERRA overlapping vs non-overlapping peaks containing >=1 group of consecutive repeats.
pdf("Fig3_ConsRepeatGroupCount_overlapVSnoOverlap_barplot.pdf", width = 8, height = 4.5)

combined_all %>% filter(Sample == "Terra") %>% filter(Peak_Group != "All")%>% filter(Cons_Repeat_Groups >= 1) %>% 
  # ggplot(.,aes(x=Cons_Repeat_Groups)) + geom_bar() +
  ggplot(.,aes(x=Cons_Repeat_Groups, fill=Peak_Group)) + 
  geom_bar(position = position_dodge2(preserve = "single")) +
  xlab("Nr. of Repeat Groups") + ylab("Nr. of Peaks") + 
  # ggtitle("Terra Peaks", subtitle =  "Groups of Consecutive (4+) Repeats") +
  labs(subtitle = "Terra Peaks") + 
  # facet_wrap(~Peak_Group) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,135), breaks=seq(0,135, by=10)) +
  scale_x_continuous(limits = c(0,10.5), breaks=seq(0,10, by=1)) + 
  
  theme(plot.subtitle = element_text(vjust = -14), legend.position = "top", legend.title = element_blank()) +
  scale_fill_manual(name = "Terra Peaks", labels = c("Non-Overlapping", "Overlapping"), values=cbPalette)

dev.off()

# # same but split into 2 facets.
# combined_all %>% filter(Sample == "Terra") %>% filter(Peak_Group != "All")%>% filter(Cons_Repeat_Groups >= 1) %>% 
#   # ggplot(.,aes(x=Cons_Repeat_Groups)) + geom_bar() +
#   ggplot(.,aes(x=Cons_Repeat_Groups, fill=Peak_Group)) + geom_bar() +
#   xlab("Nr. of Repeat Groups") + ylab("Nr. of Peaks") + 
#   # ggtitle("Terra Peaks", subtitle =  "Groups of Consecutive (4+) Repeats") +
#   labs(subtitle = "Terra Peaks") + 
#   facet_wrap(~Peak_Group) + 
#   theme_bw() + 
#   scale_y_continuous(limits = c(0,135), breaks=seq(0,135, by=10)) +
#   scale_x_continuous(limits = c(0,11), breaks=seq(0,11, by=2)) + 
#   
#   theme(strip.text.x = element_blank(), plot.subtitle = element_text(vjust = -14), legend.position = "top", legend.title = element_blank()) +
#   scale_fill_manual(name = "Terra Peaks", labels = c("Non-Overlapping", "Overlapping"), values=cbPalette)




# # TERRA overlapping vs non-overlapping peaks containing >=2 group2 of consecutive repeats.
# pdf("ConsRepeatGroupCount_overlapVSnoOverlap_min2groups_barplot.pdf", width = 8, height = 4.5)
# 
# combined_all %>% filter(Sample == "Terra") %>% filter(Peak_Group != "All")%>% filter(Cons_Repeat_Groups >= 1) %>%
#   # ggplot(.,aes(x=Cons_Repeat_Groups)) + geom_bar() +
#   ggplot(.,aes(x=Cons_Repeat_Groups, fill=Peak_Group)) + geom_bar() +
#   xlab("Nr. of Repeat Groups") + ylab("Nr. of Peaks") +
#   # ggtitle("Terra Peaks", subtitle =  "Groups of Consecutive (4+) Repeats") +
#   labs(subtitle = "Terra Peaks") +
#   facet_wrap(~Peak_Group) +
#   theme_bw() +
#   scale_y_continuous(limits = c(0,24), breaks=seq(0,24, by=2)) +
#   scale_x_continuous(limits = c(1,11), breaks=seq(0,11, by=2)) +
# 
#   theme(strip.text.x = element_blank(), plot.subtitle = element_text(vjust = -14), legend.position = "top", legend.title = element_blank()) +
#   scale_fill_manual(name = "Terra Peaks", labels = c("Non-Overlapping", "Overlapping"), values=cbPalette)
# #scale_fill_discrete(name = "Terra Peaks", labels = c("Non-Overlapping", "Overlapping"))
# 
# dev.off()





##################################################################################################################
# boxplot: distribution of number/count of consecutive repeat groups per peak type
##################################################################################################################

# overlapping vs non-overlapping peaks
pdf("Fig3_ConsRepeatGroupCount_overlapVSnoOverlap_boxplot.pdf", width = 4, height = 4.5)

combined_all %>% filter(Max_Cons_Repeats > 1) %>%  # filter to remove peaks with only 1 repeat. Max_Cons_Repeats > 1
  filter(Sample == "Terra" & Peak_Group != "All") %>%  
  # filter to keep only terra peaks, and remove 'All' peaks and only keep 'Overlapping' and 'Non-Overlapping'
   
  ggplot(.,aes(x=Peak_Group,y=Cons_Repeat_Groups,fill=Peak_Group)) + 
  geom_boxplot() + 
  xlab("") + ylab("Nr. of Repeat Groups") + 
  # ggtitle("Consecutive Repeat Groups Distribution", subtitle = "Only considers peaks with at least 2 repeats. Consecutive = at least 4 in a row") + 
  labs(caption = "Only peaks with at least 2 consecutive repeats were included") + 
  
  facet_wrap(~Sample, scales = "free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), plot.caption = element_text(hjust = 0.5, vjust = 141), legend.position = "top", legend.title = element_blank() ) +
  scale_fill_manual(values=cbPalette) +             # colour by cbPalette
  scale_y_continuous(limits = c(0,10), breaks=seq(0,10, by=1)) 
  

dev.off()

# # all peaks vs non-overlapping peaks
# pdf("ConsRepeatGroupCount_allVSnoOverlap_boxplot.pdf", width = 8, height = 4.5)
# 
# combined_all %>% filter(Max_Cons_Repeats > 1) %>%  # filter to remove peaks with only 1 repeat. Max_Cons_Repeats > 1
#   filter(Peak_Group != "Overlapping") %>%         # filter to remove 'All' peaks and only keep 'Overlapping' and 'Non-Overlapping'
#   
#   ggplot(.,aes(x=Sample,y=Cons_Repeat_Groups,fill=Peak_Group)) + 
#   geom_boxplot() + 
#   xlab("Peaks") + ylab("Groups of Consecutive Repeats ") + 
#   ggtitle("Consecutive Repeat Groups Distribution", 
#           subtitle = "Only considers peaks with at least 2 repeats. Consecutive = at least 4 in a row") + 
#   facet_wrap(~Sample, scales = "free_x") + 
#   theme_bw() + 
#   theme(axis.text.x = element_blank(),              # remove text from x-axis ticks
#         legend.title = element_blank()) +           # remove legend (fill) title
#   scale_fill_manual(values=cbPalette) +             # colour by cbPalette
#   scale_y_continuous(limits = c(0,10), breaks=seq(0,10, by=1))
# 
# dev.off()
# 

