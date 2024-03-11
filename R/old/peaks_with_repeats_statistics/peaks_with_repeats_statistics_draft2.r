# set working directory.
setwd("/proj/nb_storage/private/terra_rloop_project/results/figures/peaks_with_repeats_statistics/")

library("data.table")
library("ggplot2")
library("dplyr")

# read files.
peak_statistics<-fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peaks_with_repeats_statistics/peak_statistics.tsv",header = T,sep="\t")
peak_repeats<-fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peaks_with_repeats_statistics/peak_repeats.tsv",header = T,sep="\t")

# define color-blind palettes.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


######################################################################################################################
# bar plot: peak count of all vs overlapping peaks
######################################################################################################################

pdf("peakCount_allVsOverlapping_barplot.pdf", width = 8.74, height = 5.03)

peak_statistics %>% filter(Peak_Group!="Non-Overlapping") %>% 
  ggplot(.,aes(x=Sample, y=Peak_Count,fill=Peak_Group)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~Sample, scales="free_x") + 
  #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x") + 
  xlab("Peak Group") + ylab("Peak Count") + ggtitle("Peak Count") + 
  geom_text(aes(label=Peak_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), legend.title = element_blank()) + 
  scale_fill_manual(values=cbPalette)

dev.off()


# peak_statistics %>% filter(Sample == "Terra" & Peak_Group != "All") %>% 
#   # peak_repeats %>% filter(Sample == "Terra" & Peak_Group != "All" & Repeat_Group == "All") %>% 
#   ggplot(.,aes(x=Sample, y=Peak_Count,fill=Peak_Group)) + 
#   geom_bar(stat="identity", position="dodge") + 
#   # facet_wrap(~Sample, scales="free_x") + 
#   #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x") + 
#   xlab("Peak Group") + ylab("Peak Count") + ggtitle("Peak Count") + 
#   geom_text(aes(label=Peak_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
#   theme_bw() + 
#   theme(axis.text.x = element_blank(), legend.title = element_blank()) + 
#   scale_fill_manual(values=cbPalette)

######################################################################################################################
# bar plot: all overlapping peaks with min1repeat vs with min4consecutiveRepeats
######################################################################################################################

pdf("peakConsRepeatCount_overlapping_barplot.pdf", width = 8.74, height = 5.03)

peak_repeats %>% filter(Peak_Group=="Overlapping") %>% filter(grepl("(No_R|Min_1R|Min_4cR)",Repeat_Group,perl = T)) %>% 
  ggplot(.,aes(x = Sample, y=Ratio_Peak_Count, fill=factor(Repeat_Group, levels = c("No_R","Min_1R","Min_4cR")))) + 
  # ggplot(.,aes(x = factor(Peak_Count,levels = Peak_Count[order(Peak_Count, decreasing = TRUE)]), y=Ratio_Peak_Count, fill=factor(Repeat_Group, levels = c("No_R","Min_1R","Min_4cR")))) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_y_continuous(limits = c(0,100), breaks=seq(0,100, by=10)) +
  facet_wrap(~Sample, scales="free_x") + 
  # facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x") +
  xlab("") + ylab("% of Overlapping Peaks") + ggtitle("Count of Overlapping Peaks with Minimum Nr. of Consecutive Repeats") + 
  geom_text(aes(label=sprintf("%0.0f", Ratio_Peak_Count)), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  # geom_text(aes(label=Peak_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), legend.title = element_blank()) + 
  scale_fill_manual(values=cbPalette, labels = c('No repeats','Min 1 repeat','Min 4 consecutive repeats'))

dev.off()



# peak_repeats %>% filter(Peak_Group=="Overlapping") %>% filter(grepl("(All|Min_1R|Min_4cR)",Repeat_Group,perl = T)) %>% 
#   ggplot(.,aes(x=Sample, y=Peak_Count, fill=Repeat_Group)) + 
#   geom_bar(stat="identity", position="dodge") + 
#   facet_wrap(~Sample, scales="free_x") + 
#   #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x") + 
#   # facet_wrap(~factor(Sample, levels = unique(peak_repeats$Sample))) + 
#   xlab("Overlapping Peaks") + ylab("Peak Count") + ggtitle("Count of Overlapping Peaks with Minimum Nr. of Consecutive Repeats") + 
#   geom_text(aes(label=Peak_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
#   theme_bw() + 
#   theme(axis.text.x = element_blank(), legend.title = element_blank()) + 
#   scale_fill_manual(values=cbPalette, labels = c('All','Min 1 repeat','Min 4 consecutive repeats'))

######################################################################################################################
# bar plot: all non-overlapping peaks with min1repeat vs with min4consecutiveRepeats
######################################################################################################################

# library("scales")

pdf("peakConsRepeatCount_nonOverlapping_barplot.pdf", width = 8.74, height = 5.03)

peak_repeats %>% filter(Peak_Group=="Non-Overlapping") %>% filter(grepl("(No_R|Min_1R|Min_4cR)",Repeat_Group,perl = T)) %>% 
  ggplot(.,aes(x = Sample, y=Ratio_Peak_Count, fill=factor(Repeat_Group, levels = c("No_R","Min_1R","Min_4cR")))) + 
  # ggplot(.,aes(x = factor(Peak_Count,levels = Peak_Count[order(Peak_Count, decreasing = TRUE)]), y=Ratio_Peak_Count, fill=factor(Repeat_Group, levels = c("No_R","Min_1R","Min_4cR")))) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_y_continuous(limits = c(0,100), breaks=seq(0,100, by=10)) + 
  # scale_y_continuous(limits = c(0,100), breaks=seq(0,100, by=10), labels = scales::number_format(accuracy = 0.01)) + 
  facet_wrap(~Sample, scales="free_x") + 
  # facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x") +
  xlab("") + ylab("% of Non-Overlapping Peaks") + ggtitle("Count of Non-Overlapping Peaks with Minimum Nr. of Consecutive Repeats") + 
  geom_text(aes(label=sprintf("%0.1f", Ratio_Peak_Count)), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  # geom_text(aes(label=Peak_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), legend.title = element_blank()) + 
  scale_fill_manual(values=cbPalette, labels = c('No repeats','Min 1 repeat','Min 4 consecutive repeats'))

dev.off()

# peak_repeats %>% filter(Peak_Group=="Non-Overlapping") %>% filter(grepl("(All|Min_1R|Min_4cR)",Repeat_Group,perl = T)) %>%
#   ggplot(.,aes(x=Sample, y=Peak_Count,fill=Repeat_Group)) +
#   geom_bar(stat="identity", position="dodge") +
#   facet_wrap(~Sample, scales="free_x") +
#   #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x") +
#   xlab("Non-Overlapping Peaks") + ylab("Peak Count") + ggtitle("Count of Non-Overlapping Peaks with Minimum Nr. of Consecutive Repeats") +
#   geom_text(aes(label=Peak_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) +
#   theme_bw() +
#   theme(axis.text.x = element_blank(), legend.title = element_blank()) +
#   scale_fill_manual(values=cbPalette, labels = c('All','Min 1 repeat','Min 4 consecutive repeats'))



######################################################################################################################
# bar plot: mean nr. of repeats per peak group
######################################################################################################################

pdf("meanRepeatsPerPeak_barplot.pdf", width = 8.74, height = 5.03)

peak_statistics %>% select(Peak_Group,Sample,Min_1R,Repeat_Count) %>% mutate(ratio=Repeat_Count/Min_1R) %>%  
  ggplot(.,aes(x = factor(Min_1R, 
                          levels = Min_1R[order(Min_1R, decreasing = TRUE)]), 
               y=ratio,fill=Peak_Group)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~Sample, scales="free_x") + 
  #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x") + 
  
  xlab("Peaks with Repeats") + ylab("Average Nr. of Repeats per Peak") + 
  ggtitle("Mean Repeats per Peak", subtitle = "Numbers on bars: total repeat count") + 
  #geom_text(aes(label=Repeat_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), legend.title = element_blank()) + 
  #theme(legend.title = element_blank()) + 
  scale_fill_manual(values=cbPalette)

dev.off()

#########################################################################################################################################
# bar plot: repeats are highly enrihed in overlapping peaks. How many peaks contain how many repeats.

# only a small amount of overlapping peaks harbor almost half of all repeats
# in other words, almost half of all repeats are found within a very small amount of overlapping peaks
#########################################################################################################################################

pdf("repeatCountPerPeakType_perPeakNr.pdf", width = 8.74, height = 5.03)

peak_statistics %>% select(Peak_Group,Sample,Min_1R,Repeat_Count) %>% 
  ggplot(.,aes(x = factor(Min_1R, 
                          levels = Min_1R[order(Min_1R, decreasing = TRUE)]), 
               y=Repeat_Count,fill=Peak_Group)) + 
  geom_bar(stat="identity", position="dodge", alpha=0.6) + 
  facet_wrap(~Sample, scales = "free_x") + 
  #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x") + 
  xlab("Nr. of Peaks per Peak Group") + ylab("Repeat Count") + 
  ggtitle("Nr. Repeats in Different Peak Groups and Nr. Peaks Harboring Them",
          subtitle = "Inner bars represent the x-axis labels") + 
  geom_text(aes(label=Repeat_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  
  geom_bar(aes(x = factor(Min_1R, 
                          levels = Min_1R[order(Min_1R, decreasing = TRUE)]), 
               y=Min_1R), stat="identity", position="dodge") + 
  theme_bw() + 
  scale_fill_manual(values=cbPalette, name = "Peak Group")

dev.off()


