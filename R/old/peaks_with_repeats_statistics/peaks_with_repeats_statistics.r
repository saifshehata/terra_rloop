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

pdf("Fig1A_peakCount_barplot.pdf", width = 8, height = 4.5)

peak_statistics %>% 
  ggplot(.,aes(x=Peak_Group, y=Peak_Count,fill=Peak_Group)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("") + ylab("Nr. of Peaks") + 

  facet_wrap(~factor(Sample, levels = unique(Sample)), scales="free_x") +
  # facet_wrap(~factor(Sample, levels = c("Terra","R-loop")), scales="free_x") + 
  geom_text(aes(label=Peak_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), legend.position = "top", legend.title = element_blank()) + 
  # theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), legend.position = c(0.9, 0.85),legend.direction = "vertical", legend.background = element_rect(fill = "lightgrey",color = "lightgrey", linetype = "solid"),legend.key = element_rect(fill = "lightgrey", color = NA) ) + 
  scale_fill_manual(values=cbPalette)

dev.off()


######################################################################################################################
# bar plot: all overlapping peaks with min1repeat vs with min4consecutiveRepeats
######################################################################################################################

pdf("Fig2_percPeaksWithMin1Repeat_terra_barplot.pdf", width = 8, height = 4.5)

peak_repeats %>% filter(Sample=="Terra" & Peak_Group!="All") %>% filter(grepl("(No_R|Min_1R)",Repeat_Group,perl = T)) %>% 
  ggplot(.,aes(x = factor(Repeat_Group, levels = c("No_R","Min_1R")), y=Ratio_Peak_Count, fill=factor(Repeat_Group, levels = c("No_R","Min_1R")))) + 
  # ggplot(.,aes(x = factor(Peak_Count,levels = Peak_Count[order(Peak_Count, decreasing = TRUE)]), y=Ratio_Peak_Count, fill=factor(Repeat_Group, levels = c("No_R","Min_1R","Min_4cR")))) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_y_continuous(limits = c(0,100), breaks=seq(0,100, by=10)) +
    xlab("") + ylab("Fraction of total (%)") + 
  labs(fill = "Repeat_Group", subtitle = "Terra Peaks") + 
  
  facet_wrap(~Peak_Group) +
  # facet_wrap(~factor(Peak_Group, levels = c("Overlapping","Non-Overlapping")), scales="free_x") + 
  
  geom_text(aes(label=sprintf("%0.1f", Ratio_Peak_Count)), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  # geom_text(aes(label=Peak_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) +
  theme_bw() + 
  theme(plot.subtitle = element_text(vjust = -13), axis.text.x = element_blank(), legend.position ="top", legend.title = element_blank()) + 
  # theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), legend.position = c(0.86, 0.85), legend.direction = "vertical", legend.background = element_rect(fill = "lightgrey",color = "lightgrey", linetype = "solid"),legend.key = element_rect(fill = "lightgrey", color = NA) ) + 
  scale_fill_manual(values=cbPalette, labels = c('No Repeats','Minimum 1 Repeat'))

dev.off()


pdf("Fig2_percPeaksWithMin4Repeat_terra_barplot.pdf", width = 8, height = 4.5)

peak_repeats %>% filter(Sample=="Terra" & Peak_Group!="All") %>% filter(grepl("(No_4cR|Min_4cR)",Repeat_Group,perl = T)) %>% 
  ggplot(.,aes(x = factor(Repeat_Group, levels = c("No_4cR","Min_4cR")), y=Ratio_Peak_Count, fill=factor(Repeat_Group, levels = c("No_4cR","Min_4cR")))) + 
  # ggplot(.,aes(x = factor(Peak_Count,levels = Peak_Count[order(Peak_Count, decreasing = TRUE)]), y=Ratio_Peak_Count, fill=factor(Repeat_Group, levels = c("No_R","Min_1R","Min_4cR")))) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_y_continuous(limits = c(0,100), breaks=seq(0,100, by=10)) +
  xlab("") + ylab("Fraction of total (%)") + 
  labs(fill = "Repeat_Group", subtitle = "Terra Peaks") + 
  
  facet_wrap(~Peak_Group) +
  # facet_wrap(~factor(Peak_Group, levels = c("Overlapping","Non-Overlapping")), scales="free_x") + 
  
  geom_text(aes(label=sprintf("%0.1f", Ratio_Peak_Count)), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  # geom_text(aes(label=Peak_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_bw() + 
  theme(plot.subtitle = element_text(vjust = -13), axis.text.x = element_blank(), legend.position ="top", legend.title = element_blank()) + 
  # theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), legend.position = c(0.86, 0.85), legend.direction = "vertical", legend.background = element_rect(fill = "lightgrey",color = "lightgrey", linetype = "solid"),legend.key = element_rect(fill = "lightgrey", color = NA) ) + 
  scale_fill_manual(values=cbPalette, labels = c('Lacking 4 Consecutive Repeats','Minimum 4 Consecutive Repeats'))

dev.off()


######################################################################################################################
# bar plot: mean nr. of repeats per peak group
######################################################################################################################

pdf("Fig2_avRepeatsPerPeak_terra_barplot.pdf", width = 4, height = 4.5)

peak_statistics %>% filter(Sample=="Terra" & Peak_Group!="All") %>% 
  select(Peak_Group,Sample,Min_1R,Repeat_Count) %>% mutate(ratio=Repeat_Count/Min_1R) %>%  
  ggplot(.,aes(x = factor(Min_1R, 
                          levels = Min_1R[order(Min_1R, decreasing = TRUE)]), 
               y=ratio,fill=Peak_Group)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("") + ylab("Average Nr. of Repeats") + 

  facet_wrap(~Sample) + 

  #geom_text(aes(label=Repeat_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_bw() + 
  theme(legend.position = "top", legend.title = element_blank()) + 
  # theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), legend.position = c(0.18, 0.85),legend.direction = "vertical", legend.background = element_rect(fill = "lightgrey",color = "lightgrey", linetype = "solid"),legend.key = element_rect(fill = "lightgrey", color = NA) ) + 
  #theme(legend.title = element_blank()) + 
  scale_fill_manual(values=cbPalette, labels = c("Non-Overlapping", "Overlapping"))

dev.off()

# how to add error bars.
# ggplot(diamonds, aes(cut, price, fill = color)) +
#   stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
#   stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")

######################################################################################################################
# pie chart: percent repeat enrichment per peak group
######################################################################################################################

library(scales) # for the function percent()


pdf("Fig2_percRepeatsPerPeakGroup_terra.pdf", width = 8, height = 4.5)

peak_statistics %>% filter(Sample == "Terra" & Peak_Group != "All") %>%  
  select(Sample,Peak_Group,Min_1R,Ratio_Repeat_Count) %>% 
  # arrange(desc(Ratio_Repeat_Count)) %>% 
  mutate(prop=cumsum(Ratio_Repeat_Count/sum(Ratio_Repeat_Count)*50)) %>% 
  mutate(ypos = cumsum(prop)+ 0.05*prop ) %>%  
  
  ggplot(.,aes(x="", y=Ratio_Repeat_Count, fill=Peak_Group)) + 
  geom_bar(stat = "identity", position = "stack", color="white") + 
  ggtitle("Repeat Enrichment", subtitle = "Terra Peaks") +
  # labs(title = "Repeat Enrichment (%)", subtitle = "Terra Peaks") + 
  coord_polar("y") +
  theme_minimal() + 
  theme_void() + 
  geom_text(aes(y = ypos, label = paste(percent(Ratio_Repeat_Count/100))), color = "black", size=6 ) + 
  geom_text(aes(y = ypos, label = paste0("\n","\n", "(",as.character(Min_1R), " peaks",")" )), color = "black", size=3.5) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -110), plot.subtitle = element_text(hjust = 0.5, vjust = -15), legend.position = "top", legend.title = element_blank()) +
  #   theme(plot.title = element_text(hjust = 0.5, vjust = -110), plot.subtitle = element_text(hjust = 0.5), legend.position = c(0.5,0.97), legend.direction = "vertical") +
  scale_fill_manual(values=cbPalette, labels = c("Non-Overlapping", "Overlapping"))
  
dev.off()

