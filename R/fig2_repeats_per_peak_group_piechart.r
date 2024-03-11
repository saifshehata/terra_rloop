# set working directory.
setwd("/cluster/khiom/sshehata001/proj/terra_rloop/results/figures/Rplots")

library("data.table")
library("ggplot2")
library("dplyr")
library(scales) # for the function percent()

# read files.
peak_statistics<-fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/repeat_stats/peak_statistics.tsv",header = T,sep="\t")

# define color-blind palettes.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

######################################################################################################################
# pie chart: percent repeat enrichment per peak group
######################################################################################################################

pdf("fig2_repeats_per_peak_group_piechart.pdf", width = 4, height = 4.5)

peak_statistics[order(peak_statistics$Peak_Group)] %>% 
  filter(Sample == "TERRA" & Peak_Group != "All") %>%  
  select(Sample,Peak_Group,Min_1R,Ratio_Repeat_Count) %>% 
  mutate(prop=cumsum(Ratio_Repeat_Count/sum(Ratio_Repeat_Count)*50)) %>% 
  mutate(ypos = cumsum(prop)+ 0.05*prop ) %>% 
  
  ggplot(.,aes(x="", y=Ratio_Repeat_Count, fill=factor(Peak_Group, levels = c("Non-Intersecting", "Intersecting"))) ) + 
  geom_bar(stat = "identity", position = "stack", color="white") + 
  ggtitle("% total repeats", subtitle = "TERRA Peaks") +
  coord_polar("y") +
  theme_minimal() + 
  theme_void() + 
  geom_text(aes(y = ypos, label = paste0(round(Ratio_Repeat_Count),"%" )), color = "black", size=5.6 ) + 
  # geom_text(aes(y = ypos, label = paste(percent(Ratio_Repeat_Count/100) )), color = "black", size=5.6 ) + 
  geom_text(aes(y = ypos, label = paste0("\n","\n", "(",as.character(Min_1R), " peaks",")" )), color = "black", size=4.56) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -83, size = 16), plot.subtitle = element_text(hjust = 0.5, vjust = -13, size = 16), legend.position = "top", legend.title = element_blank(), text = element_text(size = 17.6)) +
  scale_fill_manual(values=cbPalette)

dev.off()

# peak_statistics[order(peak_statistics$Peak_Group)] %>%
#   filter(Sample == "TERRA" & Peak_Group != "All") %>%
#   select(Sample,Peak_Group,Min_1R,Ratio_Repeat_Count) %>%
#   mutate(prop=cumsum(Ratio_Repeat_Count/sum(Ratio_Repeat_Count)*50)) %>%
#   mutate(ypos = cumsum(prop)+ 0.05*prop ) %>%
# 
#   ggplot(.,aes(x="", y=Ratio_Repeat_Count, fill=factor(Peak_Group, levels = c("Non-Intersecting", "Intersecting"))) ) +
#   geom_bar(stat = "identity", position = "stack", color="white") +
#   ggtitle("% total repeats", subtitle = "TERRA Peaks") +
#   coord_polar("y") +
#   theme_minimal() +
#   theme_void() +
#   geom_text(aes(y = ypos, label = paste(percent(Ratio_Repeat_Count/100) )), color = "black", size=7 ) +
#   geom_text(aes(y = ypos, label = paste0("\n","\n", "(",as.character(Min_1R), " peaks",")" )), color = "black", size=5.7) +
#   theme(plot.title = element_text(hjust = 0.5, vjust = -67, size = 20), plot.subtitle = element_text(hjust = 0.5, vjust = -9, size = 20), legend.position = "top", legend.title = element_blank(), text = element_text(size = 22)) +
#   scale_fill_manual(values=cbPalette)


