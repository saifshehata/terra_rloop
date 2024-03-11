# set working directory.
setwd("/proj/nb_storage/private/terra_rloop_project/results/figures/")

library("data.table")
library("ggplot2")
library("dplyr")

# read files.
peak_statistics<-fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/repeat_stats/peak_statistics.tsv",header = T,sep="\t")

# define color-blind palettes.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

######################################################################################################################
# bar plot: mean nr. of repeats per peak group
######################################################################################################################

pdf("fig2_mean_repeats_per_peak_barplot.pdf", width = 4, height = 4.5)

peak_statistics %>% filter(Sample=="TERRA" & Peak_Group!="All") %>% 
  select(Peak_Group,Sample,Min_1R,Repeat_Count) %>% mutate(ratio=Repeat_Count/Min_1R) %>%  
  ggplot(.,aes(x = factor(Min_1R, levels = Min_1R[order(Min_1R, decreasing = TRUE)]), 
               y=ratio,fill=Peak_Group)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("") + ylab("Average rr. of repeats per peak") + 
  facet_wrap(~Sample) + 
  #geom_text(aes(label=Repeat_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_bw() + 
  theme(legend.position = "top", legend.title = element_blank(), text = element_text(size = 22) ) + 
  # theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), legend.position = c(0.18, 0.85),legend.direction = "vertical", legend.background = element_rect(fill = "lightgrey",color = "lightgrey", linetype = "solid"),legend.key = element_rect(fill = "lightgrey", color = NA) ) + 
  #theme(legend.title = element_blank()) + 
  scale_fill_manual(values=c(cbPalette[2], cbPalette[1])) + 
  guides(fill = guide_legend(reverse=TRUE))

dev.off()