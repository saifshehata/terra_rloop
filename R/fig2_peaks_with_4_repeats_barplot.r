# set working directory.
setwd("/cluster/khiom/sshehata001/proj/terra_rloop/results/figures/Rplots")

library("data.table")
library("ggplot2")
library("dplyr")

# read files.
peak_repeats<-fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/repeat_stats//peak_repeats.tsv",header = T,sep="\t")

# define color-blind palettes.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# create plot
pdf("fig2_peaks_with_4_repeats_barplot.pdf", width = 8, height = 4.5)

peak_repeats %>% filter(Sample=="TERRA" & Peak_Group!="All") %>% filter(grepl("(No_4tR|Min_4tR)",Repeat_Group,perl = T)) %>% 
  ggplot(.,aes(x = factor(Repeat_Group, levels = c("No_4tR","Min_4tR")), y=Ratio_Peak_Count, fill=factor(Repeat_Group, levels = c("No_4tR","Min_4tR")))) + 
  # ggplot(.,aes(x = factor(Peak_Count,levels = Peak_Count[order(Peak_Count, decreasing = TRUE)]), y=Ratio_Peak_Count, fill=factor(Repeat_Group, levels = c("No_R","Min_1R","Min_4cR")))) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_y_continuous(limits = c(0,106), breaks=seq(0,100, by=10)) +
  xlab("") + ylab("Fraction of total (%)") + 
  labs(fill = "Repeat_Group", subtitle = "TERRA Peaks") + 
  
  facet_wrap(~factor(Peak_Group, levels = unique(Peak_Group))) +
  # facet_wrap(~factor(Peak_Group, levels = c("Non-Intersecting","Intersecting")), scales="free_x") + 
  geom_text(aes(label=sprintf("%0.1f", Ratio_Peak_Count)), vjust=-0.3, color="black",position = position_dodge(0.9), size=6) + 
  # geom_text(aes(label=Peak_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  geom_text(aes(label = Peak_Count, y = 0), size=6 ) +
  theme_bw() + 
  theme(plot.subtitle = element_text(vjust = -10.5, size = 18), axis.text.x = element_blank(), legend.position ="top", legend.title = element_blank(), text = element_text(size = 22)) + 
  # theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), legend.position = c(0.86, 0.85), legend.direction = "vertical", legend.background = element_rect(fill = "lightgrey",color = "lightgrey", linetype = "solid"),legend.key = element_rect(fill = "lightgrey", color = NA) ) + 
  scale_fill_manual(values=cbPalette, labels = c('Repeats < 4','Repeats >= 4 tandem'))

dev.off()
