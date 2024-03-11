# set working directory.
setwd("/cluster/khiom/sshehata001/proj/terra_rloop/results/figures/Rplots")

# load packages
library("data.table")
library("ggplot2")
library("dplyr")

# read files.
peak_statistics<-fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/repeat_stats/peak_statistics.tsv",header = T,sep="\t")

# define color-blind palettes.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


##############################################################################################
# bar plot: peak count of all vs overlapping peaks
##############################################################################################

pdf("fig1_peak_count_barplot.pdf", width = 8, height = 4.5)

peak_statistics %>% filter(Peak_Group != "All") %>% 
  ggplot(.,aes(x=factor(Peak_Group, levels = unique(Peak_Group)), y=Peak_Count,fill=factor(Peak_Group, levels = unique(Peak_Group)))) + 
  #ggplot(.,aes(x=Peak_Group, y=Peak_Count,fill=Peak_Group)) + 
  
  geom_bar(stat="identity", position="dodge") + 
  xlab("") + ylab("Nr. of Peaks") + 
  
  facet_wrap(~factor(Sample, levels = unique(Sample)), scales="free_x") +
  # facet_wrap(~factor(Sample, levels = c("Terra","R-loop")), scales="free_x") + 
  geom_text(aes(label=Peak_Count), vjust=-0.2, color="black",position = position_dodge(0.9), size=6) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,30000), breaks=seq(0,30000, by=10000)) + # no strict need for this, but the Peak_Count labels on top of the bars will be slightly hidden otherwise. Adding the breaks=seq() part does not make any difference in this case, but keep it here in case needs to be fine-tuned later
  theme(axis.text.x = element_blank(), legend.position = "top", legend.title = element_blank(), text = element_text(size = 22)) + 
  scale_fill_manual(values=cbPalette)

dev.off()
