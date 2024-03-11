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

terra_intersect_nr_repeats <- as.numeric((peak_statistics %>% filter(Sample == "TERRA" & Peak_Group == "Intersecting")) %>% select(Ratio_Repeat_Count))
terra_no_intersect_nr_repeats <- as.numeric((peak_statistics %>% filter(Sample == "TERRA" & Peak_Group == "Non-Intersecting")) %>% select(Ratio_Repeat_Count))




######################################################################################################################
# pie chart: percent repeat enrichment per peak group
######################################################################################################################
pdf("fig2_repeats_per_peak_group_piechart_simple.pdf", width = 8, height = 4.5)

pie(c(terra_no_intersect_nr_repeats, terra_intersect_nr_repeats), labels = c("Non-Intersecting", "Intersecting"), border = "white", col = cbPalette, main = "Fraction of total repeats", radius = 1, init.angle = 90)

dev.off()


##### testing
# library(scales) # for the function percent()
# 
# 
# pdf("Fig2_percent_total_repeats_terra_pie.pdf", width = 8, height = 4.5)
# 
# peak_statistics %>% filter(Sample == "TERRA" & Peak_Group != "All") %>%  
#   select(Sample, Peak_Group, Min_1R, Ratio_Repeat_Count) %>% 
#   # arrange(desc(Ratio_Repeat_Count)) %>% 
#   mutate(prop=cumsum(Ratio_Repeat_Count/sum(Ratio_Repeat_Count)*50)) %>% 
#   mutate(ypos = cumsum(prop)+ 0.05*prop ) %>%  
#   
#   ggplot(.,aes(x="", y=Ratio_Repeat_Count, fill=Peak_Group)) + 
#   geom_bar(stat = "identity", position = "stack", color="white") + 
#   ggtitle("Repeat Enrichment", subtitle = "TERRA Peaks") +
#   # labs(title = "Repeat Enrichment (%)", subtitle = "Terra Peaks") + 
#   coord_polar("y") +
#   theme_minimal() + 
#   theme_void() + 
#   geom_text(aes(y = ypos, label = paste(percent(Ratio_Repeat_Count/100))), color = "black", size=6 ) + 
#   geom_text(aes(y = ypos, label = paste0("\n","\n", "(",as.character(Min_1R), " peaks",")" )), color = "black", size=3.5) +
#   theme(plot.title = element_text(hjust = 0.5, vjust = -110), plot.subtitle = element_text(hjust = 0.5, vjust = -15), legend.position = "top", legend.title = element_blank()) +
#   #   theme(plot.title = element_text(hjust = 0.5, vjust = -110), plot.subtitle = element_text(hjust = 0.5), legend.position = c(0.5,0.97), legend.direction = "vertical") +
#   scale_fill_manual(values=cbPalette)
# 
# dev.off()
# 
# 
# 
# peak_statistics %>% filter(Sample == "TERRA" & Peak_Group != "All") %>% 
#   ggplot(.,aes(x=factor(Peak_Group, levels = unique(Peak_Group)), y=Repeat_Count, fill=factor(Peak_Group, levels = unique(Peak_Group)))) + geom_bar(stat = "identity") + 
#   facet_wrap(~factor(Sample), scales="free_x") + 
#   theme(axis.text.x = element_blank(), legend.position = "top", legend.title = element_blank(), text = element_text(size = 22)) + 
#   scale_fill_manual(values=cbPalette)
# 
# 
# 
# 
# ##############################################################################################
# # bar plot: peak count of all vs overlapping peaks
# ##############################################################################################
# 
# # pie chart for overlapping vs non-overlapping peaks.
# peak_statistics %>% filter(Peak_Group != "All", Sample=="TERRA") %>% arrange(desc(peaks)) %>%
#   mutate(prop = Peak_Count / sum(.$Peak_Count) *100) %>%
#   mutate(ypos = cumsum(prop)- 0.5*prop )%>% ggplot(.,aes(x="", y=prop, fill=peaks)) + 
#   geom_bar(stat="identity",width=1,color="white") + 
#   coord_polar("y", start=0) +
#   theme_void()+
#   geom_text(aes(y = ypos, label = Peak_Count), color = "white", size=5) +
#   scale_fill_brewer(palette="Set1")  
# 
# 
# 
# peak_statistics %>% ggplot(.,aes(x=Sample, y=Min_1R,fill=Peak_Group)) + 
#   geom_bar(stat="identity", position="dodge") + 
#   xlab("Sample") + ylab("Peak Count") + ggtitle("test") + 
#   geom_text(aes(label=Mn_1R), vjust=-0.2, color="black",position = position_dodge(0.9), size=3.5) + 
#   theme_minimal()
# 
# 
# 
# peak_statistics %>% ggplot(.,aes(x=sample, y=min4consecutiveRepeats,fill=peaks)) + 
#   geom_bar(stat="identity", position="dodge") + 
#   xlab("Sample") + ylab("Peak Count") + ggtitle("test") + 
#   geom_text(aes(label=min4consecutiveRepeats), vjust=-0.2, color="black",position = position_dodge(0.9), size=3.5) + 
#   theme_minimal()
# 
# peak_statistics %>% ggplot(.,aes(x=sample, y=repeatCount,fill=peaks)) + 
#   geom_bar(stat="identity", position="dodge") + 
#   xlab("Sample") + ylab("Peak Count") + ggtitle("test") + 
#   geom_text(aes(label=repeatCount), vjust=-0.2, color="black",position = position_dodge(0.9), size=3.5) + 
#   theme_minimal()
# 
# + 
#   scale_fill_brewer(palette="Paired")
