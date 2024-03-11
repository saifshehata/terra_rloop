# set working directory.
setwd("/proj/nb_storage/private/terra_rloop_project/results/figures/peaks_with_repeats_statistics/")

library("data.table")
library("ggplot2")
library("dplyr")
#library(stringr) # to detect partial matches in dplyr::filter()
#library(scales) # to format labels in percent.
#library(tibble) # to add columns in between existing ones.

# filter by partial match using stringr.
#peak_repeats %>% filter((str_detect (V3, "overlap") | str_detect (V3, "overlap_rep1")) & !str_detect(V3,"rep4c"))

# add columns example using tibble.
#peak_statistics2 <- add_column(peak_statistics, noRepeats=peak_statistics$peakCount-peak_statistics$min1repeat, .after = "min1repeat")
#peak_statistics2 <- add_column(peak_statistics2, no4consecutiveRepeats=peak_statistics$peakCount-peak_statistics$min4consecutiveRepeats, .after = "min4consecutiveRepeats")

# read files.
#rloop_overlap<-fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/annotation/rloop_overlap_narrowPeak_ann_allRepeats.txt",header=F,sep="\t",skip = 1)
#rloop_overlap_repeats<-rloop_overlap %>% dplyr::select(V1,V22) %>% dplyr::filter(grepl("(TTAGGG|CCCTAA)",V22,perl = T))
peak_statistics<-fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peaks_with_repeats_statistics/peak_statistics.tsv",header = T,sep="\t")
peak_repeats<-fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peaks_with_repeats_statistics/peak_repeats.tsv",header = T,sep="\t")


#########################################################################################################################################
# define color-blind palettes.
#########################################################################################################################################

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
#scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)


#########################################################################################################################################
# bar plot: peak count of all vs overlapping
#########################################################################################################################################

# draft
peak_statistics %>% filter(Peaks!="Non-Overlapping") %>% 
  ggplot(.,aes(x=Sample, y=Peak_Count,fill=Peaks)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("") + ylab("Peak Count") + ggtitle("Count of All Peaks vs. Overlapping Peaks") + 
  geom_text(aes(label=Peak_Count), vjust=-0.5, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_minimal() + 
  scale_fill_manual(values=cbPalette)

# final
pdf("peakCount_allVsOverlapping.pdf", width = 8.74, height = 5.03)

peak_statistics %>% filter(Peaks!="Non-Overlapping") %>% 
  ggplot(.,aes(x=Sample, y=Peak_Count,fill=Peaks)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~Sample, scales="free_x") + 
  #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x") + 
  xlab("") + ylab("Peak Count") + ggtitle("Count of All Peaks vs. Overlapping Peaks") + 
  geom_text(aes(label=Peak_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=cbPalette)

dev.off()

#########################################################################################################################################
# bar plot for overlapping peaks: all vs with min1repeat vs with min4consecutiveRepeats
#########################################################################################################################################

# draft
peak_repeats %>% filter(Peaks=="Overlapping") %>% filter(grepl("(All|Min_1R|Min_4cR)",Repeats,perl = T)) %>% 
  ggplot(.,aes(x=Sample, y=Peak_Count, fill=Repeats)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("") + ylab("Peak Count") + ggtitle("Count of Overlapping Peaks with Minimum Nr. of Consecutive Repeats") + 
  geom_text(aes(label=Peak_Count), vjust=-0.5, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_minimal() + 
  scale_fill_manual(values=cbPalette)

# final
pdf("peakConsRepeatCount_overlapping.pdf", width = 8.74, height = 5.03)

peak_repeats %>% filter(Peaks=="Overlapping") %>% filter(grepl("(All|Min_1R|Min_4cR)",Repeats,perl = T)) %>% 
  ggplot(.,aes(x=Sample, y=Peak_Count, fill=Repeats)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~Sample, scales="free_x") + 
  #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x") + 
  xlab("") + ylab("Peak Count") + ggtitle("Count of Overlapping Peaks with Minimum Nr. of Consecutive Repeats") + 
  geom_text(aes(label=Peak_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=cbPalette)

dev.off()

#########################################################################################################################################
# bar plot for non-overlapping: all vs with min1repeat vs with min4consecutiveRepeats
#########################################################################################################################################

# draft
peak_repeats %>% filter(Peaks=="Non-Overlapping") %>% filter(grepl("(All|Min_1R|Min_4cR)",Repeats,perl = T)) %>% 
  ggplot(.,aes(x=Sample, y=Peak_Count,fill=Repeats)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("") + ylab("Peak Count") + ggtitle("Count of Non-Overlapping Peaks with Minimum Nr. of Consecutive Repeats") + 
  geom_text(aes(label=Peak_Count), vjust=-0.5, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_minimal() + 
  scale_fill_manual(values=cbPalette)

# final
pdf("peakConsRepeatCount_nonOverlapping.pdf", width = 8.74, height = 5.03)

peak_repeats %>% filter(Peaks=="Non-Overlapping") %>% filter(grepl("(All|Min_1R|Min_4cR)",Repeats,perl = T)) %>% 
  ggplot(.,aes(x=Sample, y=Peak_Count,fill=Repeats)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~Sample, scales="free_x") + 
  #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x") + 
  xlab("") + ylab("Peak Count") + ggtitle("Count of Non-Overlapping Peaks with Minimum Nr. of Consecutive Repeats") + 
  geom_text(aes(label=Peak_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=cbPalette)

dev.off()

#########################################################################################################################################
# bar plot for total repeat count in overlapping vs non-overlapping peaks
# shows that repeats are highly enrihed in overlapping peaks
# only a small amount of overlapping peaks harbor almost half of all repeats
# in other words, almost half of all repeats are found within a very small amount of overlapping peaks
#########################################################################################################################################

# draft
peak_statistics %>% select(Peaks,Sample,Min_1R,Repeat_Count) %>% mutate(Ratio=Repeat_Count/Min_1R) %>%  
  ggplot(.,aes(x=Sample, y=Ratio,fill=Peaks)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("") + ylab("Average Nr. of Repeats per Peak") + 
  ggtitle("Mean Repeats per Peak", subtitle = "Black numbers: total repeat count. White numbers: Nr. of peaks harboring the repeats") + 
  geom_text(aes(label=Repeat_Count), vjust=-0.5, color="black",position = position_dodge(0.9), size=3.5) + 
  geom_text(aes(label=Min_1R), vjust=2, color="white",position = position_dodge(0.9), size=3.5, fontface="bold") + 
  theme_minimal()  + 
  scale_fill_manual(values=cbPalette)

# final

pdf("meanRepeatsPerPeak.pdf", width = 8.74, height = 5.03)

peak_statistics %>% select(Peaks,Sample,Min_1R,Repeat_Count) %>% mutate(Ratio=Repeat_Count/Min_1R) %>%  
  ggplot(.,aes(x = factor(Min_1R, 
                          levels = Min_1R[order(Min_1R, decreasing = TRUE)]), 
               y=Ratio,fill=Peaks)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~Sample, scales="free_x") + 
  #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x") + 
  
  xlab("Nr. of Peaks Containing the Repeats") + ylab("Average Nr. of Repeats per Peak") + 
  ggtitle("Mean Repeats per Peak", subtitle = "Numbers on bars: total repeat count") + 
  geom_text(aes(label=Repeat_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_bw() + 
  scale_fill_manual(values=cbPalette)

dev.off()


# draft 1
peak_statistics %>% select(Peaks,Sample,Min_1R,Repeat_Count) %>%  
  ggplot(.,aes(x=Sample, y=Repeat_Count,fill=Peaks)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("") + ylab("Average Nr. of Repeats per Peak") + 
  ggtitle("Mean Repeats per Peak", subtitle = "Black numbers: total repeat count. White numbers: Nr. of peaks harboring the repeats") + 
  geom_text(aes(label=Repeat_Count), vjust=-0.5, color="black",position = position_dodge(0.9), size=3.5) + 
  geom_text(aes(label=Min_1R), vjust=2, color="white",position = position_dodge(0.9), size=3.5, fontface="bold") + 
  theme_minimal()  + 
  scale_fill_manual(values=cbPalette)

# draft 2
peak_statistics %>% select(Peaks,Sample,Min_1R,Repeat_Count) %>%  
  ggplot(.,aes(x=Sample, y=Repeat_Count,fill=Peaks)) + 
  geom_bar(stat="identity", position="dodge", alpha=0.6) + 
  xlab("") + ylab("Average Nr. of Repeats per Peak") + 
  ggtitle("Mean Repeats per Peak", subtitle = "Black numbers: total repeat count. White numbers: Nr. of peaks harboring the repeats") + 
  geom_text(aes(label=Repeat_Count), vjust=-0.5, color="black",position = position_dodge(0.9), size=3.5) + 
  geom_bar(aes(x=Sample,y=Min_1R), stat="identity", position="dodge") + 
  geom_text(aes(x=Sample,y=Min_1R,label=Min_1R), vjust=-0.3, color="white",position = position_dodge(0.9), size=3.5, fontface="bold") + 
  theme_minimal()  + 
  scale_fill_manual(values=cbPalette)

# draft 3
peak_statistics %>% select(Peaks,Sample,Min_1R,Repeat_Count) %>% 
  ggplot(.,aes(x=Sample, y=Repeat_Count,fill=Peaks)) + 
  geom_bar(stat="identity", position="dodge", alpha=0.6) + 
  xlab("") + ylab("Repeat Count") + 
  ggtitle("Nr. Repeats in Different Peak Groups and Nr. Peaks Harboring Them", subtitle = "Black numbers: total repeat count. White numbers: Nr. of peaks harboring the repeats") + 
  geom_text(aes(label=Repeat_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) +
  geom_bar(aes(x=Sample,y=Min_1R), stat="identity", position="dodge") + 
  geom_text(aes(x=Sample,y=Min_1R,label=Min_1R), vjust=-0.3, color="white",position = position_dodge(0.9), size=3.5, fontface="bold") + 
  facet_wrap(~Sample, scales="free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values=cbPalette)

# draft 4
peak_statistics %>% select(Peaks,Sample,Min_1R,Repeat_Count) %>% 
  ggplot(.,aes(x=as.character(Min_1R), y=Repeat_Count,fill=Peaks)) + 
  geom_bar(stat="identity", position="dodge", alpha=0.6) + 
  facet_wrap(~Sample, scales="free_x") + 
  xlab("Nr. of Peaks Containing the Repeats") + ylab("Repeat Count") + 
  ggtitle("Nr. Repeats in Different Peak Groups and Nr. Peaks Harboring Them") + 
  geom_text(aes(label=Repeat_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) +
  geom_bar(aes(x=as.character(Min_1R),y=Min_1R), stat="identity", position="dodge") + 
  theme_bw() + 
  scale_fill_manual(values=cbPalette)

# draft 5
peak_statistics %>% select(Peaks,Sample,Min_1R,Repeat_Count) %>% 
  ggplot(.,aes(x = factor(Min_1R, levels = c("2718","2454","264","3669","3433","236")), y=Repeat_Count,fill=Peaks)) + 
  geom_bar(stat="identity", position="dodge", alpha=0.6) + 
  facet_wrap(~Sample, scales="free_x") + 
  xlab("Nr. of Peaks Containing the Repeats") + ylab("Repeat Count") + 
  ggtitle("Nr. Repeats in Different Peak Groups and Nr. Peaks Harboring Them") + 
  geom_text(aes(label=Repeat_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) +
  geom_bar(aes(x = as.character(factor(Min_1R, levels = c("2718","2454","264","3669","3433","236"))),y=Min_1R), stat="identity", position="dodge") + 
  theme_bw() + 
  scale_fill_manual(values=cbPalette)


# final

pdf("repeatCountPerPeakType_perPeakNr.pdf", width = 8.74, height = 5.03)

peak_statistics %>% select(Peaks,Sample,Min_1R,Repeat_Count) %>% 
  ggplot(.,aes(x = factor(Min_1R, 
                          levels = Min_1R[order(Min_1R, decreasing = TRUE)]), 
               y=Repeat_Count,fill=Peaks)) + 
  geom_bar(stat="identity", position="dodge", alpha=0.6) + 
  facet_wrap(~Sample, scales = "free_x") + 
  #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x") + 
  
  xlab("Nr. of Peaks Per Peak Group") + ylab("Repeat Count") + 
  ggtitle("Nr. Repeats in Different Peak Groups and Nr. Peaks Harboring Them") + 
  geom_text(aes(label=Repeat_Count), vjust=-0.3, color="black",position = position_dodge(0.9), size=3.5) + 
  
  geom_bar(aes(x = factor(Min_1R, 
                          levels = Min_1R[order(Min_1R, decreasing = TRUE)]), 
               y=Min_1R), stat="identity", position="dodge") + 
  theme_bw() + 
  scale_fill_manual(values=cbPalette)

dev.off()

#facet_wrap(facets = ~reorder(Sample, -Repeat_Count), scales="free_x") + 
#facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x") + 
# facet_wrap(~Sample)
#ordered(name, levels = rev(sort(unique(name))))

#########################################################################################################################################
# END
#########################################################################################################################################


###

#position = position_fill(vjust=0,reverse=T), size=3.5

#x = factor(V1, levels = unique(terra_all$V1))
#factor(data3$x,                                    # Factor levels in decreasing order
#       levels = data3$x[order(data3$y, decreasing = TRUE)])
#x = as.character(factor(Min_1R, levels = order(Min_1R,decreasing = T)))



# to add colours.
#scale_fill_manual(values=c("red","green","blue")) + 

##################################### TESTING ####################################################################

peak_statistics %>% select(Peaks,Sample,Min_1R,Repeat_Count) %>% mutate(Ratio=Min_1R/Repeat_Count*100) %>%  
  ggplot(.,aes(x=Sample, y=Ratio,fill=Peaks)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("") + ylab("% Nr. Peaks / Nr. Repeats") + ggtitle("Proportion of Peaks with Repeats") + 
  geom_text(aes(label=Repeat_Count), vjust=-0.5, color="black",position = position_dodge(0.9), size=3.5) + 
  geom_text(aes(label=Min_1R), vjust=2, color="white",position = position_dodge(0.9), size=3.5) + 
  theme_minimal() + 
  scale_fill_manual(values=cbPalette)


peak_statistics %>% select(Peaks,Sample,Min_1R,Repeat_Count) %>% 
  ggplot(.,aes(x=Peaks, y=Repeat_Count ,fill=Sample)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("") + ylab("Peak Count") + ggtitle("Overlapping Peaks with Repeats") + 
  geom_text(aes(label=Repeat_Count), vjust=-0.2, color="black",position = position_dodge(0.9), size=3.5) +
  geom_text(aes(label=Min_1R), vjust=2, color="white",position = position_dodge(0.9), size=3.5) + 
  theme_minimal() + 
  scale_fill_manual(values=cbPalette)

peak_statistics %>% select(Peaks,Sample,Min_1R,Repeat_Count) %>% 
  ggplot(.,aes(x=Peaks, y=Repeat_Count ,fill=Peaks)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~Sample) + 
  xlab("") + ylab("Peak Count") + ggtitle("Overlapping Peaks with Repeats") + 
  geom_text(aes(label=Repeat_Count), vjust=-0.2, color="black",position = position_dodge(0.9), size=3.5) +
  geom_text(aes(label=Min_1R), vjust=-1, color="white",position = position_fill(vjust=0,reverse=T), size=3.5, fontface="bold") + 
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=cbPalette)


peak_statistics %>% select(Peaks,Sample,Min_1R,Repeat_Count) %>% 
  ggplot(.,aes(x=Peaks, y=Repeat_Count ,fill=Sample)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("") + ylab("Peak Count") + ggtitle("Overlapping Peaks with Repeats") + 
  geom_text(aes(label=Repeat_Count), vjust=-0.2, color="black",position = position_dodge(0.9), size=3.5) +
  geom_text(aes(label=Min_1R), vjust=-2.5, color="white",position = position_fill(vjust=0,reverse=T), size=3.5) + 
  theme_minimal() 
+ 
  facet_wrap(~sample)



# testing.
# prepare data to use. filter and select columns, then transpose the data as a matrix and re-convert it to a data frame.
overlapping_min1repeat <- peak_statistics2 %>% filter(peaks=="overlap") %>% select(.,sample,min1repeat,noRepeats) %>% 
  t() %>% as.data.frame() %>% .[-1,]

# convert 1st row to column names then remove 1st row.
names(overlapping_min1repeat) <- lapply(overlapping_min1repeat[1, ], as.character)
overlapping_min1repeat <- overlapping_min1repeat[-1,]


# testing.

ggplot(overlapping_min1repeat,aes(x=c("Terra","R-loop"), y=terra)) + 
  geom_bar(stat="identity", position="stack") 
+ 
  xlab("") + ylab("Peak Count") + ggtitle("Overlapping vs. Non-overlapping Peaks") + 
  geom_text(aes(label=peakCount), vjust=0, color="black",position = position_fill(vjust=0,reverse=T), size=3.5) + 
  theme_minimal()







# stacked bar plot.
peak_statistics %>% filter(peaks!="allPeaks") %>% 
  ggplot(.,aes(x=sample, y=peakCount,fill=peaks)) + 
  geom_bar(stat="identity", position=position_fill(reverse=T)) + 
  xlab("") + ylab("Peak Count") + ggtitle("Overlapping vs. Non-overlapping Peaks") + 
  geom_text(aes(label=peakCount), vjust=0, color="black",position = position_fill(vjust=0,reverse=T), size=3.5) + 
  theme_minimal()

# pie chart for overlapping vs non-overlapping peaks.
peak_statistics %>% filter(peaks!="allPeaks",sample=="terra") %>% arrange(desc(peaks)) %>%
  mutate(prop = peakCount / sum(.$peakCount) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )%>% ggplot(.,aes(x="", y=prop, fill=peaks)) + 
  geom_bar(stat="identity",width=1,color="white") + 
  coord_polar("y", start=0) +
  theme_void()+
  geom_text(aes(y = ypos, label = peakCount), color = "white", size=5) +
  scale_fill_brewer(palette="Set1")  



peak_statistics %>% ggplot(.,aes(x=sample, y=min1repeat,fill=peaks)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("Sample") + ylab("Peak Count") + ggtitle("test") + 
  geom_text(aes(label=min1repeat), vjust=-0.2, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_minimal()

peak_statistics %>% ggplot(.,aes(x=sample, y=min4consecutiveRepeats,fill=peaks)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("Sample") + ylab("Peak Count") + ggtitle("test") + 
  geom_text(aes(label=min4consecutiveRepeats), vjust=-0.2, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_minimal()

peak_statistics %>% ggplot(.,aes(x=sample, y=repeatCount,fill=peaks)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("Sample") + ylab("Peak Count") + ggtitle("test") + 
  geom_text(aes(label=repeatCount), vjust=-0.2, color="black",position = position_dodge(0.9), size=3.5) + 
  theme_minimal()

+ 
  scale_fill_brewer(palette="Paired")
