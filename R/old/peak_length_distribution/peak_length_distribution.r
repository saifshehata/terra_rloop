# set working directory.
setwd("/proj/nb_storage/private/terra_rloop_project/results/figures/peak_length_distribution/")

# load packages.
library("data.table")
library("ggplot2")
library("dplyr")

# read bed files with peak coordinates.
terra_all <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_peaks_peakLength.bed", header = F,sep="\t")
rloop_all <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/rloop_peaks_peakLength.bed", header = F,sep="\t")
terra_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_overlap_peakLength.bed", header = F,sep="\t")
rloop_overlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/rloop_overlap_peakLength.bed", header = F,sep="\t")
terra_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_noOverlap_peakLength.bed", header = F,sep="\t")
rloop_noOverlap <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/rloop_noOverlap_peakLength.bed", header = F,sep="\t")

terra_overlap_min4repeats <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_overlap_min4repeats_peakLength.bed", header = F,sep="\t")
terra_overlap_noConsRepeats <- fread("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_overlap_noConsRepeats_peakLength.bed", header = F,sep="\t")


##########################################################################################################
# Statistical test to check if the overlap between TERRA and R-loop peaks (~800 peaks) would be expected to occur by chance
##########################################################################################################
# load rtracklayer to use import() in order to import bed files directly as GRanges objects.
## try if this works with normal bed files using read_csv() or read.table().
## organize filesbed file nameing (e.g. noChrUn) before finalizing.
library('rtracklayer')
terra_all <- import("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/terra_peaks_peakLength_noChrUn.bed", format = 'bed')
rloop_all <- import("/proj/nb_storage/private/terra_rloop_project/results/peaks/Rfiles/peak_length_distribution/rloop_peaks_peakLength_noChrUn.bed", format = 'bed')

# install regioneR package to use its overlapPermTest() for the statistical test.
#BiocManager::install('regioneR', update=FALSE)
library('regioneR')

# count real number of overlaps of regions in one file with regions from the other.
#filter(!('chr_Un' in V1))
numOverlaps(terra_all, rloop_all, count.once=TRUE)
numOverlaps(rloop_all, terra_all, count.once=TRUE)

#numOverlaps(as.data.frame(terra_all), as.data.frame(rloop_all), count.once=TRUE)
#numOverlaps(as.data.frame(rloop_all), as.data.frame(terra_all), count.once=TRUE)
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10",update=FALSE)

pt <- overlapPermTest(A=terra_all, B=rloop_all, ntimes=100, genome="mm10", alternative = 'auto')
pt
plot(pt)


# terra_all <- terra_all %>%
#   filter(!grepl("chrUn",V1, perl=TRUE))%>%
#   select(V1,V2,V3)
# 
# rloop_all <- rloop_all %>%
#   filter(!grepl("chrUn",V1, perl=TRUE))%>%
#   select(V1,V2,V3)


# traceback()

##########################################################################################################
# End of statistical test to check if the overlap between TERRA and R-loop peaks (~800 peaks) would be expected to occur by chance
##########################################################################################################

# add new column to existing data frames to be able to create boxplots.
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

terra_overlap_min4repeats$Sample="Terra"
terra_overlap_min4repeats$Peak_Group="Overlapping-4CR"

terra_overlap_noConsRepeats$Sample="Terra"
terra_overlap_noConsRepeats$Peak_Group="Overlapping-no4CR"

# combine all data frames.
combined_all <- rbind(terra_all,terra_overlap,terra_noOverlap,rloop_all,rloop_overlap,rloop_noOverlap)
combined_terra <- rbind(terra_all,terra_overlap,terra_noOverlap)
combined_rloop <- rbind(rloop_all,rloop_overlap,rloop_noOverlap)
combined_terra_overlap <- rbind(terra_overlap_min4repeats,terra_overlap_noConsRepeats)
  
# define color-blind palettes.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# do some statistics.
# install.packages("car")
library("car")
leveneTest(V5 ~ Peak_Group, data = combined_terra)
oneway.test(V5 ~ Peak_Group, data = combined_terra)
t_test <- pairwise.t.test(combined_terra$V5, combined_terra$Peak_Group, p.adjust.method = "BH", pool.sd = F)
t_test
str(t_test)
##################################################################################################################
# boxplot: peak length distribution per peak group (divided by all_peaks, overlap and noOverlap), 
# and separated by sample (terra vs rloop)
##################################################################################################################

pdf("Fig1_peakLengthDistr_boxplot.pdf", width = 8, height = 4.5)
combined_all %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  # mutate(outlier = V5 > median(V5) + IQR(V5) * 1.5) %>%
  ggplot(.,aes(x = Peak_Group, y = V5, fill = Peak_Group)) + 
  geom_boxplot(outlier.size = 0) + 
  # geom_jitter(aes(x=Peak_Group, y= ifelse(outlier==TRUE,V5,NA)), pch=21) +
  facet_wrap(~factor(Sample, levels = c("Terra","R-loop")), scales="free_x") + 
  #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x", strip.position = "bottom") + 
  xlab("") + ylab("Peak Length (bp)") + 
  scale_y_continuous(limits = c(100,3000), breaks=seq(0,3000, by=500)) + 
  
  theme_bw() + 
  theme(axis.text.x = element_blank(), legend.position = "top") + 
  scale_fill_manual(values=cbPalette, name = "")
dev.off()

##################################################################################################################
# density plot: same as above with limit peaks of length 3000bp or lower.
##################################################################################################################

pdf("Fig1_peakLengthDistr_densityplot.pdf", width = 8, height = 4.5)
combined_all %>% filter(Sample == "Terra") %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  filter(V5 <= 3000) %>%
  # mutate(outlier = V5 > median(V5) + IQR(V5) * 1.5) %>%
  ggplot(.,aes(x = V5, fill = Peak_Group)) + 
  geom_density(outline.type = "upper",  alpha=0.5) + 
  # geom_jitter(aes(x=Peak_Group, y= ifelse(outlier==TRUE,V5,NA)), pch=21) +
  facet_wrap(~factor(Sample, levels = c("Terra","R-loop")), scales="free_x") + 
  #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x", strip.position = "bottom") + 
  xlab("Peak Length (bp)") + ylab("Density") +
  scale_x_continuous(breaks=seq(0,3000, by=200)) + 
  # xlim(0,3000) +
  theme_bw() + 
  theme(legend.position = "top") + 
  scale_fill_manual(values=cbPalette, name = "")
dev.off()
##################################################################################################################
# boxplot: same as above with outliers plotted separately
##################################################################################################################

combined_all %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  # group_by(Peak_Group) %>% 
  # make new columns only containing TERRA peak lengths for each Peak_Group.
  mutate(overlap_lengths_terra=ifelse(Peak_Group=="Overlapping" & Sample=="Terra",V5,NA), 
         noOverlap_lengths_terra=ifelse(Peak_Group=="Non-Overlapping" & Sample=="Terra",V5,NA),
         all_lengths_terra=ifelse(Peak_Group=="All" & Sample=="Terra",V5,NA)) %>% 
         
  # make new columns only containing R-LOOP peak lengths for each Peak_Group.
  mutate(overlap_lengths_rloop=ifelse(Peak_Group=="Overlapping" & Sample=="R-loop",V5,NA), 
         noOverlap_lengths_rloop=ifelse(Peak_Group=="Non-Overlapping" & Sample=="R-loop",V5,NA),
         all_lengths_rloop=ifelse(Peak_Group=="All" & Sample=="R-loop",V5,NA)) %>% 
 
  # make new columns for TERRA outliers split by Peak_Group.
  mutate(outlier_overlap_terra = overlap_lengths_terra > quantile(overlap_lengths_terra, .75, na.rm = T) + 1.50*IQR(overlap_lengths_terra, na.rm = T),
         outlier_noOverlap_terra = noOverlap_lengths_terra > quantile(noOverlap_lengths_terra, .75, na.rm = T) + 1.50*IQR(noOverlap_lengths_terra, na.rm = T),
         outlier_all_terra = all_lengths_terra > quantile(all_lengths_terra, .75, na.rm = T) + 1.50*IQR(all_lengths_terra, na.rm = T)) %>%
  # make new columns for R-LOOP outliers split by Peak_Group
  mutate(outlier_overlap_rloop = overlap_lengths_rloop > quantile(overlap_lengths_rloop, .75, na.rm = T) + 1.50*IQR(overlap_lengths_rloop, na.rm = T),
         outlier_noOverlap_rloop = noOverlap_lengths_rloop > quantile(noOverlap_lengths_rloop, .75, na.rm = T) + 1.50*IQR(noOverlap_lengths_rloop, na.rm = T),
         outlier_all_rloop = all_lengths_rloop > quantile(all_lengths_rloop, .75, na.rm = T) + 1.50*IQR(all_lengths_rloop, na.rm = T)) %>%
  
    ggplot(.,aes(x = Peak_Group, y = V5, fill = Peak_Group)) + 
  geom_boxplot(outlier.size = 0.2, outlier.shape = NA) + 
  
  # plot TERRA outliers
  geom_jitter(aes(x=Peak_Group, y= ifelse(outlier_overlap_terra==TRUE,V5,NA)), pch=21, alpha=0.3) + 
  geom_jitter(aes(x=Peak_Group, y= ifelse(outlier_noOverlap_terra==TRUE,V5,NA)), pch=21, alpha=0.3) + 
  geom_jitter(aes(x=Peak_Group, y= ifelse(outlier_all_terra==TRUE,V5,NA)), pch=21, alpha=0.3) + 
  
  # plot R-loop outliers
  geom_jitter(aes(x=Peak_Group, y= ifelse(outlier_overlap_rloop==TRUE,V5,NA)), pch=21, alpha=0.3) + 
  geom_jitter(aes(x=Peak_Group, y= ifelse(outlier_noOverlap_rloop==TRUE,V5,NA)), pch=21, alpha=0.3) + 
  geom_jitter(aes(x=Peak_Group, y= ifelse(outlier_all_rloop==TRUE,V5,NA)), pch=21, alpha=0.3) + 
  
  facet_wrap(~factor(Sample, levels = c("Terra","R-loop")), scales="free_x") + 
  #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x", strip.position = "bottom") + 
  xlab("") + ylab("Peak Length (bp)") + ylim(0,3000) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), legend.position = "top") + 
  scale_fill_manual(values=cbPalette, name = "")


# combined_all %>% filter(Peak_Group != "All") %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
#   group_by(Peak_Group) %>%
#   mutate(outlier.high = V5 > quantile(V5, .75) + 1.50*IQR(V5),
#          outlier.low = V5 < quantile(V5, .25) - 1.50*IQR(V5)) %>% 
#   mutate(outlier.color = case_when(outlier.high ~ "red",
#                                    outlier.low ~ "steelblue")) %>% 
#   
#   ggplot(.,aes(x = Peak_Group, y = V5, fill = Peak_Group)) + 
#   geom_boxplot(outlier.shape = NA)  + 
#   geom_jitter(color = combined_all$outlier.color, width = .2) + 
#   # geom_point(aes(outlier.high), pch=21, position = position_jitterdodge(), alpha=0.2) +
#   facet_wrap(~factor(Sample, levels = c("Terra","R-loop")), scales="free_x") + 
#   #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x", strip.position = "bottom") + 
#   xlab("") + ylab("Peak Length (bp)") + ylim(0,3000) + 
#   theme_bw() + 
#   theme(axis.text.x = element_blank(), legend.position = "top") + 
#   scale_fill_manual(values=cbPalette, name = "")



# combined_all$overlapping_peakLengths <- with(combined_all, ifelse(Peak_Group=="Overlapping", V5, NA))
# combined_all$non_overlapping_peakLengths <- with(combined_all, ifelse(Peak_Group=="Non-Overlapping", V5, NA))
# 
# 
# pdf("Fig1_peakLengths_byPeaks_boxplot.pdf", width = 8, height = 4.5)
# combined_all %>% filter(Peak_Group != "All") %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
#   # mutate(overlapping_peakLengths = if_else(Peak_Group=="Overlapping", V5, NA))
#   # transform(combined_all, overlapping_peakLengths= ifelse(Peak_Group=="Overlapping", NA, V5) )
#   # mutate(overlapping_peakLengths = case_when(Peak_Group %in% c("Overlapping")~ 'yes') )
#   # mutate(outlier = V5 > median(V5) + IQR(V5) * 1.5) %>% 
#   mutate(outlier.high.overlapping = overlapping_peakLengths > quantile(overlapping_peakLengths, .75, na.rm = TRUE) + 1.50*IQR(overlapping_peakLengths) ) %>% 
#   mutate(outlier.high.non_overlapping = non_overlapping_peakLengths > quantile(non_overlapping_peakLengths, .75, na.rm = TRUE) + 1.50*IQR(non_overlapping_peakLengths)) %>% 
#   
#   # mutate(outlier.high.non_overlapping = Peak_Group=="NOn-Overlapping" & V5 > quantile(V5, .75) + 1.50*IQR(V5)) %>% 
#   
#   ggplot(.,aes(x = Sample, y = V5, fill = Peak_Group)) + 
#   geom_boxplot(outlier.shape = NA) + 
#   # geom_point(aes(y=V5), pch=21, position = position_jitterdodge(), alpha = 0.5) +
#   # geom_point(data = function(x) dplyr::filter_(x, ~ outlier.low), position = 'jitter') + 
#   geom_jitter(data = function(x) filter(x, outlier.high.overlapping ==T | outlier.high.non_overlapping ==T), position = position_jitterdodge(), pch=21, alpha=0.2) +
#   facet_wrap(~factor(Sample, levels = c("Terra","R-loop")), scales="free_x") + 
#   #facet_wrap(~factor(ordered(Sample, levels = rev(sort(unique(Sample))))), scales="free_x", strip.position = "bottom") + 
#   xlab("") + ylab("Peak Length (bp)") + ylim(0,3000) + 
#   theme_bw() + 
#   theme(axis.text.x = element_blank(), legend.position = "top") + 
#   scale_fill_manual(values=cbPalette, name = "")
# dev.off()
# 
# mutate(outlier = V5 > median(V5) + IQR(V5) * 1.5)
# mutate(outlier.high = V5 > quantile(V5, .75) + 1.50*IQR(V5), outlier.low = V5 < quantile(V5, .25) - 1.50*IQR(V5))
# geom_jitter(data = filter(a, outlier.high ==T | outlier.low == T), color = "red", width = .2)
# transform(data_frame, col4= ifelse(col1==col3, col1+col2, col1+col3)
          

##################################################################################################################
# boxplot of peak length distribution showing TERRA overlapping peaks are generally smaller than non-overlapping
##################################################################################################################

pdf("FigS1B_peakLengths_chr_terraOverlapLarge_boxplot.pdf", width = 8, height = 4.5)
combined_all %>% filter(Peak_Group != "All" & Sample == "Terra") %>% 
  filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Peak_Group)) + 
  geom_boxplot(outlier.size = 0) + 
  xlab("") + ylab("Peak Length (bp)") + 
  scale_y_continuous(limits = c(100,3000), breaks=seq(0,3000, by=500)) + 
  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.title = element_blank(), legend.position = "top") + 
  facet_wrap(~Sample) +
  scale_fill_manual(values=cbPalette)
dev.off()


##################################################################################################################
# boxplot of peak length distribution showing R-LOOP overlapping peaks are generally smaller than non-overlapping
##################################################################################################################

pdf("FigS1C_peakLengths_chr_rloopOverlapSmall_boxplot.pdf", width = 8, height = 4.5)
combined_all %>% filter(Peak_Group != "All" & Sample == "R-loop") %>% 
  filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Peak_Group)) + 
  geom_boxplot(outlier.size = 0) + 
  xlab("") + ylab("Peak Length (bp)") + 
  scale_y_continuous(limits = c(100,3000), breaks=seq(0,3000, by=500)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.title = element_blank(), legend.position = "top") + 
  facet_wrap(~Sample) +
  scale_fill_manual(values=cbPalette)
dev.off()


##################################################################################################################
# boxplots of peak length distribution among all_peaks, overlapping_peaks and non-overlapping_peaks 
# expanding upon the above boxplot to compare peak length distr. among individual chromosomes between 
# terra and rloop.
##################################################################################################################

# for Overlapping Peaks.
pdf("Fig1D_peakLengths_chr_overlap_boxplot.pdf", width = 8, height = 4.5)
combined_all %>% filter(Peak_Group == "Overlapping") %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Sample)) + 
  geom_boxplot(outlier.size = 0) + 
  xlab("") + ylab("Peak Length (bp)") + 
  scale_y_continuous(limits = c(100,3000), breaks=seq(0,3000, by=500)) + 
  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.title = element_blank(), legend.position = "top") + 
  facet_wrap(~Peak_Group) + 
  #guides(fill = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values=cbPalette)
dev.off()

# # for All Peaks.
# pdf("peakLengths_chr_allPeaks_boxplot.pdf", width = 8, height = 4.5)
# combined_all %>% filter(Peak_Group == "All") %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>%
#   ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Sample)) +
#   geom_boxplot(position = "dodge") +
#   xlab("") + ylab("Peak Length (bp)") + ylim(0,3000) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "top", legend.title = element_blank()) +
#   facet_wrap(~Peak_Group) +
#   scale_fill_manual(values=cbPalette)
# dev.off()
# 
# # for Non-Overlapping Peaks.
# pdf("peakLengths_chr_noOverlap_boxplot.pdf", width = 8, height = 4.5)
# combined_all %>% filter(Peak_Group =="Non-Overlapping") %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>%
#   ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = Sample)) +
#   geom_boxplot(position = "dodge") +
#   xlab("") + ylab("Peak Length (bp)") + ylim(0,3000) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "top", legend.title = element_blank()) +
#   facet_wrap(~Peak_Group) +
#   scale_fill_manual(values=cbPalette)
# dev.off()





###########################################################################################################################
# peak distribution histogram showing terra overlapping peaks are longer than rloop ones
# i.e. there is a larger amount of long terra peaks, and a larger amount of short rloop peaks
###########################################################################################################################

pdf("Supp1_overlappingPeakLengths_histogram.pdf", width = 8, height = 4.5)
combined_all %>% filter(Peak_Group == "Overlapping") %>% 
  ggplot(.,aes(x=V5, fill=Sample)) + 
  geom_histogram(binwidth=40, position = "dodge") + 
  scale_x_continuous(limits = c(100,1500), breaks=seq(100,1500, by=80)) + 
  xlab("Peak Length (in bins of 40bp)") + ylab("Nr. of Peaks") + 
  facet_wrap(~Peak_Group) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values=cbPalette)
dev.off()

pdf("Supp1_overlappingPeakLengths_freqpoly.pdf", width = 8, height = 4.5)
combined_all %>% filter(Peak_Group == "Overlapping") %>% 
  ggplot(.,aes(x=V5, fill=Sample, color=Sample)) + 
  # geom_histogram(binwidth=40, position = position_dodge2(preserve = "single")) + 
  geom_freqpoly(size=1) + 
  scale_x_continuous(limits = c(100,1500), breaks=seq(100,1500, by=100)) + 
  xlab("Peak Length (bp)") + ylab("Nr. of Peaks") + 
  facet_wrap(~Peak_Group) + 
  theme_bw() + 
  theme(legend.title = element_blank(), legend.position = "top") + 
  scale_color_manual (values=cbPalette)
dev.off()



###########################################################################################################################
# bar chart of peak distribution per chromosome (nr of peaks per chromosome)
###########################################################################################################################

pdf("Fig1B_peakCount_chr_barplot.pdf", width = 8, height = 4.5)
combined_all %>% filter(!grepl("(random|chrUn|chrM)",V1,perl = T)) %>% 
  ggplot(.,aes(y = factor(V1, levels = unique(terra_all$V1)), fill = Sample)) + 
  geom_bar(position=position_dodge()) +
  facet_grid(Peak_Group ~., scales = "free") + 
  xlab("Nr. of Peaks") + ylab("") + 
  coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.title = element_blank(), legend.position = "top") + 
  scale_fill_manual(values=cbPalette)
dev.off()
