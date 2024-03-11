# set working directory.
setwd("/cluster/khiom/sshehata001/proj/terra_rloop/results/figures/Rplots")

# load packages.
library("data.table")
library("ggplot2")
library("dplyr")

# read bed files with intersecting peak coordinates.
terra_intron <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/genes/terra_intersect_intron.bed", header = F,sep="\t")
terra_intergenic <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/genes/terra_intersect_intergenic.bed", header = F,sep="\t")
terra_no4t_intron <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/genes/terra_intersect_no4tandem_intron.bed", header = F,sep="\t")
terra_no4t_intergenic <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/genes/terra_intersect_no4tandem_intergenic.bed", header = F,sep="\t")
terra_4t_intron <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/genes/terra_intersect_4tandem_intron.bed", header = F,sep="\t")
terra_4t_intergenic <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/genes/terra_intersect_4tandem_intergenic.bed", header = F,sep="\t")

rloop_intron <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/genes/rloop_intersect_intron.bed", header = F,sep="\t")
rloop_intergenic <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/genes/rloop_intersect_intergenic.bed", header = F,sep="\t")
rloop_no4t_intron <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/genes/rloop_intersect_no4tandem_intron.bed", header = F,sep="\t")
rloop_no4t_intergenic <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/genes/rloop_intersect_no4tandem_intergenic.bed", header = F,sep="\t")
rloop_4t_intron <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/genes/rloop_intersect_4tandem_intron.bed", header = F,sep="\t")
rloop_4t_intergenic <- fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/genes/rloop_intersect_4tandem_intergenic.bed", header = F,sep="\t")



# add new column to existing data frames to be able to create plots.
terra_intron$Sample="TERRA"
terra_intron$Repeat_Group="All"
terra_intergenic$Sample="TERRA"
terra_intergenic$Repeat_Group="All"

terra_no4t_intron$Sample="TERRA"
terra_no4t_intron$Repeat_Group="No_4_tandem"
terra_no4t_intergenic$Sample="TERRA"
terra_no4t_intergenic$Repeat_Group="No_4_tandem"

terra_4t_intron$Sample="TERRA"
terra_4t_intron$Repeat_Group="4_tandem"
terra_4t_intergenic$Sample="TERRA"
terra_4t_intergenic$Repeat_Group="4_tandem"

rloop_intron$Sample="R-loop"
rloop_intron$Repeat_Group="All"
rloop_intergenic$Sample="R-loop"
rloop_intergenic$Repeat_Group="All"

rloop_no4t_intron$Sample="R-loop"
rloop_no4t_intron$Repeat_Group="No_4_tandem"
rloop_no4t_intergenic$Sample="R-loop"
rloop_no4t_intergenic$Repeat_Group="No_4_tandem"

rloop_4t_intron$Sample="R-loop"
rloop_4t_intron$Repeat_Group="4_tandem"
rloop_4t_intergenic$Sample="R-loop"
rloop_4t_intergenic$Repeat_Group="4_tandem"

# combine all data frames.
combined_all <- rbind(terra_intron, terra_intergenic, terra_no4t_intron, terra_no4t_intergenic, terra_4t_intron, terra_4t_intergenic, rloop_intron, rloop_intergenic, rloop_no4t_intron, rloop_no4t_intergenic, rloop_4t_intron, rloop_4t_intergenic)

# set column names for combined data frame
colnames(combined_all) <- c("chr", "start", "end", "peak_name", "annotation", "Sample", "Repeat_Group")

# get peak lenghts into a new column
combined_all$Peak_Length <- combined_all$end - combined_all$start

# sort by Sample, then Repeat_Group, then annotation, then chr
# issue with sorting by chr is that it starts: chr1, chr10, chr11..., rather than chr1, chr2...
combined_all <- combined_all[order(Sample, Repeat_Group, annotation, chr),]
# method = c("auto", "shell", "quick", "radix")

# define color-blind palettes.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##################################################################################################################
# boxplots of peak length distribution among all_peaks, overlapping_peaks and non-overlapping_peaks 
# expanding upon the above boxplot to compare peak length distr. among individual chromosomes between 
# terra and rloop.
##################################################################################################################

# for Overlapping Peaks.
pdf("fig3_draft_peak_lengths_per_annotation_and_repeats_boxplot.pdf", width = 8, height = 4.5)
combined_all %>% filter(Sample == "R-loop" & Repeat_Group == "4_tandem") %>% 
  #filter(!grepl("(random|chrUn|chrM)",V1, perl = T)) %>% 
  #ggplot(.,aes(x = factor(V1, levels = unique(terra_all$V1)), y = V5, fill = factor(Sample, levels = unique(Sample)))) + 
  ggplot(.,aes(x = factor(chr, levels = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY')), y = Peak_Length, fill = annotation)) + 
  # levels = unique(combined_all$chr)
  geom_boxplot(outlier.size = 0) + 
  xlab("") + ylab("Peak Length (bp)") + 
  # scale_y_continuous(limits = c(100,3000), breaks=seq(0,3000, by=500)) + 
  
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.title = element_blank(), legend.position = "top", text = element_text(size = 22)) + 
  facet_wrap(~Sample) + 
  #guides(fill = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values=cbPalette)
dev.off()
