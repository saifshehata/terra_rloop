################################################################################
# Important Note:
# Venn diagrams are not the best idea to visualize overlaps/intersections in 
# genomic regions, because some regions from one dataset might intersect with 
# more than one region from anohter, and vice versa, so the exact numbers will 
# not be correctly reflected in the venn diagram. 
# See https://www.biostars.org/p/176316/ for some details.
# Nonetheless, this is a proof of concept script to show hoe to do it.
################################################################################

# set working directory.
setwd("/cluster/khiom/sshehata001/proj/terra_rloop/results/figures/Rplots")

# load packages
library("data.table")
library("dplyr")
library("VennDiagram")
# library("ggvenn")



# read files.
peak_statistics<-fread("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/repeat_stats/peak_statistics.tsv",header = T,sep="\t")

# define color-blind palettes.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# get numbers for venn diagram
nr_terra_peaks <- as.numeric(peak_statistics %>% filter(Sample == "TERRA" & Peak_Group == "All") %>% select(Peak_Count))
nr_rloop_peaks <- as.numeric(peak_statistics %>% filter(Sample == "R-loop" & Peak_Group == "All") %>% select(Peak_Count))

nr_terra_intersect <- as.numeric(peak_statistics %>% filter(Sample == "TERRA" & Peak_Group == "Intersecting") %>% select(Peak_Count))
nr_rloop_intersect <- as.numeric(peak_statistics %>% filter(Sample == "R-loop" & Peak_Group == "Intersecting") %>% select(Peak_Count))

##############################################################################################
# venn diagram
##############################################################################################

pdf("venn_intersecting_peaks", width = 8, height = 4.5)

grid.newpage()
draw.pairwise.venn(area1 = nr_terra_peaks, area2 = nr_rloop_peaks, cross.area = nr_terra_intersect, fill = cbPalette[1:2], lty = 'blank', category = c("TERRA", "R-loop"))  

dev.off()


grid.newpage()
draw.pairwise.venn(area1 = 589, area2 = 199, cross.area = 8, fill = cbPalette[1:2], lty = 'blank', category = c("Intersecting", "TERRA-dependent"))

pdf("venn_intersecting_genes.pdf", width = 3.9, height = 2)
grid.newpage()
draw.pairwise.venn(area1 = 171, area2 = 199, cross.area = 4, fill = cbPalette[1:2], lty = 'blank', category = c("Intersecting", "TERRA-dependent"))  
dev.off()

