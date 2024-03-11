# set working directory.
setwd("/cluster/khiom/sshehata001/proj/terra_rloop/results/figures/Rplots")

library(Gviz)
library(rtracklayer)

# define color-blind palettes.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# read bigwig files
bwTerra <- import.bw("/cluster/khiom/sshehata001/proj/terra_rloop/results/mapped_reads/SRR2062968_pe_sort_rmdup.bw", as = "GRanges")
bwRloop <- import.bw("/cluster/khiom/sshehata001/proj/terra_rloop/results/mapped_reads/SRR2075686_se_sort_rmdup.bw", as = "GRanges")
bwAtrx <- import.bw("/cluster/khiom/sshehata001/proj/terra_rloop/results/mapped_reads/SRR057567_se_sort_rmdup.bw", as = "GRanges")


# create datatrack and specify chromosome
bwTerraDT <- DataTrack(bwTerra, chromosome = "chr2", name = "TERRA")
bwRloopDT <- DataTrack(bwRloop, chromosome = "chr2", name = "R-loop")
bwAtrxDT <- DataTrack(bwAtrx, chromosome = "chr2", name = "ATRX")

# create genome axis track
genomeAxis <- GenomeAxisTrack(name="Genome") 

# create ideogram track
itrack <- IdeogramTrack(genome = "mm9")

# create gene track
library(GenomicFeatures) 
gtf <- "/cluster/khiom/sshehata001/proj/terra_rloop/raw_data/genome/mm10.ncbiRefSeq_filtered.gtf"
txdb_from_gtf <- makeTxDbFromGFF(file = gtf) 
genetrack <- GeneRegionTrack(txdb_from_gtf, name = "Genes") 

# create sequence track
library(BSgenome.Mmusculus.UCSC.mm10)
strack <- SequenceTrack(Mmusculus)

# create motif track
# library(rtracklayer)
telrep_positions <- import("/cluster/khiom/sshehata001/proj/terra_rloop/results/annotate/terra_peaks_motif_positions.bed")
telrep_track <- AnnotationTrack(telrep_positions, name = "Telomere Repeats", stacking = "dense")

# create highlight track
# htrack <- HighlightTrack(trackList = list(bwTerraDT, bwRloopDT, bwAtrxDT), start = c(95423378), width = 5000, chromosome = "chr9")


# plot visualizations
# plotTracks(list(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack, telrep_track), chromosome = "chr9", from = 95389778, to = 95514233, type = "h", labelPos="above", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "brown", title.width = 1.25, rotation.title=0, sizes = c(0.5,0.5,1.5,1.5,1.5,0.75,0.95))
# plotTracks(strack, chromosome = "chr9", from = 95425898, to = 95426034)
pdf("4tR_intergenic.pdf", width = 3.9, height = 2)
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr9")
plotTracks(list(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr9", from = 95385778, to = 95514233, type = "h", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "#999999", scale = 0.5, margin = 1, innerMargin = 1, col = NULL, labelPos = "beside", sizes = c(0.3,0.3,1,1,1,0.4))
dev.off()

# pde7b
pdf("4tR_intron_pde7b.pdf", width = 3.9, height = 2)
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr10")
plotTracks(list(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr10", from = 20243977, to = 20730071, type = "h", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "#999999", scale = 0.5, margin = 1, innerMargin = 1, col = NULL, labelPos = "beside", sizes = c(0.3,0.3,1,1,1,1.2))
dev.off()

# dapk1
pdf("4tR_intron_dapk1.pdf", width = 3.9, height = 2)
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr13")
plotTracks(list(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr13", from = 60560258, to = 60774759, type = "h", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "#999999", scale = 0.5, margin = 1, innerMargin = 1, col = NULL, labelPos = "beside", sizes = c(0.3,0.3,1,1,1,0.6))
dev.off()

# map2k6 chr11:110380000-110529404
pdf("4tR_intron_map2k6.pdf", width = 3.9, height = 2)
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr11")
plotTracks(list(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr11", from = 110380000, to = 110529404, type = "h", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "#999999", scale = 0.5, margin = 1, innerMargin = 1, col = NULL, labelPos = "beside", sizes = c(0.3,0.3,1,1,1,0.2))
dev.off()

# bard1 chr1:71017604-71108481
pdf("4tR_intron_bard1.pdf", width = 3.9, height = 2)
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr1")
plotTracks(list(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr1", from = 71017604, to = 71108481, type = "h", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "#999999", scale = 0.5, margin = 1, innerMargin = 1, col = NULL, labelPos = "beside", sizes = c(0.3,0.3,1,1,1,0.2))
dev.off()

# a1cf chr19:31855850-31992726 with terra-independent 2 rloop peaks in AG-rich regions upstream within same gene
pdf("4tR_intron_a1cf_ag_rich.pdf", width = 3.9, height = 2)
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr19")
plotTracks(list(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr19", from = 31855850, to = 31992726, type = "h", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "#999999", scale = 0.5, margin = 1, innerMargin = 1, col = NULL, labelPos = "beside", sizes = c(0.3,0.3,1,1,1,0.2))
dev.off()

# no4tR intergenic chr9:77615835-77799007
pdf("no4tR_intergenic_1.pdf", width = 3.9, height = 2)
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr9")
plotTracks(list(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr9", from = 77615835, to = 77799007, type = "h", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "#999999", scale = 0.5, margin = 1, innerMargin = 1, col = NULL, labelPos = "beside", sizes = c(0.3,0.3,1,1,1,0.2))
dev.off()

# no4tR intergenic chr19:42104666-42161786 chr19:42069314-42190324
pdf("no4tR_intergenic_2.pdf", width = 3.9, height = 2)
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr19")
plotTracks(list(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr19", from = 42069314, to = 42190324, type = "h", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "#999999", scale = 0.5, margin = 1, innerMargin = 1, col = NULL, labelPos = "beside", sizes = c(0.3,0.3,1,1,1,0.6))
dev.off()

# no4tR intron chr2:73578110-73788220
pdf("no4tR_intron.pdf", width = 3.9, height = 2)
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr2")
plotTracks(list(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr2", from = 73578110, to = 73788220, type = "h", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "#999999", scale = 0.5, margin = 1, innerMargin = 1, col = NULL, labelPos = "beside", sizes = c(0.3,0.3,1,1,1,1.2))
dev.off()

###

# 4tR intergenic tnks chr8:34745000-35200000
pdf("4tR_intron_tnks.pdf", width = 3.9, height = 2)
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr8")
plotTracks(list(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr8", from = 34745000, to = 35200000, type = "h", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "#999999", scale = 0.5, margin = 1, innerMargin = 1, col = NULL, labelPos = "beside", sizes = c(0.3,0.3,1,1,1,0.4))
dev.off()

# 4tR intergenic tnks2 chr19:36677421-36905903
pdf("4tR_intron_tnks2.pdf", width = 3.9, height = 2)
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr19")
plotTracks(list(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr19", from = 36677421, to = 36905903, type = "h", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "#999999", scale = 0.5, margin = 1, innerMargin = 1, col = NULL, labelPos = "beside", sizes = c(0.3,0.3,1,1,1,0.2))
dev.off()



