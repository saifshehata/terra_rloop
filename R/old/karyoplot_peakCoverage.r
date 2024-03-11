# set working directory.
setwd("/cluster/khiom/sshehata001/proj/terra_rloop/results/figures/Rplots")

#BiocManager::install("ChIPseeker")
#BiocManager::install("Gviz")
#BiocManager::install("trackViewer")
#BiocManager::install("profileplyr", update = F) # did not work. non-zero exit status
#install.packages("tiff", dependencies = T) # needed for profileplyr, but also did not work

library(ChIPseeker)

# load peaks
rloop_intersect <- readPeakFile("/cluster/khiom/sshehata001/proj/terra_rloop/results/overlaps/rloop_intersect.narrowPeak")
terra_intersect <- readPeakFile("/cluster/khiom/sshehata001/proj/terra_rloop/results/overlaps/terra_intersect.narrowPeak")
#peaks <- readPeakFile("/cluster/khiom/sshehata001/proj/terra_rloop/results/peaks/terra_sense_peaks_filtered.narrowPeak")
rloop_intersect_filtered <- readPeakFile("/cluster/khiom/sshehata001/proj/terra_rloop/results/overlaps/rloop_intersect_10fold_enrichment.narrowPeak")
#terra_peaks_filtered <- readPeakFile("/cluster/khiom/sshehata001/proj/terra_rloop/results/peaks/terra_sense_peaks_filtered_10fold_enrichment.narrowPeak")

#peaks
# plot coverage over all peaks
pdf("rloop_intersect_peak_coverage_on_chromosomes.pdf", width = 8, height = 4.5)
covplot(rloop_intersect, weightCol = "V7", ylab = "Fold Enrichment", title = "Intersecting peaks over chromosomes")
dev.off()

covplot(rloop_intersect_filtered, weightCol = "V7", ylab = "Fold Enrichment")

################################################ Testing Gviz

library(Gviz)

# read bam files
bamTerra <- AlignmentsTrack("/cluster/khiom/sshehata001/proj/terra_rloop/results/mapped_reads/SRR2062968_pe_sort_rmdup.bam", name = "TERRA", isPaired = T)
bamRloop <- AlignmentsTrack("/cluster/khiom/sshehata001/proj/terra_rloop/results/mapped_reads/SRR2075686_se_sort_rmdup.bam", name = "R-loop", isPaired = F)
bamAtrx <- AlignmentsTrack("/cluster/khiom/sshehata001/proj/terra_rloop/results/mapped_reads/SRR057567_se_sort_rmdup.bam", name = "ATRX", isPaired = F)

# create genome axis track
genomeAxis <- GenomeAxisTrack(name="MyAxis") 
genomeAxis

# plot coverage
plotTracks(c(bamTerra, bamRloop, bamAtrx, genomeAxis), chromosome = "chr2", from = 57628851, to = 57630549, type = "coverage", add53=T,add35=T, littleTicks = TRUE, labelPos="below", scale=1500)


library(rtracklayer)
# read bigwig files
bwTerra <- import.bw("/cluster/khiom/sshehata001/proj/terra_rloop/results/mapped_reads/SRR2062968_pe_sort_rmdup.bw", as = "GRanges")
bwRloop <- import.bw("/cluster/khiom/sshehata001/proj/terra_rloop/results/mapped_reads/SRR2075686_se_sort_rmdup.bw", as = "GRanges")
bwAtrx <- import.bw("/cluster/khiom/sshehata001/proj/terra_rloop/results/mapped_reads/SRR057567_se_sort_rmdup.bw", as = "GRanges")


# create datatrack and specify chromosome
bwTerraDT <- DataTrack(bwTerra, chromosome = "chr2", name = "TERRA")
bwRloopDT <- DataTrack(bwRloop, chromosome = "chr2", name = "R-loop")
bwAtrxDT <- DataTrack(bwAtrx, chromosome = "chr2", name = "ATRX")

# create genome axis track
genomeAxis <- GenomeAxisTrack(name="MyAxis") 

# clearSessionCache()
# data(ucscItems)
# ucscTables("mm10", "knownGene")
# track <- UcscTrack(genome = "mm10", chromosome = "chr2", track = "knownGene")
# ideogram track
#Gviz.ucscUrl="https://genome.ucsc.edu/cgi-bin/"

options(Gviz.ucscUrl="http://genome-euro.ucsc.edu/cgi-bin/")
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr2")


# BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene')
# BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')


# gtf <-fread("/cluster/khiom/sshehata001/proj/terra_rloop/raw_data/genome/mm10.refGene.gtf")
# gtf <- GTFFile("/cluster/khiom/sshehata001/proj/terra_rloop/raw_data/genome/mm10.refGene.gtf")
# gtf <- GTFFile('TxDb.Mmusculus.UCSC.mm10.knownGene')
# GCA_000001635.5_GRCm38.p3_full_analysis_set.refseq_annotation.gff.gz
# mm10.refGene.gtf
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# txdb <- loadDb('TxDb.Mmusculus.UCSC.mm10.knownGene')
# data()
# data("geneModels")
# head()
# gtrack <- GeneRegionTrack("TxDb.Mmusculus.UCSC.mm10.knownGene", genome = "mm10", chromosome = "chr13", name = "Gene Model")
# AnnotationTrack(TxDb.Mmusculus.UCSC.mm10.knownGene)


library(GenomicFeatures) 
ensembleGTF <- "/cluster/khiom/sshehata001/proj/terra_rloop/raw_data/genome/mm10.ncbiRefSeq_filtered.gtf"
# ensembleGTF <- "GCA_000001635.5_GRCm38.p3_full_analysis_set.refseq_annotation.gff.gz"
txdbFromGFF <- makeTxDbFromGFF(file = ensembleGTF) 
customFromTxDb <- GeneRegionTrack(txdbFromGFF, name = "Genes") 
plotTracks(customFromTxDb, 
           from=60600607,to=60763901, 
           transcriptAnnotation="gene")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
gene_track <- GeneRegionTrack(txdb, name = "Genes")

# bgrTrack <- BiomartGeneRegionTrack(genome="mm10",
#                                    start=60675338,
#                                    end=60678336,
#                                    chromosome = "chr13",
#                                    name="ENSEMBL",              
#                                    filter=list(source="ensembl_havana"))
# listDatasets()

library(BSgenome.Mmusculus.UCSC.mm10)
strack <- SequenceTrack(Mmusculus, chromosome = 13, )
library(rtracklayer)
telrep_positions <- import("/cluster/khiom/sshehata001/proj/terra_rloop/results/annotate/terra_peaks_motif_positions.bed")
telrep_track <- AnnotationTrack(telrep_positions, name = "Telomere Repeats", stacking = "dense")

# define color-blind palettes.
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

###
# intergenic 4tR
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr9")
plotTracks(c(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, customFromTxDb, telrep_track), chromosome = "chr9", from = 95389778, to = 95514233, type = "h", labelPos="above", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "brown", title.width = 1.25, rotation.title=0, sizes = c(0.5,0.5,1.5,1.5,1.5,0.75,0.95))

###
# intron 4tR chr10:20243977-20730071 pde7b
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr10")
plotTracks(c(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, customFromTxDb), chromosome = "chr10", from = 20243977, to = 20730071, type = "h", labelPos="above", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "brown")
# intron 4tR
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr5")
plotTracks(c(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr5", from = 147343217, to = 147976143, type = "h", labelPos="above", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "brown")
# 4tR intron  with 2 rloop peaks in AG-rich regions upstream within same gene
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr19")
plotTracks(c(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr19", from = 31855850, to = 31992726, type = "h", labelPos="above", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "brown")
###
# 4tR intron dapk1
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr13")
plotTracks(c(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr13", from = 60560258, to = 60774759, type = "h", labelPos="above", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "brown")

# no4tR intergenic chr12:69305436-69417134
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr12")
plotTracks(c(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr12", from = 69305436, to = 69417134, type = "h", labelPos="above", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "brown")
###
# no4tR intergenic chr9:77632835-77799007
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr9")
plotTracks(c(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr9", from = 77632835, to = 77799007, type = "h", labelPos="above", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "brown")
# no4tR intergenic chr19:42104666-42161786
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr19")
plotTracks(c(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr19", from = 42104666, to = 42161786, type = "h", labelPos="above", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "brown")

###
# no4tR intron chr2:73578110-73788220
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr2")
plotTracks(c(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr2", from = 73578110, to = 73788220, type = "h", labelPos="above", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "brown")

# very interesting Rn4.5s oscilating peaks region: chr6:47506921-47834912
itrack <- IdeogramTrack(genome = "mm9", chromosome = "chr6")
plotTracks(c(itrack, genomeAxis, bwTerraDT, bwRloopDT, bwAtrxDT, genetrack), chromosome = "chr6", from = 47506921, to = 47834912, type = "h", labelPos="above", transcriptAnnotation = "gene", background.panel = "#FFFEDB", background.title = "brown")


# programmed cell death genes with intersecting peaks min4t repeats
# ENSMUSG00000004637
# ENSMUSG00000021559
# ENSMUSG00000022150
# ENSMUSG00000020623
# ENSMUSG00000028284
# ENSMUSG00000024500
# ENSMUSG00000044468
# cat terra_intersect_4tandem_annotate_motif_counts.txt | egrep '(ENSMUSG00000004637|ENSMUSG00000021559|ENSMUSG00000022150|ENSMUSG00000020623|ENSMUSG00000028284|ENSMUSG00000024500|ENSMUSG00000044468)' | cut -f16

# dna damage
# ENSMUSG00000026196
# ENSMUSG00000040359
# ENSMUSG00000027242
# ENSMUSG00000046295
# ENSMUSG00000020326
# cat terra_intersect_4tandem_annotate_motif_counts.txt | egrep '(ENSMUSG00000004637|ENSMUSG00000021559|ENSMUSG00000022150|ENSMUSG00000020623|ENSMUSG00000028284|ENSMUSG00000024500|ENSMUSG00000044468)' | cut -f16

# telomere regulation
# ENSMUSG00000031529
# cat terra_intersect_4tandem_annotate_motif_counts.txt | grep ENSMUSG00000031529 | cut -f16
chr13:60560258-60774759
chr13:60600607-60763901
chr13:60675338-60678336
chr13:60676532-60676748
chr13:60676588-60676687
chr8 34826000 35148000
chr8 35134877 35135032
chr2 57628851 57630549


################################## heatmap from deeptools
library(profileplyr)
# import deeptools matrix
proplyr_object <- import_deepToolsMat(con = "/cluster/khiom/sshehata001/proj/terra_rloop/results/profile/terra_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz")














###
afrom=2960000
ato=3160000
#bam file
alTrack <- AlignmentsTrack(system.file(package = "Gviz", "extdata", "gapped.bam"), isPaired = TRUE)
bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr12",
                              start = afrom, end = ato, filter = list(with_ox_refseq_mrna = TRUE),
                              stacking = "dense")
plotTracks(c(bmt, alTrack), from = afrom, to = ato, chromosome = "chr12")
###


plotTracks(c(alTrack, bmt), from = afrom, to = ato, chromosome = "chr12", type = "coverage")


library(trackViewer)
library(profileplyr)


library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

