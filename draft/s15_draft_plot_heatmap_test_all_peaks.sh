#! /bin/bash
#$ -N plot_profile
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 2

# go to directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
mapped_reads=$project_dir/results/mapped_reads
peaks=$project_dir/results/peaks
overlaps=$project_dir/results/overlaps
annotate=$project_dir/results/annotate
meme=$project_dir/results/meme
genes=$project_dir/results/genes
profile=$project_dir/results/profile
figures=$project_dir/results/figures

mkdir -p $profile $figures
cd $profile


######################################################################################################################################################
# plot profile for terra reads on terra OVERLAPPING_4CR (INTERGENIC & INTRON) VS OVERLAPPING_NO4CR (INTERGENIC & INTRON)
######################################################################################################################################################

# define variables. copied from above.
TERRA_AS_BW=$mapped_reads/SRR2062968_pe_sort_rmdup.bw
RLOOP_BW=$mapped_reads/SRR2075686_se_sort_rmdup.bw
ATRX_BW=$mapped_reads/SRR057567_se_sort_rmdup.bw
TERRA_INTERSECT=$overlaps/terra_intersect.narrowPeak
TERRA_INTERSECT_4tR=$overlaps/terra_intersect_4tandem.narrowPeak
TERRA_INTERSECT_no4tR=$overlaps/terra_intersect_no4tandem.narrowPeak

TERRA_ALL_PEAKS=$peaks/terra_sense_peaks_filtered.narrowPeak
RLOOP_ALL_PEAKS=$peaks/rloop_peaks_filtered.narrowPeak
ATRX_ALL_PEAKS=$peaks/atrx_input_peaks.narrowPeak



# terra read coverage on all peak regions
computeMatrix scale-regions --scoreFileName $TERRA_AS_BW \
--regionsFileName $TERRA_ALL_PEAKS \
--startLabel "peak start" \
--endLabel "peak end" \
--samplesLabel "TERRA" \
--upstream 1000 \
--regionBodyLength 1000 \
--downstream 1000 \
--binSize 10 \
--outFileName terra_reads_on_all_peaks.matrix.gz \
--blackListFileName mm10-blacklist.bed \
--numberOfProcessors max

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/terra_reads_on_all_peaks.matrix.gz \
--outFileName $figures/terra_reads_on_all_peaks.pdf \
--outFileSortedRegions terra_reads_on_all_peaks.bed \
--kmeans 3

# atrx read coverage on all peak regions
computeMatrix scale-regions --scoreFileName $ATRX_BW \
--regionsFileName $ATRX_ALL_PEAKS \
--startLabel "peak start" \
--endLabel "peak end" \
--samplesLabel "TERRA" \
--upstream 1000 \
--regionBodyLength 1000 \
--downstream 1000 \
--binSize 10 \
--outFileName atrx_reads_on_all_atrx_peaks.matrix.gz \
--blackListFileName mm10-blacklist.bed \
--numberOfProcessors max

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/atrx_reads_on_all_atrx_peaks.matrix.gz \
--outFileName $figures/atrx_reads_on_all_atrx_peaks.pdf \
--outFileSortedRegions atrx_reads_on_all_atrx_peaks.bed \
--kmeans 3

# atrx read coverage on all terra peak regions
computeMatrix scale-regions --scoreFileName $ATRX_BW \
--regionsFileName $TERRA_ALL_PEAKS \
--startLabel "peak start" \
--endLabel "peak end" \
--samplesLabel "TERRA" \
--upstream 1000 \
--regionBodyLength 1000 \
--downstream 1000 \
--binSize 10 \
--outFileName atrx_reads_on_all_peaks.matrix.gz \
--blackListFileName mm10-blacklist.bed \
--numberOfProcessors max

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/atrx_reads_on_all_peaks.matrix.gz \
--outFileName $figures/atrx_reads_on_all_peaks.pdf \
--outFileSortedRegions atrx_reads_on_all_peaks.bed \
--kmeans 3