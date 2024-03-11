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

# terra read coverage on intersecting peak regions woth or without 4 tandem telomeric repeats
computeMatrix scale-regions --scoreFileName $TERRA_AS_BW \
--regionsFileName $TERRA_INTERSECT_4tR $TERRA_INTERSECT_no4tR \
--startLabel "peak start" \
--endLabel "peak end" \
--samplesLabel "TERRA" \
--upstream 1000 \
--regionBodyLength 1000 \
--downstream 1000 \
--binSize 10 \
--outFileName terra_reads_on_intersecting_peaks_with_without_4tandem.matrix.gz \
--blackListFileName mm10-blacklist.bed \
--numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile $profile/terra_reads_on_intersecting_peaks_with_without_4tandem.matrix.gz \
--outFileName $figures/terra_reads_on_intersecting_peaks_with_without_4tandem.pdf \


# rloop read coverage on intersecting peak regions woth or without 4 tandem telomeric repeats
computeMatrix scale-regions --scoreFileName $RLOOP_BW \
--regionsFileName $TERRA_INTERSECT_4tR $TERRA_INTERSECT_no4tR \
--startLabel "peak start" \
--endLabel "peak end" \
--samplesLabel "R-loop" \
--upstream 1000 \
--regionBodyLength 1000 \
--downstream 1000 \
--binSize 10 \
--outFileName rloop_reads_on_intersecting_peaks_with_without_4tandem.matrix.gz \
--blackListFileName mm10-blacklist.bed \
--numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile $profile/rloop_reads_on_intersecting_peaks_with_without_4tandem.matrix.gz \
--outFileName $figures/rloop_reads_on_intersecting_peaks_with_without_4tandem.pdf \


# atrx read coverage on intersecting peak regions woth or without 4 tandem telomeric repeats
computeMatrix scale-regions --scoreFileName $ATRX_BW \
--regionsFileName $TERRA_INTERSECT_4tR $TERRA_INTERSECT_no4tR \
--startLabel "peak start" \
--endLabel "peak end" \
--samplesLabel "ATRX" \
--upstream 1000 \
--regionBodyLength 1000 \
--downstream 1000 \
--binSize 10 \
--outFileName atrx_reads_on_intersecting_peaks_with_without_4tandem.matrix.gz \
--blackListFileName mm10-blacklist.bed \
--numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile $profile/atrx_reads_on_intersecting_peaks_with_without_4tandem.matrix.gz \
--outFileName $figures/atrx_reads_on_intersecting_peaks_with_without_4tandem.pdf \


###
# must plot profile of terra reads on all terra peaks, and plot profile of atrx reads on all atrx reads, and check if locations of highest coverage/cluster_1 in both are within the overlapping peaks. i.e. do the overlapping peaks have highest coverage compared to ALL terra/atrx reads?
###



# # terra read coverage on intersecting peak regions woth or without 4 tandem telomeric repeats
# computeMatrix scale-regions --scoreFileName $TERRA_AS_BW $RLOOP_BW $ATRX_BW \
# --regionsFileName $TERRA_INTERSECT_4tR $TERRA_INTERSECT_no4tR \
# --startLabel "peak start" \
# --endLabel "peak end" \
# --samplesLabel "TERRA" "R-loop" "ATRX" \
# --upstream 1000 \
# --regionBodyLength 1000 \
# --downstream 1000 \
# --binSize 10 \
# --outFileName terra_reads_on_intersecting_peaks_with_without_4tandem_no_blacklist.matrix.gz \
# --numberOfProcessors max

# # create a profile plot for scores over sets of genomic regions.
# plotProfile --matrixFile $profile/terra_reads_on_intersecting_peaks_with_without_4tandem_no_blacklist.matrix.gz \
# --outFileName $figures/terra_reads_on_intersecting_peaks_with_without_4tandem_no_blacklist.pdf \

