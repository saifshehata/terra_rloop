#! /bin/bash
#$ -N plot_profile_on_rloop
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 2

######################################################################################################################################################
# plot profile for terra, rloop and atrx reads on terra intersecting peaks with or without 4 tandem repeats, each split into intergenic and intron peaks
# note: peaks intersecting with blacklisted genomic regions known for high non-specific signals in chip-seq experiments are filtered out
######################################################################################################################################################

# go to directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
mapped_reads=$project_dir/results/mapped_reads
peaks=$project_dir/results/peaks
overlaps=$project_dir/results/overlaps
annotate=$project_dir/results/annotate
meme=$project_dir/results/meme
genes=$project_dir/results/genes
genome=$project_dir/raw_data/genome
profile=$project_dir/results/profile
figures=$project_dir/results/figures

mkdir -p $profile $figures
cd $profile

# define variables
TERRA_AS_BW=$mapped_reads/SRR2062968_pe_sort_rmdup.bw
RLOOP_BW=$mapped_reads/SRR2075686_se_sort_rmdup.bw
ATRX_BW=$mapped_reads/SRR057567_se_sort_rmdup.bw

RLOOP_INTERSECT=$overlaps/rloop_intersect.narrowPeak
RLOOP_INTERSECT_4tR=$overlaps/rloop_intersect_4tandem.narrowPeak
RLOOP_INTERSECT_no4tR=$overlaps/rloop_intersect_no4tandem.narrowPeak

RLOOP_INTERSECT_INTERGENIC=$genes/rloop_intersect_intergenic.bed
RLOOP_INTERSECT_INTRON=$genes/rloop_intersect_intron.bed

RLOOP_INTERSECT_4tR_INTERGENIC=$genes/rloop_intersect_4tandem_intergenic.bed
RLOOP_INTERSECT_4tR_INTRON=$genes/rloop_intersect_4tandem_intron.bed
RLOOP_INTERSECT_no4tR_INTERGENIC=$genes/rloop_intersect_no4tandem_intergenic.bed
RLOOP_INTERSECT_no4tR_INTRON=$genes/rloop_intersect_no4tandem_intron.bed

# ATRX_INTERSECT_TERRA=$overlaps/atrx_intersect_terra.narrowPeak
# TERRA_INTERSECT_ATRX=$overlaps/terra_intersect_atrx.narrowPeak

# download blacklisted regions
# wget https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz | gunzip
wget --timestamping http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
# unzip but keep original zipped file. This is to avoid repeat wget downloads if the script is rerun
gunzip --stdout mm10.blacklist.bed.gz > mm10.blacklist.bed


# add subtelomeric regions to blacklisted regions to avoid their very high terra signal skewing results from other non-telomeric peaks.
# I defined subtelomeric regions as the last 500kb of each chromosome (referenced in text) to avoid any peaks near the chromosome ends and focus on non-tlomeric loci.
# The first 3M bases at the start of each chromosomes are all N's, so 3.5M bases (3M + 500kb) at the start of each chromosome were removed.
# Chromosome Y has different values added for start and end due to its short length and checking start and end manually on IGV. I removed 150kb from the start, and 1050kb from the end of chrY.
cat $profile/mm10.blacklist.bed > $profile/mm10.blacklist_telomeres.bed
egrep -v '(random|chrM|chrUn|chrY)' $genome/mm10.chrom.sizes | awk 'OFS="\t"{print $1,$2-500000,$2}'| sort -V >> $profile/mm10.blacklist_telomeres.bed
echo -e chrY\\t90694698\\t91680079 >> $profile/mm10.blacklist_telomeres.bed
egrep -v '(random|chrM|chrUn|chrY)' $genome/mm10.chrom.sizes | awk 'OFS="\t"{print $1,1,3500000}'| sort -V >> $profile/mm10.blacklist_telomeres.bed
echo -e chrY\\t1\\t150000 >> $profile/mm10.blacklist_telomeres.bed

# Just for clarity, there were very few overlaps with blacklisted regions (only 13 overlaps between shared TERRA peaks lacking 4 tandem repeats, and no overlaps at all in shared TERRA peaks with 4 tandem repeats).
# There were 20 and 24 peaks in shared TERRA peaks with and without 4 tandem repeats, respectively, that overlaped with subtelomeric regions.
# The main takeaways from the figures did not change with or without removal of blacklist/subtelomere regions, but the signal coming from non-telomeric peaks was clearer when they were removed.
# Removal was achieved by adding the mm10.blacklist_telomeres.bed file to the computeMatrix command below using the '--blackListFileName' parameter.


####################### TERRA ######################################
# terra read coverage on intersecting peak regions split into peaks with or without 4 tandem telomeric repeats
computeMatrix scale-regions --scoreFileName $TERRA_AS_BW \
--regionsFileName $RLOOP_INTERSECT_4tR $RLOOP_INTERSECT_no4tR \
--samplesLabel "TERRA" \
--upstream 3000 --regionBodyLength 3000 --downstream 3000 --binSize 10 \
--outFileName terra_reads_on_rloop_intersecting_peaks_4tR_no4tR.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile terra_reads_on_rloop_intersecting_peaks_4tR_no4tR.matrix.gz \
--startLabel "Peak Start" --endLabel "Peak End" \
--outFileName $figures/terra_reads_on_rloop_intersecting_peaks_4tR_no4tR.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/terra_reads_on_rloop_intersecting_peaks_4tR_no4tR.matrix.gz \
--outFileSortedRegions terra_reads_on_rloop_intersecting_peaks_4tR_no4tR.bed \
--outFileName $figures/terra_reads_on_rloop_intersecting_peaks_4tR_no4tR_heatmap.pdf \
--startLabel "Peak
Start" --endLabel "Peak
End"  --whatToShow "plot and heatmap" --heatmapHeight 15


####################### TERRA ######################################
# terra read coverage on intersecting peak regions split into intron and intergenic peak regions
computeMatrix scale-regions --scoreFileName $TERRA_AS_BW \
--regionsFileName $RLOOP_INTERSECT_INTERGENIC $RLOOP_INTERSECT_INTRON \
--samplesLabel "TERRA" \
--upstream 3000 --regionBodyLength 3000 --downstream 3000 --binSize 10 \
--outFileName terra_reads_on_rloop_intersecting_peaks_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile terra_reads_on_rloop_intersecting_peaks_intergenic_intron.matrix.gz \
--startLabel "Peak Start" --endLabel "Peak End" \
--outFileName $figures/terra_reads_on_rloop_intersecting_peaks_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/terra_reads_on_rloop_intersecting_peaks_intergenic_intron.matrix.gz \
--outFileSortedRegions terra_reads_on_rloop_intersecting_peaks_intergenic_intron.bed \
--outFileName $figures/terra_reads_on_rloop_intersecting_peaks_intergenic_intron_heatmap.pdf \
--startLabel "Peak
Start" --endLabel "Peak
End"  --whatToShow "plot and heatmap" --heatmapHeight 15

####################### TERRA ######################################
# terra read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $TERRA_AS_BW \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--samplesLabel "TERRA" \
--upstream 3000 --regionBodyLength 3000 --downstream 3000 --binSize 10 \
--outFileName terra_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile terra_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--startLabel "Peak Start" --endLabel "Peak End" \
--outFileName $figures/terra_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/terra_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions terra_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/terra_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf \
--startLabel "Peak
Start" --endLabel "Peak
End"  --whatToShow "plot and heatmap" --heatmapHeight 15

####################### R-loop ######################################
# rloop read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $RLOOP_BW \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--samplesLabel "R-loop" \
--upstream 3000 --regionBodyLength 3000 --downstream 3000 --binSize 10 \
--outFileName rloop_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile rloop_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--startLabel "Peak Start" --endLabel "Peak End" \
--outFileName $figures/rloop_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/rloop_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions rloop_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/rloop_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf \
--startLabel "Peak
Start" --endLabel "Peak
End"  --whatToShow "plot and heatmap" --heatmapHeight 15


####################### ATRX ######################################
# atrx read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $ATRX_BW \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--samplesLabel "ATRX" \
--upstream 3000 --regionBodyLength 3000 --downstream 3000 --binSize 10 \
--outFileName atrx_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile atrx_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--startLabel "Peak Start" --endLabel "Peak End" \
--outFileName $figures/atrx_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/atrx_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions atrx_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/atrx_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf \
--startLabel "Peak
Start" --endLabel "Peak
End"  --whatToShow "plot and heatmap" --heatmapHeight 15




################################### extra heatmaps for clearer figures for publication ################################

# terra read coverage on intersecting peak regions with 4 tandem repeats split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $TERRA_AS_BW \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON \
--samplesLabel "TERRA" \
--upstream 3000 --regionBodyLength 3000 --downstream 3000 --binSize 10 \
--outFileName terra_reads_on_rloop_intersecting_peaks_4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile terra_reads_on_rloop_intersecting_peaks_4tR_intergenic_intron.matrix.gz \
--startLabel "Peak Start" --endLabel "Peak End" \
--outFileName $figures/terra_reads_on_rloop_intersecting_peaks_4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/terra_reads_on_rloop_intersecting_peaks_4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions terra_reads_on_rloop_intersecting_peaks_4tR_intergenic_intron.bed \
--outFileName $figures/terra_reads_on_rloop_intersecting_peaks_4tR_intergenic_intron_heatmap.pdf \
--startLabel "Peak
Start" --endLabel "Peak
End"  --whatToShow "heatmap only" --heatmapHeight 15


# terra read coverage on intersecting peak regions without 4 tandem repeats split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $TERRA_AS_BW \
--regionsFileName $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--samplesLabel "TERRA" \
--upstream 3000 --regionBodyLength 3000 --downstream 3000 --binSize 10 \
--outFileName terra_reads_on_rloop_intersecting_peaks_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile terra_reads_on_rloop_intersecting_peaks_no4tR_intergenic_intron.matrix.gz \
--startLabel "Peak Start" --endLabel "Peak End" \
--outFileName $figures/terra_reads_on_rloop_intersecting_peaks_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/terra_reads_on_rloop_intersecting_peaks_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions terra_reads_on_rloop_intersecting_peaks_no4tR_intergenic_intron.bed \
--outFileName $figures/terra_reads_on_rloop_intersecting_peaks_no4tR_intergenic_intron_heatmap.pdf \
--startLabel "Peak
Start" --endLabel "Peak
End"  --whatToShow "heatmap only" --heatmapHeight 15


### ATRX ### 

# atrx read coverage on intersecting peak regions with 4 tandem repeats split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $ATRX_BW \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON \
--samplesLabel "ATRX" \
--upstream 3000 --regionBodyLength 3000 --downstream 3000 --binSize 10 \
--outFileName atrx_reads_on_rloop_intersecting_peaks_4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile atrx_reads_on_rloop_intersecting_peaks_4tR_intergenic_intron.matrix.gz \
--startLabel "Peak Start" --endLabel "Peak End" \
--outFileName $figures/atrx_reads_on_rloop_intersecting_peaks_4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/atrx_reads_on_rloop_intersecting_peaks_4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions atrx_reads_on_rloop_intersecting_peaks_4tR_intergenic_intron.bed \
--outFileName $figures/atrx_reads_on_rloop_intersecting_peaks_4tR_intergenic_intron_heatmap.pdf \
--startLabel "Peak
Start" --endLabel "Peak
End"  --whatToShow "heatmap only" --heatmapHeight 15


# atrx read coverage on intersecting peak regions without 4 tandem repeats split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $ATRX_BW \
--regionsFileName $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--samplesLabel "ATRX" \
--upstream 3000 --regionBodyLength 3000 --downstream 3000 --binSize 10 \
--outFileName atrx_reads_on_rloop_intersecting_peaks_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile atrx_reads_on_rloop_intersecting_peaks_no4tR_intergenic_intron.matrix.gz \
--startLabel "Peak Start" --endLabel "Peak End" \
--outFileName $figures/atrx_reads_on_rloop_intersecting_peaks_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/atrx_reads_on_rloop_intersecting_peaks_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions atrx_reads_on_rloop_intersecting_peaks_no4tR_intergenic_intron.bed \
--outFileName $figures/atrx_reads_on_rloop_intersecting_peaks_no4tR_intergenic_intron_heatmap.pdf \
--startLabel "Peak
Start" --endLabel "Peak
End"  --whatToShow "heatmap only" --heatmapHeight 15
