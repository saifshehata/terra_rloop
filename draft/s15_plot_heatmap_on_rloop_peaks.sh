#! /bin/bash
#$ -N plot_heatmap_on_rloop
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 3

######################################################################################################################################################
# plot heatmap profile for terra, rloop and atrx read coverage on terra intersecting peaks
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

TERRA_ALL_PEAKS=$peaks/terra_sense_peaks_filtered.narrowPeak
RLOOP_ALL_PEAKS=$peaks/rloop_peaks_filtered.narrowPeak
ATRX_ALL_PEAKS=$peaks/atrx_input_peaks_filtered.narrowPeak


# download blacklisted regions
# wget https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz | gunzip
wget --timestamping http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
# unzip but keep original zipped file. This is to avoid repeat wget downloads if the script is rerun
gunzip --stdout mm10.blacklist.bed.gz > mm10.blacklist.bed

# add telomeric regions to blacklisted regions to avoid their very high terra signal obscuring/skewing results from other nont-telomeric peaks. I defined telomeric regions as the last 500kb of each chromosome and the first 3.5M bases at the start of each chromosome, just to avoid any peaks near the chromosome ends. The first 3M bases at the start of each chromosomes are all N's. Chromosome Y has different values added for start and end due to its short length and checking start and end manually on IGV. I removed 150kb from the start, and 1050kb from the end of chrY
cat $profile/mm10.blacklist.bed > $profile/mm10.blacklist_telomeres.bed
egrep -v '(random|chrM|chrUn|chrY)' $genome/mm10.chrom.sizes | awk 'OFS="\t"{print $1,$2-500000,$2}'| sort -V >> $profile/mm10.blacklist_telomeres.bed
echo -e chrY\\t90694698\\t91680079 >> $profile/mm10.blacklist_telomeres.bed
egrep -v '(random|chrM|chrUn|chrY)' $genome/mm10.chrom.sizes | awk 'OFS="\t"{print $1,1,3500000}'| sort -V >> $profile/mm10.blacklist_telomeres.bed
echo -e chrY\\t1\\t150000 >> $profile/mm10.blacklist_telomeres.bed


####################### TERRA ######################################
# terra read coverage on intersecting peak regions
computeMatrix scale-regions --scoreFileName $TERRA_AS_BW --regionsFileName $RLOOP_INTERSECT \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "TERRA" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName terra_reads_on_rloop_intersecting_peaks.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a heatmap for scores over sets of genomic regions. split into clusters using kmeans clustering
plotHeatmap --matrixFile $profile/terra_reads_on_rloop_intersecting_peaks.matrix.gz \
--outFileName $figures/terra_reads_on_rloop_intersecting_peaks_heatmap.pdf \
--outFileSortedRegions terra_reads_on_rloop_intersecting_peaks.bed \
--kmeans 3

####################### TERRA ######################################
# terra read coverage on intersecting peak regions with 4 tandem telomeric repeats
computeMatrix scale-regions --scoreFileName $TERRA_AS_BW --regionsFileName $RLOOP_INTERSECT_4tR \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "TERRA" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName terra_reads_on_rloop_intersecting_peaks_4tandem.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a heatmap for scores over sets of genomic regions. split into clusters using kmeans clustering
plotHeatmap --matrixFile $profile/terra_reads_on_rloop_intersecting_peaks_4tandem.matrix.gz \
--outFileName $figures/terra_reads_on_rloop_intersecting_peaks_4tandem_heatmap.pdf \
--outFileSortedRegions terra_reads_on_rloop_intersecting_peaks_4tandem.bed \
--kmeans 3

####################### TERRA ######################################
# terra read coverage on intersecting peak regions without 4 tandem telomeric repeats
computeMatrix scale-regions --scoreFileName $TERRA_AS_BW --regionsFileName $RLOOP_INTERSECT_no4tR \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "TERRA" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName terra_reads_on_rloop_intersecting_peaks_no4tandem.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/terra_reads_on_rloop_intersecting_peaks_no4tandem.matrix.gz \
--outFileName $figures/terra_reads_on_rloop_intersecting_peaks_no4tandem_heatmap.pdf \
--outFileSortedRegions terra_reads_on_rloop_intersecting_peaks_no4tandem.bed \
--kmeans 3


####################### R-loop ######################################
# rloop read coverage on intersecting peak regions
computeMatrix scale-regions --scoreFileName $RLOOP_BW --regionsFileName $RLOOP_INTERSECT \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "R-loop" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName rloop_reads_on_rloop_intersecting_peaks.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/rloop_reads_on_rloop_intersecting_peaks.matrix.gz \
--outFileName $figures/rloop_reads_on_rloop_intersecting_peaks_heatmap.pdf \
--outFileSortedRegions rloop_reads_on_rloop_intersecting_peaks.bed \
--kmeans 3

####################### R-loop ######################################
# rloop read coverage on intersecting peak regions with 4 tandem telomeric repeats
computeMatrix scale-regions --scoreFileName $RLOOP_BW --regionsFileName $RLOOP_INTERSECT_4tR \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "R-loop" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName rloop_reads_on_rloop_intersecting_peaks_4tandem.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/rloop_reads_on_rloop_intersecting_peaks_4tandem.matrix.gz \
--outFileName $figures/rloop_reads_on_rloop_intersecting_peaks_4tandem_heatmap.pdf \
--outFileSortedRegions rloop_reads_on_rloop_intersecting_peaks_4tandem.bed \
--kmeans 3

####################### R-loop ######################################
# rloop read coverage on intersecting peak regions without 4 tandem telomeric repeats
computeMatrix scale-regions --scoreFileName $RLOOP_BW --regionsFileName $RLOOP_INTERSECT_no4tR \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "R-loop" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName rloop_reads_on_rloop_intersecting_peaks_no4tandem.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/rloop_reads_on_rloop_intersecting_peaks_no4tandem.matrix.gz \
--outFileName $figures/rloop_reads_on_rloop_intersecting_peaks_no4tandem_heatmap.pdf \
--outFileSortedRegions rloop_reads_on_rloop_intersecting_peaks_no4tandem.bed \
--kmeans 3


####################### ATRX ######################################
# atrx read coverage on intersecting peak regions
computeMatrix scale-regions --scoreFileName $ATRX_BW --regionsFileName $RLOOP_INTERSECT \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "ATRX" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName atrx_reads_on_rloop_intersecting_peaks.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/atrx_reads_on_rloop_intersecting_peaks.matrix.gz \
--outFileName $figures/atrx_reads_on_rloop_intersecting_peaks_heatmap.pdf \
--outFileSortedRegions atrx_reads_on_rloop_intersecting_peaks.bed \
--kmeans 3

####################### ATRX ######################################
# atrx read coverage on intersecting peak regions with 4 tandem telomeric repeats
computeMatrix scale-regions --scoreFileName $ATRX_BW --regionsFileName $RLOOP_INTERSECT_4tR \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "ATRX" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName atrx_reads_on_rloop_intersecting_peaks_4tandem.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/atrx_reads_on_rloop_intersecting_peaks_4tandem.matrix.gz \
--outFileName $figures/atrx_reads_on_rloop_intersecting_peaks_4tandem_heatmap.pdf \
--outFileSortedRegions atrx_reads_on_rloop_intersecting_peaks_4tandem.bed \
--kmeans 3

####################### ATRX ######################################
# atrx read coverage on intersecting peak regions without 4 tandem telomeric repeats
computeMatrix scale-regions --scoreFileName $ATRX_BW --regionsFileName $RLOOP_INTERSECT_no4tR \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "ATRX" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName atrx_reads_on_rloop_intersecting_peaks_no4tandem.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/atrx_reads_on_rloop_intersecting_peaks_no4tandem.matrix.gz \
--outFileName $figures/atrx_reads_on_rloop_intersecting_peaks_no4tandem_heatmap.pdf \
--outFileSortedRegions atrx_reads_on_rloop_intersecting_peaks_no4tandem.bed \
--kmeans 3

