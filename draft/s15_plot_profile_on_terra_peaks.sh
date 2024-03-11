#! /bin/bash
#$ -N plot_profile_on_terra
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 3

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

TERRA_INTERSECT=$overlaps/terra_intersect.narrowPeak
TERRA_INTERSECT_4tR=$overlaps/terra_intersect_4tandem.narrowPeak
TERRA_INTERSECT_no4tR=$overlaps/terra_intersect_no4tandem.narrowPeak

TERRA_INTERSECT_INTERGENIC=$genes/terra_intersect_intergenic.bed
TERRA_INTERSECT_INTRON=$genes/terra_intersect_intron.bed

TERRA_INTERSECT_4tR_INTERGENIC=$genes/terra_intersect_4tandem_intergenic.bed
TERRA_INTERSECT_4tR_INTRON=$genes/terra_intersect_4tandem_intron.bed
TERRA_INTERSECT_no4tR_INTERGENIC=$genes/terra_intersect_no4tandem_intergenic.bed
TERRA_INTERSECT_no4tR_INTRON=$genes/terra_intersect_no4tandem_intron.bed

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
# terra read coverage on intersecting peak regions split into peaks with or without 4 tandem telomeric repeats
computeMatrix scale-regions --scoreFileName $TERRA_AS_BW \
--regionsFileName $TERRA_INTERSECT_4tR $TERRA_INTERSECT_no4tR \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "TERRA" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName terra_reads_on_terra_intersecting_peaks_4tR_no4tR.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile terra_reads_on_terra_intersecting_peaks_4tR_no4tR.matrix.gz \
--outFileName $figures/terra_reads_on_terra_intersecting_peaks_4tR_no4tR.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/terra_reads_on_terra_intersecting_peaks_4tR_no4tR.matrix.gz \
--outFileName $figures/terra_reads_on_terra_intersecting_peaks_4tR_no4tR_heatmap.pdf


####################### TERRA ######################################
# terra read coverage on intersecting peak regions split into intron and intergenic peak regions
computeMatrix scale-regions --scoreFileName $TERRA_AS_BW \
--regionsFileName $TERRA_INTERSECT_INTERGENIC $TERRA_INTERSECT_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "TERRA" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName terra_reads_on_terra_intersecting_peaks_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile terra_reads_on_terra_intersecting_peaks_intergenic_intron.matrix.gz \
--outFileName $figures/terra_reads_on_terra_intersecting_peaks_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/terra_reads_on_terra_intersecting_peaks_intergenic_intron.matrix.gz \
--outFileName $figures/terra_reads_on_terra_intersecting_peaks_intergenic_intron_heatmap.pdf


####################### TERRA ######################################
# terra read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $TERRA_AS_BW \
--regionsFileName $TERRA_INTERSECT_4tR_INTERGENIC $TERRA_INTERSECT_4tR_INTRON $TERRA_INTERSECT_no4tR_INTERGENIC $TERRA_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "TERRA" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName terra_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile terra_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/terra_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/terra_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/terra_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf


####################### R-loop ######################################
# rloop read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $RLOOP_BW \
--regionsFileName $TERRA_INTERSECT_4tR_INTERGENIC $TERRA_INTERSECT_4tR_INTRON $TERRA_INTERSECT_no4tR_INTERGENIC $TERRA_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "R-loop" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName rloop_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile rloop_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/rloop_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/rloop_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/rloop_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf


####################### ATRX ######################################
# atrx read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $ATRX_BW \
--regionsFileName $TERRA_INTERSECT_4tR_INTERGENIC $TERRA_INTERSECT_4tR_INTRON $TERRA_INTERSECT_no4tR_INTERGENIC $TERRA_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "ATRX" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName atrx_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile atrx_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/atrx_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/atrx_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/atrx_reads_on_terra_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf

