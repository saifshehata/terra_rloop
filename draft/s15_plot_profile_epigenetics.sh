#! /bin/bash
#$ -N plot_profile_histones
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

H3K27me3=$mapped_reads/SRR7088942_H3K27me3.bwa_mm10.sorted.rg.rmdup.bw
H3K9me3=$mapped_reads/SRR7088946_H3K9me3.bwa_mm10.sorted.rg.rmdup.bw
H3K4me3=$mapped_reads/SRR7088938_H3K4me3.bwa_mm10.sorted.rg.rmdup.bw

H3K27me3_2=$mapped_reads/SRR944121_H3K27me3.bwa_mm10.sorted.rg.rmdup.bw
H3K4me3_2=$mapped_reads/SRR944122_H3K4me3.bwa_mm10.sorted.rg.rmdup.bw
H3K36me3=$mapped_reads/SRR944123_H3K36me3.bwa_mm10.sorted.rg.rmdup.bw
RNAPII=$mapped_reads/SRR944124_RNA-POLII-Ser5P.bwa_mm10.sorted.rg.rmdup.bw
EZH2=$mapped_reads/SRR944120_Ezh2.bwa_mm10.sorted.rg.rmdup.bw

H3K4me1=$mapped_reads/SRR6159353_H3K4me1.bwa_mm10.sorted.rg.rmdup.bw
H3K4me2=$mapped_reads/SRR6159354_H3K4me2.bwa_mm10.sorted.rg.rmdup.bw
H3K4me3_3=$mapped_reads/SRR6159355_H3K4me3.bwa_mm10.sorted.rg.rmdup.bw
H3K27ac=$mapped_reads/SRR6159358_H3K27ac.bwa_mm10.sorted.rg.rmdup.bw

H3K4me1_2=$mapped_reads/SRR6159363_H3K4me1.bwa_mm10.sorted.rg.rmdup.bw
H3K4me2_2=$mapped_reads/SRR6159372_H3K27ac.bwa_mm10.sorted.rg.rmdup.bw
H3K4me3_3=$mapped_reads/SRR6159377_H3K4me2.bwa_mm10.sorted.rg.rmdup.bw
H3K27ac_2=$mapped_reads/SRR6159378_H3K4me3.bwa_mm10.sorted.rg.rmdup.bw



RLOOP_INTERSECT=$overlaps/rloop_intersect.narrowPeak
RLOOP_INTERSECT_4tR=$overlaps/rloop_intersect_4tandem.narrowPeak
RLOOP_INTERSECT_no4tR=$overlaps/rloop_intersect_no4tandem.narrowPeak

RLOOP_INTERSECT_INTERGENIC=$genes/rloop_intersect_intergenic.bed
RLOOP_INTERSECT_INTRON=$genes/rloop_intersect_intron.bed

RLOOP_INTERSECT_4tR_INTERGENIC=$genes/rloop_intersect_4tandem_intergenic.bed
RLOOP_INTERSECT_4tR_INTRON=$genes/rloop_intersect_4tandem_intron.bed
RLOOP_INTERSECT_no4tR_INTERGENIC=$genes/rloop_intersect_no4tandem_intergenic.bed
RLOOP_INTERSECT_no4tR_INTRON=$genes/rloop_intersect_no4tandem_intron.bed


####################### H3K27me3 ######################################
# H3K27me3 read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $H3K27me3 \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "H3K27me3" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName H3K27me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile H3K27me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/H3K27me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/H3K27me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions H3K27me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/H3K27me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf


####################### H3K9me3 ######################################
# H3K9me3 read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $H3K9me3 \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "H3K9me3" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName H3K9me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile H3K9me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/H3K9me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/H3K9me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions H3K9me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/H3K9me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf



####################### H3K4me3 ######################################
# H3K4me3 read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $H3K4me3 \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "TERRA" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName H3K4me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile H3K4me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/H3K4me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/H3K4me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions H3K4me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/H3K4me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf




# H3K27me3_2=$mapped_reads/SRR944121_H3K27me3.bwa_mm10.sorted.rg.rmdup.bw
# H3K4me3_2=$mapped_reads/SRR944122_H3K4me3.bwa_mm10.sorted.rg.rmdup.bw
# H3K36me3=$mapped_reads/SRR944123_H3K36me3.bwa_mm10.sorted.rg.rmdup.bw

####################### H3K27me3_2 ######################################
# H3K27me3_2 read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $H3K27me3_2 \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "H3K27me3_2" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName H3K27me3_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile H3K27me3_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/H3K27me3_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/EZH2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions H3K27me3_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/H3K27me3_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf


####################### H3K4me3_2 ######################################
# H3K4me3_2 read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $H3K4me3_2 \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "H3K4me3_2" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName H3K4me3_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile H3K4me3_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/H3K4me3_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/H3K4me3_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions H3K4me3_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/H3K4me3_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf


####################### H3K36me3 ######################################
# H3K36me3 read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $H3K36me3 \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "H3K36me3" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName H3K36me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile H3K36me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/H3K36me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/H3K36me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions H3K36me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/H3K36me3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf


# RNAPII=$mapped_reads/SRR944124_RNA-POLII-Ser5P.bwa_mm10.sorted.rg.rmdup.bw
# EZH2=$mapped_reads/SRR944120_Ezh2.bwa_mm10.sorted.rg.rmdup.bw


####################### RNAPII ######################################
# RNAPII read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $RNAPII \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "RNAPII" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName RNAPII_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile RNAPII_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/RNAPII_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/RNAPII_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions RNAPII_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/RNAPII_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf



####################### EZH2 ######################################
# EZH2 read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $EZH2 \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "EZH2" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName EZH2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile EZH2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/EZH2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/EZH2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions EZH2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/EZH2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf


# H3K4me1=$mapped_reads/SRR6159353_H3K4me1.bwa_mm10.sorted.rg.rmdup.bw
# H3K4me2=$mapped_reads/SRR6159354_H3K4me2.bwa_mm10.sorted.rg.rmdup.bw
# H3K4me3_3=$mapped_reads/SRR6159355_H3K4me3.bwa_mm10.sorted.rg.rmdup.bw
# H3K27ac=$mapped_reads/SRR6159358_H3K27ac.bwa_mm10.sorted.rg.rmdup.bw

####################### H3K4me1 ######################################
# H3K4me1 read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $H3K4me1 \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "H3K4me1" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName H3K4me1_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile H3K4me1_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/H3K4me1_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/H3K4me1_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions H3K4me1_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/H3K4me1_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf


####################### H3K4me2 ######################################
# H3K4me2 read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $H3K4me2 \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "H3K4me2" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName H3K4me2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile H3K4me2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/H3K4me2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/H3K4me2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions H3K4me2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/H3K4me2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf


####################### H3K4me3_3 ######################################
# H3K4me3_3 read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $H3K4me3_3 \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "H3K4me3_3" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName H3K4me3_3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile H3K4me3_3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/H3K4me3_3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/H3K4me3_3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions H3K4me3_3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/H3K4me3_3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf


####################### H3K27ac ######################################
# H3K27ac read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $H3K27ac \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "H3K27ac" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName H3K27ac_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile H3K27ac_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/H3K27ac_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/H3K27ac_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions H3K27ac_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/H3K27ac_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf






# H3K4me1_2=$mapped_reads/SRR6159363_H3K4me1.bwa_mm10.sorted.rg.rmdup.bw
# H3K4me2_2=$mapped_reads/SRR6159372_H3K27ac.bwa_mm10.sorted.rg.rmdup.bw
# H3K4me3_3=$mapped_reads/SRR6159377_H3K4me2.bwa_mm10.sorted.rg.rmdup.bw
# H3K27ac_2=$mapped_reads/SRR6159378_H3K4me3.bwa_mm10.sorted.rg.rmdup.bw



####################### H3K4me1_2 ######################################
# H3K4me1_2 read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $H3K4me1_2 \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "H3K4me1_2" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName H3K4me1_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile H3K4me1_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/H3K4me1_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/H3K4me1_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions H3K4me1_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/H3K4me1_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf



####################### H3K4me2_2 ######################################
# H3K4me2_2 read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $H3K4me2_2 \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "H3K4me2_2" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName H3K4me2_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile H3K4me2_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/H3K4me2_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/H3K4me2_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions H3K4me2_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/H3K4me2_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf



####################### H3K4me3_3 ######################################
# H3K4me3_3 read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $H3K4me3_3 \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "H3K4me3_3" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName H3K4me3_3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile H3K4me3_3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/H3K4me3_3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/H3K4me3_3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions H3K4me3_3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/H3K4me3_3_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf



####################### EZH2 ######################################
# H3K27ac_2 read coverage on intersecting peak regions with or without 4 tandem repeats, split into intergenic and intronic peaks
computeMatrix scale-regions --scoreFileName $H3K27ac_2 \
--regionsFileName $RLOOP_INTERSECT_4tR_INTERGENIC $RLOOP_INTERSECT_4tR_INTRON $RLOOP_INTERSECT_no4tR_INTERGENIC $RLOOP_INTERSECT_no4tR_INTRON \
--startLabel "peak start" --endLabel "peak end" --samplesLabel "H3K27ac_2" \
--upstream 2000 --regionBodyLength 3000 --downstream 2000 --binSize 10 \
--outFileName H3K27ac_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--blackListFileName mm10.blacklist_telomeres.bed --numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile H3K27ac_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileName $figures/H3K27ac_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.pdf

# create a heatmap for scores over sets of genomic regions.
plotHeatmap --matrixFile $profile/H3K27ac_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.matrix.gz \
--outFileSortedRegions H3K27ac_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed \
--outFileName $figures/H3K27ac_2_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron_heatmap.pdf



