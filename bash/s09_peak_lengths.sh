#! /bin/bash
#$ -N peak_lengths
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 1

# go to directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
peaks=$project_dir/results/peaks
overlaps=$project_dir/results/overlaps

peak_lengths=$project_dir/results/repeat_analysis/peak_lengths
mkdir -p $peak_lengths
cd $peak_lengths

################################################################################################ 
# create FILES with peak length of each peak for figure generation with ggplot2 in R
################################################################################################ 

# define variables for files with all peaks, intersecting peaks, and non-intersecting peaks
### original peaks ###
TERRA_PEAKS=$peaks/terra_sense_peaks_filtered.narrowPeak
RLOOP_PEAKS=$peaks/rloop_peaks_filtered.narrowPeak

### intersecting peaks ###
TERRA_INTERSECT=$overlaps/terra_intersect.narrowPeak
RLOOP_INTERSECT=$overlaps/rloop_intersect.narrowPeak

### non-intersecting peaks ###
TERRA_NO_INTERSECT=$overlaps/terra_no_intersect.narrowPeak
RLOOP_NO_INTERSECT=$overlaps/rloop_no_intersect.narrowPeak

# extract chr, start, end, and peakID from MACS2 narrowPeak files, then calculate peak lengths and sort by version to get ch1,chr2,chr3,... instead of chr1,chr10,chr11,...
cat $TERRA_PEAKS | cut -f1-4 | awk 'OFS="\t" {$5=$3-$2}1' | sort -k1 -V  > terra_peaks_peak_lengths.bed
cat $RLOOP_PEAKS | cut -f1-4 | awk 'OFS="\t" {$5=$3-$2}1' | sort -k1 -V  > rloop_peaks_peak_lengths.bed

cat $TERRA_INTERSECT | cut -f1-4 | awk 'OFS="\t" {$5=$3-$2}1' | sort -k1 -V  > terra_intersect_peak_lengths.bed
cat $RLOOP_INTERSECT | cut -f1-4 | awk 'OFS="\t" {$5=$3-$2}1' | sort -k1 -V  > rloop_intersect_peak_lengths.bed

cat $TERRA_NO_INTERSECT | cut -f1-4 | awk 'OFS="\t" {$5=$3-$2}1' | sort -k1 -V  > terra_no_intersect_peak_lengths.bed
cat $RLOOP_NO_INTERSECT | cut -f1-4 | awk 'OFS="\t" {$5=$3-$2}1' | sort -k1 -V  > rloop_no_intersect_peak_lengths.bed


# then run the peak_length_distribution.r script to get box plot of peak sizes per chromosome.

# ### extra datasets after doing some analysis ###

# # define variable for overlapping peaks with/without minimum 4 consecutive repeats.
# TERRA_INTERSECT_4tR=$overlaps/terra_intersect_4tandem.narrowPeak
# TERRA_INTERSECT_NO4tR=$overlaps/terra_intersect_no4tandem.narrowPeak

# # extract chr, start, end, and peakID from MACS2 narrowPeak files, then calculate peak lengths and sort by version to get ch1,chr2,chr3,... instead of chr1,chr10,chr11,...
# cat $TERRA_INTERSECT_4tR | cut -f1-4 | awk 'OFS="\t" {$5=$3-$2}1' | sort -k1 -V  > terra_intersect_4tandem_peak_lengths.bed
# cat $TERRA_INTERSECT_NO4tR | cut -f1-4 | awk 'OFS="\t" {$5=$3-$2}1' | sort -k1 -V  > terra_intersect_no4tandem_peak_lengths.bed
