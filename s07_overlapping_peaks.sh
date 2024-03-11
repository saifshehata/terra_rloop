#! /bin/bash
#$ -N intersecting_peaks
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 1

# go to directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
peaks=$project_dir/results/peaks
overlaps=$project_dir/results/overlaps
mkdir -p $overlaps
cd $overlaps

# define variables
TERRA_PEAKS=$peaks/terra_sense_peaks.narrowPeak
RLOOP_PEAKS=$peaks/rloop_peaks.narrowPeak
ATRX_PEAKS=$peaks/atrx_input_peaks.narrowPeak


# remove peaks on mitochondrial, random and unknown chromosomes. Will be useful when counting peaks and for the permutation test R script. Also do this for ATRX peaks in case it is needed e.g. for plotProfile or plotHeatmap with deepTools. There are 82 terra peaks, 165 rloop peaks, and 188 atrx peaks that will be removed during this step. 
cat $TERRA_PEAKS | egrep -v '(random|chrUn|chrM)' > $peaks/terra_sense_peaks_filtered.narrowPeak
cat $RLOOP_PEAKS | egrep -v '(random|chrUn|chrM)' > $peaks/rloop_peaks_filtered.narrowPeak
cat $ATRX_PEAKS | egrep -v '(random|chrUn|chrM)' > $peaks/atrx_input_peaks_filtered.narrowPeak
# redefine variables
TERRA_PEAKS=$peaks/terra_sense_peaks_filtered.narrowPeak
RLOOP_PEAKS=$peaks/rloop_peaks_filtered.narrowPeak
ATRX_PEAKS=$peaks/atrx_input_peaks_filtered.narrowPeak

# get intersecting peaks between terra and rloops. -wo keeps original peak coordinates from both peak files and reports the number of overlapping bases.
intersectBed -a $TERRA_PEAKS -b $RLOOP_PEAKS -wo > $overlaps/intersect_terra_rloop_original_peaks.narrowPeak
# intersectBed -a $TERRA_PEAKS -b $RLOOP_PEAKS | cut -f1-3 > $project_dir/results/overlap/intersecting_regions.bed

# extract terra/rloop peaks (with original coordinates) into 2 separate files. sort -u (unique)removes duplicate lines arizing from more than 1 rloop peak overlapping with the same terra peak, resulting in the terra peak being written on more than 1 line, or vice versa. must resort to use -u, otherwize will sort differently. -k1,2 -k2,2n keeps the sorting as it is.
cat $overlaps/intersect_terra_rloop_original_peaks.narrowPeak | awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | sort -k1,1 -k2,2n -u  > $overlaps/terra_intersect.narrowPeak
cat $overlaps/intersect_terra_rloop_original_peaks.narrowPeak | awk 'OFS="\t"{print $11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' | sort -k1,1 -k2,2n -u  > $overlaps/rloop_intersect.narrowPeak

# filter intersecting peaks from total macs2 peaks to get a control non-intersecting peak file (for later terra repeat motif analysis/conparison)
intersectBed -a $TERRA_PEAKS -b $RLOOP_PEAKS -v > $overlaps/terra_no_intersect.narrowPeak
intersectBed -a $RLOOP_PEAKS -b $TERRA_PEAKS -v > $overlaps/rloop_no_intersect.narrowPeak

# # another way to get non-intersecting peaks, but much slower
# grep -v -xf <(cat terra_intersect.narrowPeak) $TERRA_PEAKS > $overlaps/terra_no_intersect2.narrowPeak
# grep -v -xf <(cat rloop_intersect.narrowPeak) $RLOOP_PEAKS > $overlaps/rloop_no_intersect2.narrowPeak

# # get original atrx peaks that intersect with terra peaks. -wa writes the original entry in -a for each overlap. the sort sommand ensures duplicate entries (where 1 peak in -b intersects with >1 peak in -a) are removed
# intersectBed -a $ATRX_PEAKS -b $TERRA_PEAKS -wa  | sort -k1,1 -k2,2n -u > $overlaps/atrx_intersect_terra.narrowPeak

# # do the same but with vice versa
# intersectBed -a $TERRA_PEAKS -b $ATRX_PEAKS -wa  | sort -k1,1 -k2,2n -u > $overlaps/terra_intersect_atrx.narrowPeak
