#! /bin/bash
#$ -N telomeric_repeat_counts
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 1

# go to directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
annotate=$project_dir/results/annotate
repeat_counts=$project_dir/results/repeat_analysis/repeat_counts
mkdir -p $repeat_counts
cd $repeat_counts

################################################################################################ 
# create FILES with repeat counts per peak (not necessarily in tandem) for figure generation with ggplot2 in R
################################################################################################ 


# define variables for peak annotation files stating the number of telomeric repeats per peak.
### annotated original peaks ###
TERRA_PEAKS_ANN_COUNTS=$annotate/terra_peaks_annotate_motif_counts.txt
RLOOP_PEAKS_ANN_COUNTS=$annotate/rloop_peaks_annotate_motif_counts.txt

### annotated intersecting peaks ###
TERRA_INTERSECT_ANN_COUNTS=$annotate/terra_intersect_annotate_motif_counts.txt
RLOOP_INTERSECT_ANN_COUNTS=$annotate/rloop_intersect_annotate_motif_counts.txt

### annotated non-intersecting peaks ###
TERRA_NO_INTERSECT_ANN_COUNTS=$annotate/terra_no_intersect_annotate_motif_counts.txt
RLOOP_NO_INTERSECT_ANN_COUNTS=$annotate/rloop_no_intersect_annotate_motif_counts.txt

# get 2 columns with peak name and number of motifs per peak. Filter out peaks with 0 repeats. i.e. only keep peaks with >= 1 repeat. No need to sort, but looks nice when file is opened in Rstudio or using head
cat $TERRA_PEAKS_ANN_COUNTS | cut -f1,22 | sed "1 d" | awk '$2>=1 {print $0}' | sort -nk2,2 -r > terra_peaks_motif_counts.txt
cat $RLOOP_PEAKS_ANN_COUNTS | cut -f1,22 | sed "1 d" | awk '$2>=1 {print $0}' | sort -nk2,2 -r > rloop_peaks_motif_counts.txt

cat $TERRA_INTERSECT_ANN_COUNTS | cut -f1,22 | sed "1 d" | awk '$2>=1 {print $0}' | sort -nk2,2 -r > terra_intersect_motif_counts.txt
cat $RLOOP_INTERSECT_ANN_COUNTS | cut -f1,22 | sed "1 d" | awk '$2>=1 {print $0}' | sort -nk2,2 -r > rloop_intersect_motif_counts.txt

cat $TERRA_NO_INTERSECT_ANN_COUNTS | cut -f1,22 | sed "1 d" | awk '$2>=1 {print $0}' | sort -nk2,2 -r > terra_no_intersect_motif_counts.txt
cat $RLOOP_NO_INTERSECT_ANN_COUNTS | cut -f1,22 | sed "1 d" | awk '$2>=1 {print $0}' | sort -nk2,2 -r > rloop_no_intersect_motif_counts.txt


# create repeat count distribution plots in R using the files just created.
#r repeat_count_distribution.r

# interesting peak in rloop_no_intersect_narrowPeak_ann_motifs.txt: SRR2075686_rloop_peak_9586. contains several telomeric repeats spaced exactly 40bp apart!
