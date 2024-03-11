### check out this peak in IGV, as it has interesting TTAGGG pattern of 6 11 6 spaced apart
(ngs) [sshehata001@login2 annotate]$ cat terra_intersect_annotate_all_repeats.txt | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | head -2
terra_sense_peak_27704 23 137 30 8 78 52 7 225 72 7 16 99 375 115 50 128 58 188 52 150 102 931 6 11 6 6 5 6 6 5 6 5 6 5 6 12 5 6 9 5 6 9 5 6 6 9 5 6 9 5 6 6 6 9 5 6 9 5 6 9 6 6 11 6 9 6 10 6 6 73 187 48 11 6 136 12 113 80 11 6 76 6 11 6 47 306 345 54 87 90 6 134 12 57 6 1273 108 135 29 59 115 42 1 9 1 145 46 13 1 19 232 1 136 13 1 19 42 1 31 108 13 41 84 65 21 1 10 246 29 267 76 160 112 483 6 11 6 79 6 11 6 63 6 11 6 77 6 11 6 77 6 11 6 79 6 11 6 69 6 11 6 47 306 345 54 87 90 6 134 12 57 6 1273 108 135 29 59 115 42 1 9 1 145 46 13 1 19 232 1 136 13 1 19 42 1 31 108 13 41 84 65 21 1 10 246 29 267 76 160 112 563 45 36 3228 178 92 156

###################################################################################################
# TESTING
###################################################################################################

# load some modules
module load bioinfo-tools
module load BEDOPS
module load BEDTools
module load HOMER

# go to directory.
cd $project_dir/results/peaks/overlap/

# define variables for terra and rloop narrowPeak files.
TERRA_PEAKS=$project_dir/results/peaks/terra/terra_sense_bampe_q0.01_peaks.narrowPeak
RLOOP_PEAKS=$project_dir/results/peaks/rloop/rloop_q0.01_peaks.narrowPeak

# find overlapping peaks between terra and rloops. -wo keeps original peak coordinates from both peak files and reports the number of overlapping bases.
bedtools intersect -a $TERRA_PEAKS -b $RLOOP_PEAKS -wo > $project_dir/results/peaks/overlap/terraAS-S_rloop_overlap_originalPeaks.bed

# extract terra/rloop peaks (with original coordinates) into 2 separate files. sort -u removes duplicate lines arizing from more than 1 rloop peak overlapping with the same terra peak, resulting in the terra peak being written on more than 1 line, or vice versa. must resort to use -u, otherwize will sort differently. -k1,2 -k2,2n keeps the sorting as it is.
cat terraAS-S_rloop_overlap_originalPeaks.bed | awk 'OFS="\t"{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | sort -k1,1 -k2,2n -u  > terra_overlap_terraAS-S.narrowPeak
cat terraAS-S_rloop_overlap_originalPeaks.bed | awk 'OFS="\t"{print $11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' | sort -k1,1 -k2,2n -u  > rloop_overlap_terraAS-S.narrowPeak


# annotate peaks.
cd $project_dir/results/peaks/annotation/

TERRA_OVERLAP=$project_dir/results/peaks/overlap/terra_overlap_terraAS-S.narrowPeak
RLOOP_OVERLAP=$project_dir/results/peaks/overlap/rloop_overlap_terraAS-S.narrowPeak

annotatePeaks.pl $TERRA_OVERLAP mm10 -size given -m telomeric.repeat -cpu 10 > terraAS-S_overlap_narrowPeak_ann_allRepeats.txt
annotatePeaks.pl $RLOOP_OVERLAP mm10 -size given -m telomeric.repeat -cpu 10 > rloop_terraAS-S_overlap_narrowPeak_ann_allRepeats.txt

# what is the total number of overlapping peaks for terra and rloop.
wc -l $project_dir/results/peaks/overlap/terra_overlap_terraAS-S.narrowPeak $project_dir/results/peaks/overlap/rloop_overlap_terraAS-S.narrowPeak

# how many of original terra/rloop peaks with OVERLAP (with at least 1 repeat in column 22) contain at least 4 consecutive telomeric repeats (i.e. three 6's or 7's)?
#cat terraAS-S_overlap_narrowPeak_ann_allRepeats.txt | awk 'BEGIN {FS="\t"}; $22~/TTAGGG/ || $22~/CCCTAA/ {print}' | cut -f1,22 | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l # 174
cat terraAS-S_overlap_narrowPeak_ann_allRepeats.txt | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l # 174

cat rloop_terraAS-S_overlap_narrowPeak_ann_allRepeats.txt | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l # 173

# comparing between how many lines contain telomeric repeats vs how many repeats there are total (grep -o). E.g. the noOverlap file containes 2454 lines with telReps, while the overlap file only contains 264 lines with telReps. However, the total number of repeats in both is very similar (4875 vs 4029)!!! therefore, repeats are much more enriched in the overlapping peaks.
cat terraAS-S_overlap_narrowPeak_ann_allRepeats.txt | cut -f22 | grep -E "(TTAGGG|CCCTAA)" | wc -l # 264
cat terraAS-S_overlap_narrowPeak_ann_allRepeats.txt | cut -f22 | grep -o -E "(TTAGGG|CCCTAA)" | wc -l # 4029
cat rloop_terraAS-S_overlap_narrowPeak_ann_allRepeats.txt | cut -f22 | grep -E "(TTAGGG|CCCTAA)" | wc -l # 236
cat rloop_terraAS-S_overlap_narrowPeak_ann_allRepeats.txt | cut -f22 | grep -o -E "(TTAGGG|CCCTAA)" | wc -l # 3443
################################################################################################ 
################################################################################################ 
################################################################################################ 
################################################################################################ 


################################################################################################ 
# CALL ATRX PEAKS WITH INPUT AS CONTROL
################################################################################################ 

# load some modules
module load bioinfo-tools
module load MACS

# define terra-AS and terra-input alignment bam file variables.
B1=$project_dir/results/alignments/law_paper/SRR057567_Atrx.bwa_mm10.sorted.rg.rmdup.bam
B2=$project_dir/results/alignments/law_paper/SRR080724_Atrx_input.bwa_mm10.sorted.rg.rmdup.bam

# go to terra peaks directory
cd $project_dir/results/peaks/atrx/

# call peaks using MACS2 tool. version 2.2.6. -g: effective genome size or organism (mm for mus musculus). -s: read length. -B/--bdg: store the fragment pileup, control lambda, -log10pvalue and -log10qvalue scores in bedGraph files. -n: The prefix string for output files. -q: q-value (minimum FDR) cutoff, default 0.05. Time 10min.
macs2 callpeak -t $B1 -c $B2 -f BAM -g mm -n SRR057567_Atrx

################################################################################################ 
# GET ATRX PEAKS OVERLAPPING WITH TERRA, THEN CHECK FRACTION OF THESE THAT OVERLAP WITH RLOOP
################################################################################################ 

# load modules.
module load bioinfo-tools
module load BEDOPS
module load BEDTools

# go to directory.
cd $project_dir/results/peaks/overlap/

# define variables for terra and rloop narrowPeak files.
TERRA_PEAKS=$project_dir/results/peaks/terra/SRR2062968_terra_peaks.narrowPeak
RLOOP_PEAKS=$project_dir/results/peaks/rloop/SRR2075686_rloop_peaks.narrowPeak
ATRX_PEAKS=$project_dir/results/peaks/atrx/SRR057567_Atrx_peaks.narrowPeak

# extract terra peaks that overlap with Rloops and Atrx vs Terra peaks that only overlap with Atrx.
intersectBed -wa -a terra_overlap.narrowPeak -b $ATRX_PEAKS | uniq > terra-rloop-overlap_atrx.narrowPeak
intersectBed -wa -a terra_noOverlap.narrowPeak -b $ATRX_PEAKS | uniq > terra-rloop-noOverlap_atrx.narrowPeak

# define variables.
TERRA_RLOOP_OVERLAP_ATRX=$project_dir/results/peaks/overlap/terra-rloop-overlap_atrx.narrowPeak
TERRA_RLOOP_NO_OVERLAP_ATRX=$project_dir/results/peaks/overlap/terra-rloop-noOverlap_atrx.narrowPeak


# get number of R-loop-overlapping/non-overlapping TERRA peaks that overlap/intersect with ATRX peaks.
intersectBed -wa -a terra_overlap.narrowPeak -b $ATRX_PEAKS | uniq | wc # 341
intersectBed -wa -a terra_noOverlap.narrowPeak -b $ATRX_PEAKS | uniq |wc # 231
intersectBed -wa -a $ATRX_PEAKS -b terra_noOverlap.narrowPeak | uniq |wc # 270
intersectBed -wa -a $ATRX_PEAKS -b terra_overlap.narrowPeak | uniq |wc # 554
