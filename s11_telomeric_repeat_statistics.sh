#! /bin/bash
#$ -N telomeric_repeat_stats
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 1

# go to directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
peaks=$project_dir/results/peaks
overlaps=$project_dir/results/overlaps
annotate=$project_dir/results/annotate
repeat_stats=$project_dir/results/repeat_analysis/repeat_stats
mkdir -p $repeat_stats
cd $repeat_stats

################################################################################################ 
# create DATA FRAMES with peak and repeat values for figure generation with ggplot2 in R
################################################################################################ 

# define variables for files with all peaks, intersecting peaks, and non-intersecting peaks, and their annotation files.
### original peaks ###
TERRA_PEAKS=$peaks/terra_sense_peaks_filtered.narrowPeak
RLOOP_PEAKS=$peaks/rloop_peaks_filtered.narrowPeak

### intersecting peaks ###
TERRA_INTERSECT=$overlaps/terra_intersect.narrowPeak
RLOOP_INTERSECT=$overlaps/rloop_intersect.narrowPeak

### non-intersecting peaks ###
TERRA_NO_INTERSECT=$overlaps/terra_no_intersect.narrowPeak
RLOOP_NO_INTERSECT=$overlaps/rloop_no_intersect.narrowPeak

### annotated original peaks ###
TERRA_PEAKS_ANN_POSITIONS=$annotate/terra_peaks_annotate_motif_positions.txt
RLOOP_PEAKS_ANN_POSITIONS=$annotate/rloop_peaks_annotate_motif_positions.txt

### annotated intersecting peaks ###
TERRA_INTERSECT_ANN_POSITIONS=$annotate/terra_intersect_annotate_motif_positions.txt
RLOOP_INTERSECT_ANN_POSITIONS=$annotate/rloop_intersect_annotate_motif_positions.txt

### annotated non-intersecting peaks ###
TERRA_NO_INTERSECT_ANN_POSITIONS=$annotate/terra_no_intersect_annotate_motif_positions.txt
RLOOP_NO_INTERSECT_ANN_POSITIONS=$annotate/rloop_no_intersect_annotate_motif_positions.txt

################################################################################################ 
# assign TERRA peak statistics to variables.
################################################################################################ # go to directory.

NR_TERRA_PEAKS=$(cat $TERRA_PEAKS | wc -l)
NR_TERRA_PEAKS_1REP=$(cat $TERRA_PEAKS_ANN_POSITIONS | cut -f22 | grep -E "(TTAGGG|CCCTAA)" | wc -l)
NR_TERRA_PEAKS_4tREP=$(cat $TERRA_PEAKS_ANN_POSITIONS | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l)
NR_REP_TERRA_PEAKS=$(cat $TERRA_PEAKS_ANN_POSITIONS | cut -f22 | grep -o -E "(TTAGGG|CCCTAA)" | wc -l)
NR_TERRA_PEAKS_1SREP=$(cat $TERRA_PEAKS_ANN_POSITIONS | cut -f9 | grep -E "(TTAGGG|CCCTAA)" | wc -l)
NR_TERRA_PEAKS_4tSREP=$(cat $TERRA_PEAKS_ANN_POSITIONS | awk 'BEGIN {FS="\t"}; $9~/TTAGGG/ || $9~/CCCTAA/ {print}' | cut -f1,22 | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l)

NR_TERRA_INTERSECT=$(cat $TERRA_INTERSECT | wc -l)
NR_TERRA_INTERSECT_1REP=$(cat $TERRA_INTERSECT_ANN_POSITIONS | cut -f22 | grep -E "(TTAGGG|CCCTAA)" | wc -l)
NR_TERRA_INTERSECT_4tREP=$(cat $TERRA_INTERSECT_ANN_POSITIONS | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l)
NR_REP_TERRA_INTERSECT=$(cat $TERRA_INTERSECT_ANN_POSITIONS | cut -f22 | grep -o -E "(TTAGGG|CCCTAA)" | wc -l)
NR_TERRA_INTERSECT_1SREP=$(cat $TERRA_INTERSECT_ANN_POSITIONS | cut -f9 | grep -E "(TTAGGG|CCCTAA)" | wc -l)
NR_TERRA_INTERSECT_4tSREP=$(cat $TERRA_INTERSECT_ANN_POSITIONS | awk 'BEGIN {FS="\t"}; $9~/TTAGGG/ || $9~/CCCTAA/ {print}' | cut -f1,22 | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l)

NR_TERRA_NO_INTERSECT=$(cat $TERRA_NO_INTERSECT | wc -l)
NR_TERRA_NO_INTERSECT_1REP=$(cat $TERRA_NO_INTERSECT_ANN_POSITIONS | cut -f22 | grep -E "(TTAGGG|CCCTAA)" | wc -l)
NR_TERRA_NO_INTERSECT_4tREP=$(cat $TERRA_NO_INTERSECT_ANN_POSITIONS | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l)
NR_REP_TERRA_NO_INTERSECT=$(cat $TERRA_NO_INTERSECT_ANN_POSITIONS | cut -f22 | grep -o -E "(TTAGGG|CCCTAA)" | wc -l)
NR_TERRA_NO_INTERSECT_1SREP=$(cat $TERRA_NO_INTERSECT_ANN_POSITIONS | cut -f9 | grep -E "(TTAGGG|CCCTAA)" | wc -l)
NR_TERRA_NO_INTERSECT_4tSREP=$(cat $TERRA_NO_INTERSECT_ANN_POSITIONS | awk 'BEGIN {FS="\t"}; $9~/TTAGGG/ || $9~/CCCTAA/ {print}' | cut -f1,22 | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l)

################################################################################################ 
# assign RLOOP peak statistics to variables.
################################################################################################ # go to directory.
NR_RLOOP_PEAKS=$(cat $RLOOP_PEAKS | wc -l)
NR_RLOOP_PEAKS_1REP=$(cat $RLOOP_PEAKS_ANN_POSITIONS | cut -f22 | grep -E "(TTAGGG|CCCTAA)" | wc -l)
NR_RLOOP_PEAKS_4tREP=$(cat $RLOOP_PEAKS_ANN_POSITIONS | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l)
NR_REP_RLOOP_PEAKS=$(cat $RLOOP_PEAKS_ANN_POSITIONS | cut -f22 | grep -o -E "(TTAGGG|CCCTAA)" | wc -l)
NR_RLOOP_PEAKS_1SREP=$(cat $RLOOP_PEAKS_ANN_POSITIONS | cut -f9 | grep -E "(TTAGGG|CCCTAA)" | wc -l)
NR_RLOOP_PEAKS_4tSREP=$(cat $RLOOP_PEAKS_ANN_POSITIONS | awk 'BEGIN {FS="\t"}; $9~/TTAGGG/ || $9~/CCCTAA/ {print}' | cut -f1,22 | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l)

NR_RLOOP_INTERSECT=$(cat $RLOOP_INTERSECT | wc -l)
NR_RLOOP_INTERSECT_1REP=$(cat $RLOOP_INTERSECT_ANN_POSITIONS | cut -f22 | grep -E "(TTAGGG|CCCTAA)" | wc -l)
NR_RLOOP_INTERSECT_4tREP=$(cat $RLOOP_INTERSECT_ANN_POSITIONS | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l)
NR_REP_RLOOP_INTERSECT=$(cat $RLOOP_INTERSECT_ANN_POSITIONS | cut -f22 | grep -o -E "(TTAGGG|CCCTAA)" | wc -l)
NR_RLOOP_INTERSECT_1SREP=$(cat $RLOOP_INTERSECT_ANN_POSITIONS | cut -f9 | grep -E "(TTAGGG|CCCTAA)" | wc -l)
NR_RLOOP_INTERSECT_4tSREP=$(cat $RLOOP_INTERSECT_ANN_POSITIONS | awk 'BEGIN {FS="\t"}; $9~/TTAGGG/ || $9~/CCCTAA/ {print}' | cut -f1,22 | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l)

NR_RLOOP_NO_INTERSECT=$(cat $RLOOP_NO_INTERSECT | wc -l)
NR_RLOOP_NO_INTERSECT_1REP=$(cat $RLOOP_NO_INTERSECT_ANN_POSITIONS | cut -f22 | grep -E "(TTAGGG|CCCTAA)" | wc -l)
NR_RLOOP_NO_INTERSECT_4tREP=$(cat $RLOOP_NO_INTERSECT_ANN_POSITIONS | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l)
NR_REP_RLOOP_NO_INTERSECT=$(cat $RLOOP_NO_INTERSECT_ANN_POSITIONS | cut -f22 | grep -o -E "(TTAGGG|CCCTAA)" | wc -l)
NR_RLOOP_NO_INTERSECT_1SREP=$(cat $RLOOP_NO_INTERSECT_ANN_POSITIONS | cut -f9 | grep -E "(TTAGGG|CCCTAA)" | wc -l)
NR_RLOOP_NO_INTERSECT_4tSREP=$(cat $RLOOP_NO_INTERSECT_ANN_POSITIONS | awk 'BEGIN {FS="\t"}; $9~/TTAGGG/ || $9~/CCCTAA/ {print}' | cut -f1,22 | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | wc -l)


################################################################################################ 
# create PEAK STATISTICS data frame
################################################################################################ # go to directory.

# key.
#rep1=min1repeat
#rep0=noRepeats
#rep4c=min4consecutiveRepeats
#rep4c0=no4consecutiveRepeats
#srep=simpleRepeat
#srep4c=simpleRepeatMin4ConsecutiveRepeats

# you can switch the following two expressions for a more simple calculation of Ratio_Repeat_Count without decimal places
# $(printf %.1f "$((10**3 * NR_REP_TERRA_INTERSECT/NR_REP_TERRA_PEAKS*100))e-3")"\t"
# $((10**2 * NR_REP_TERRA_INTERSECT/NR_REP_TERRA_PEAKS))"\t" \
# $(echo | awk '{ print t2/t1*100 }' t2=$NR_REP_TERRA_INTERSECT t1=$NR_REP_TERRA_PEAKS)"\t"

# create header line of data frame/table.
echo -e "Sample\t""Peak_Group\t""Peak_Count\t""No_R\t""Min_1R\t""No_4tR\t""Min_4tR\t""Repeat_Count\t""Ratio_Repeat_Count\t""Silple_Repeat\t""Simple_Repeat_Min_4tR" > peak_statistics.tsv

# create second line.
echo -e "TERRA\t""All\t" \
$NR_TERRA_PEAKS"\t" \
$((NR_TERRA_PEAKS-NR_TERRA_PEAKS_1REP))"\t" \
$NR_TERRA_PEAKS_1REP"\t" \
$((NR_TERRA_PEAKS-NR_TERRA_PEAKS_4tREP))"\t" \
$NR_TERRA_PEAKS_4tREP"\t" \
$NR_REP_TERRA_PEAKS"\t" \
$(echo | awk '{ print t2/t1*100 }' t2=$NR_REP_TERRA_PEAKS t1=$NR_REP_TERRA_PEAKS)"\t" \
$NR_TERRA_PEAKS_1SREP"\t" \
$NR_TERRA_PEAKS_4tSREP >> peak_statistics.tsv

# create third line.
echo -e "TERRA\t""Non-Intersecting\t" \
$NR_TERRA_NO_INTERSECT"\t" \
$((NR_TERRA_NO_INTERSECT-NR_TERRA_NO_INTERSECT_1REP))"\t" \
$NR_TERRA_NO_INTERSECT_1REP"\t" \
$((NR_TERRA_NO_INTERSECT-NR_TERRA_NO_INTERSECT_4tREP))"\t" \
$NR_TERRA_NO_INTERSECT_4tREP"\t" \
$NR_REP_TERRA_NO_INTERSECT"\t" \
$(echo | awk '{ print t2/t1*100 }' t2=$NR_REP_TERRA_NO_INTERSECT t1=$NR_REP_TERRA_PEAKS)"\t" \
$NR_TERRA_NO_INTERSECT_1SREP"\t" \
$NR_TERRA_NO_INTERSECT_4tSREP >> peak_statistics.tsv

# create fourth line.
echo -e "TERRA\t""Intersecting\t" \
$NR_TERRA_INTERSECT"\t" \
$((NR_TERRA_INTERSECT-NR_TERRA_INTERSECT_1REP))"\t" \
$NR_TERRA_INTERSECT_1REP"\t" \
$((NR_TERRA_INTERSECT-NR_TERRA_INTERSECT_4tREP))"\t" \
$NR_TERRA_INTERSECT_4tREP"\t" \
$NR_REP_TERRA_INTERSECT"\t" \
$(echo | awk '{ print t2/t1*100 }' t2=$NR_REP_TERRA_INTERSECT t1=$NR_REP_TERRA_PEAKS)"\t" \
$NR_TERRA_INTERSECT_1SREP"\t" \
$NR_TERRA_INTERSECT_4tSREP >> peak_statistics.tsv

# create fifth line.
echo -e "R-loop\t""All\t" \
$NR_RLOOP_PEAKS"\t" \
$((NR_RLOOP_PEAKS-NR_RLOOP_PEAKS_1REP))"\t" \
$NR_RLOOP_PEAKS_1REP"\t" \
$((NR_RLOOP_PEAKS-NR_RLOOP_PEAKS_4tREP))"\t" \
$NR_RLOOP_PEAKS_4tREP"\t" \
$NR_REP_RLOOP_PEAKS"\t" \
$(echo | awk '{ print t2/t1*100 }' t2=$NR_REP_RLOOP_PEAKS t1=$NR_REP_RLOOP_PEAKS)"\t" \
$NR_RLOOP_PEAKS_1SREP"\t" \
$NR_RLOOP_PEAKS_4tSREP >> peak_statistics.tsv

# create sixth line.
echo -e "R-loop\t""Non-Intersecting\t" \
$NR_RLOOP_NO_INTERSECT"\t" \
$((NR_RLOOP_NO_INTERSECT-NR_RLOOP_NO_INTERSECT_1REP))"\t" \
$NR_RLOOP_NO_INTERSECT_1REP"\t" \
$((NR_RLOOP_NO_INTERSECT-NR_RLOOP_NO_INTERSECT_4tREP))"\t" \
$NR_RLOOP_NO_INTERSECT_4tREP"\t" \
$NR_REP_RLOOP_NO_INTERSECT"\t" \
$(echo | awk '{ print t2/t1*100 }' t2=$NR_REP_RLOOP_NO_INTERSECT t1=$NR_REP_RLOOP_PEAKS)"\t" \
$NR_RLOOP_NO_INTERSECT_1SREP"\t" \
$NR_RLOOP_NO_INTERSECT_4tSREP >> peak_statistics.tsv

# create seventh line.
echo -e "R-loop\t""Intersecting\t" \
$NR_RLOOP_INTERSECT"\t" \
$((NR_RLOOP_INTERSECT-NR_RLOOP_INTERSECT_1REP))"\t" \
$NR_RLOOP_INTERSECT_1REP"\t" \
$((NR_RLOOP_INTERSECT-NR_RLOOP_INTERSECT_4tREP))"\t" \
$NR_RLOOP_INTERSECT_4tREP"\t" \
$NR_REP_RLOOP_INTERSECT"\t" \
$(echo | awk '{ print t2/t1*100 }' t2=$NR_REP_RLOOP_INTERSECT t1=$NR_REP_RLOOP_PEAKS)"\t" \
$NR_RLOOP_INTERSECT_1SREP"\t" \
$NR_RLOOP_INTERSECT_4tSREP >> peak_statistics.tsv



################################################################################################ 
# create PEAK REPEATS data frame
################################################################################################ 

# create transposed and modified data frame only containing subset of first data frame. to be used for showing repeat statistics.
cat peak_statistics.tsv | sed "1 d" |\
    awk 'BEGIN {print "Sample\t""Peak_Group\t""Peak_Count\t""Ratio_Peak_Count\t""Repeat_Group"}
    {OFS="\t";
            {print $1,$2,$3,$3/$3*100,"All"} 
            {print $1,$2,$4,$4/$3*100,"No_R"} 
            {print $1,$2,$5,$5/$3*100,"Min_1R"} 
            {print $1,$2,$6,$6/$3*100,"No_4tR"} 
            {print $1,$2,$7,$7/$3*100,"Min_4tR"}
    }' \
> peak_repeats.tsv


# add the following to reduce the Ratio column to one decimal place: print the 1st line of the stdin, then do stuff on column 5 ($5) but skip the first line of that output, then print it.
#     awk 'FNR == 1 {print} FNR == 1 {next}{$5=sprintf("%.1f",$5)}1' \


# example.
#echo -n -e $()"\t" >> peak_statistics.tsv
#printf "%s \n %s" "$HEADER $(cat peak_repeats.tsv)"
#echo -e $HEADER "\n" "$(cat peak_repeats.tsv)" > peak_repeats.tsv







