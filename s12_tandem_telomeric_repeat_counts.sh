#! /bin/bash
#$ -N tandem_repeats
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 1

# go to directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
peaks=$project_dir/results/peaks
# overlaps=$project_dir/results/overlaps
annotate=$project_dir/results/annotate
tandem_repeats=$project_dir/results/repeat_analysis/tandem_repeats
mkdir -p $tandem_repeats
cd $tandem_repeats

################################################################################################ 
# create FILES with maximum consecutive repeat counts per peak and nr of groups of at least 4 consecutive repeats per peak
# for figure generation with ggplot2 in R
################################################################################################ 

# define variables for peaks annotation files (all peaks, intersecting peaks, and non-intersecting peaks) containing all individual telomeric repeats with their positions
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
# create FILES with maximum repeat counts per peak and nr of groups of at least 4 consecutive repeats per peak
################################################################################################ 

# key
# column 1: peakID
# column 2: total number of grouped consecutive repeats (i.e. separated by more than 6 or 7 bp in the genome)
# column 3: maximum number of consecutive repeats per peak
# column 4: second largest number of consecutive repeats in that peak
# column 5: third largest number of ...

# using file with ALL TERRA peaks.
cat $TERRA_PEAKS_ANN_POSITIONS | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" |\
    sed -r 's|(\(.{13}\))||g'| tr "," "\t" |\
    awk '{c=1; 
        for (i=2;i<=NF;i++) {
            if(($(i+1)-$i==6 || $(i+1)-$i==7)) {c++; $i=""} 
            else if(i==NF) {$i=c} 
            else {$i=c","; c=1}
        }
    }1' |\
    awk '{for(i=1;i<=NF;i++) printf $i""FS ; print ""}' |\
    sed 's| |\t|' | sed 's| ||g' |\
    awk '{printf $1}{
        split($2,a,","); 
        asort(a,a,"@val_num_desc"); 
        for(i=1; i<=length(a); i++) printf("\t%s", a[i]); printf("\n"); 
        }' |\
    sort -k2 -r -V |\
    awk '{OFS="\t"; c=0;
        for(i=2;i<=NF;i++) {
            if($i >= 4) {c++}  
        } 
        {$(NF+1)=c} {print $NF,$0}
    }' |\
    awk 'OFS="\t" {$NF=""}1' |\
    awk '{ temp = $1; $1 = $2; $2 = temp } 1' OFS='\t' |\
    awk 'BEGIN {print "Peak_ID\t""Tandem_Repeat_Groups\t""Max_Tandem_Repeats"}{OFS="\t"; print $1,$2,$3}' \
> terra_peaks_tandem_motif_counts.txt

# using file with INTERSECTING TERRA peaks.
cat $TERRA_INTERSECT_ANN_POSITIONS | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" |\
    sed -r 's|(\(.{13}\))||g'| tr "," "\t" |\
    awk '{c=1; 
        for (i=2;i<=NF;i++) {
            if(($(i+1)-$i==6 || $(i+1)-$i==7)) {c++; $i=""} 
            else if(i==NF) {$i=c} 
            else {$i=c","; c=1}
        }
    }1' |\
    awk '{for(i=1;i<=NF;i++) printf $i""FS ; print ""}' |\
    sed 's| |\t|' | sed 's| ||g' |\
    awk '{printf $1}{
        split($2,a,","); 
        asort(a,a,"@val_num_desc"); 
        for(i=1; i<=length(a); i++) printf("\t%s", a[i]); printf("\n"); 
        }' |\
    sort -k2 -r -V |\
    awk '{OFS="\t"; c=0;
        for(i=2;i<=NF;i++) {
            if($i >= 4) {c++}  
        } 
        {$(NF+1)=c} {print $NF,$0}
    }' |\
    awk 'OFS="\t" {$NF=""}1' |\
    awk '{ temp = $1; $1 = $2; $2 = temp } 1' OFS='\t' |\
    awk 'BEGIN {print "Peak_ID\t""Tandem_Repeat_Groups\t""Max_Tandem_Repeats"}{OFS="\t"; print $1,$2,$3}' \
> terra_intersect_tandem_motif_counts.txt

# using file with NON_INTERSECTING TERRA peaks.
cat $TERRA_NO_INTERSECT_ANN_POSITIONS | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" |\
    sed -r 's|(\(.{13}\))||g'| tr "," "\t" |\
    awk '{c=1; 
        for (i=2;i<=NF;i++) {
            if(($(i+1)-$i==6 || $(i+1)-$i==7)) {c++; $i=""} 
            else if(i==NF) {$i=c} 
            else {$i=c","; c=1}
        }
    }1' |\
    awk '{for(i=1;i<=NF;i++) printf $i""FS ; print ""}' |\
    sed 's| |\t|' | sed 's| ||g' |\
    awk '{printf $1}{
        split($2,a,","); 
        asort(a,a,"@val_num_desc"); 
        for(i=1; i<=length(a); i++) printf("\t%s", a[i]); printf("\n"); 
        }' |\
    sort -k2 -r -V |\
    awk '{OFS="\t"; c=0;
        for(i=2;i<=NF;i++) {
            if($i >= 4) {c++}  
        } 
        {$(NF+1)=c} {print $NF,$0}
    }' |\
    awk 'OFS="\t" {$NF=""}1' |\
    awk '{ temp = $1; $1 = $2; $2 = temp } 1' OFS='\t' |\
    awk 'BEGIN {print "Peak_ID\t""Tandem_Repeat_Groups\t""Max_Tandem_Repeats"}{OFS="\t"; print $1,$2,$3}' \
> terra_no_intersect_tandem_motif_counts.txt

# using file with ALL RLOOP peaks.
cat $RLOOP_PEAKS_ANN_POSITIONS | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" |\
    sed -r 's|(\(.{13}\))||g'| tr "," "\t" |\
    awk '{c=1; 
        for (i=2;i<=NF;i++) {
            if(($(i+1)-$i==6 || $(i+1)-$i==7)) {c++; $i=""} 
            else if(i==NF) {$i=c} 
            else {$i=c","; c=1}
        }
    }1' |\
    awk '{for(i=1;i<=NF;i++) printf $i""FS ; print ""}' |\
    sed 's| |\t|' | sed 's| ||g' |\
    awk '{printf $1}{
        split($2,a,","); 
        asort(a,a,"@val_num_desc"); 
        for(i=1; i<=length(a); i++) printf("\t%s", a[i]); printf("\n"); 
        }' |\
    sort -k2 -r -V |\
    awk '{OFS="\t"; c=0;
        for(i=2;i<=NF;i++) {
            if($i >= 4) {c++}  
        } 
        {$(NF+1)=c} {print $NF,$0}
    }' |\
    awk 'OFS="\t" {$NF=""}1' |\
    awk '{ temp = $1; $1 = $2; $2 = temp } 1' OFS='\t' |\
    awk 'BEGIN {print "Peak_ID\t""Tandem_Repeat_Groups\t""Max_Tandem_Repeats"}{OFS="\t"; print $1,$2,$3}' \
> rloop_peaks_tandem_motif_counts.txt

# using file with INTERSECTING RLOOP peaks.
cat $RLOOP_INTERSECT_ANN_POSITIONS | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" |\
    sed -r 's|(\(.{13}\))||g'| tr "," "\t" |\
    awk '{c=1; 
        for (i=2;i<=NF;i++) {
            if(($(i+1)-$i==6 || $(i+1)-$i==7)) {c++; $i=""} 
            else if(i==NF) {$i=c} 
            else {$i=c","; c=1}
        }
    }1' |\
    awk '{for(i=1;i<=NF;i++) printf $i""FS ; print ""}' |\
    sed 's| |\t|' | sed 's| ||g' |\
    awk '{printf $1}{
        split($2,a,","); 
        asort(a,a,"@val_num_desc"); 
        for(i=1; i<=length(a); i++) printf("\t%s", a[i]); printf("\n"); 
        }' |\
    sort -k2 -r -V |\
    awk '{OFS="\t"; c=0;
        for(i=2;i<=NF;i++) {
            if($i >= 4) {c++}  
        } 
        {$(NF+1)=c} {print $NF,$0}
    }' |\
    awk 'OFS="\t" {$NF=""}1' |\
    awk '{ temp = $1; $1 = $2; $2 = temp } 1' OFS='\t' |\
    awk 'BEGIN {print "Peak_ID\t""Tandem_Repeat_Groups\t""Max_Tandem_Repeats"}{OFS="\t"; print $1,$2,$3}' \
> rloop_intersect_tandem_motif_counts.txt

# using file with NON-INTERSECTING RLOOP peaks.
cat $RLOOP_NO_INTERSECT_ANN_POSITIONS | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" |\
    sed -r 's|(\(.{13}\))||g'| tr "," "\t" |\
    awk '{c=1; 
        for (i=2;i<=NF;i++) {
            if(($(i+1)-$i==6 || $(i+1)-$i==7)) {c++; $i=""} 
            else if(i==NF) {$i=c} 
            else {$i=c","; c=1}
        }
    }1' |\
    awk '{for(i=1;i<=NF;i++) printf $i""FS ; print ""}' |\
    sed 's| |\t|' | sed 's| ||g' |\
    awk '{printf $1}{
        split($2,a,","); 
        asort(a,a,"@val_num_desc"); 
        for(i=1; i<=length(a); i++) printf("\t%s", a[i]); printf("\n"); 
        }' |\
    sort -k2 -r -V |\
    awk '{OFS="\t"; c=0;
        for(i=2;i<=NF;i++) {
            if($i >= 4) {c++}  
        } 
        {$(NF+1)=c} {print $NF,$0}
    }' |\
    awk 'OFS="\t" {$NF=""}1' |\
    awk '{ temp = $1; $1 = $2; $2 = temp } 1' OFS='\t' |\
    awk 'BEGIN {print "Peak_ID\t""Tandem_Repeat_Groups\t""Max_Tandem_Repeats"}{OFS="\t"; print $1,$2,$3}' \
> rloop_no_intersect_tandem_motif_counts.txt






################################################################################################ 
# SAME CODE IN EXPANDED VIEW FOR CLARITY. EXPLANATION BELOW CODE
################################################################################################ 
# cat $VARIABLE_NAME | cut -f1,22 |grep -E "(TTAGGG|CCCTAA)" |\

# cat $RLOOP_INTERSECT_ANN_POSITIONS | cut -f1,22 |grep -E "(TTAGGG|CCCTAA)" |\
#     sed -r 's|(\(.{13}\))||g'| tr "," "\t" |\
#     awk '{c=1; 
#         for (i=2;i<=NF;i++) {
#             if(($(i+1)-$i==6 || $(i+1)-$i==7)) {c++; $i=""} 
#             else if(i==NF) {$i=c} 
#             else {$i=c","; c=1}
#         }
#     }1' |\
#     awk '{for(i=1;i<=NF;i++) printf $i""FS ; print ""}' |\
#     sed 's| |\t|' | sed 's| ||g' |\
#     awk '{printf $1}{
#         split($2,a,","); 
#         asort(a,a,"@val_num_desc"); 
#         for(i=1; i<=length(a); i++) printf("\t%s", a[i]); printf("\n"); 
#         }' |\
#     sort -k2 -r -V |\
#     awk '{OFS="\t"; c=0;
#         for(i=2;i<=NF;i++) {
#             if($i >= 4) {c++}  
#         } 
#         {$(NF+1)=c} {print $NF,$0}
#     }' |\
#     awk 'OFS="\t" {$NF=""}1' |\
#     awk '{ temp = $1; $1 = $2; $2 = temp } 1' OFS='\t' |\
#     awk 'BEGIN {print "Peak_ID\t""Tandem_Repeat_Groups\t""Max_Tandem_Repeats"}{OFS="\t"; print $1,$2,$3}' \
# > file_tandem_motif_counts.txt

################################################################################################ 
# EXPLANATIONS
################################################################################################ 

# cat: open file to stout.
# cut: extract specific columns.
# grep: only keep lines that contain the repeat or its reverse complement.
# sed: from the annotation in column 22, change e.g. `93(TTAGGG,+,0.00)` to `93`.
# tr: replace all commas with tabs.
# awk: count and report how many consecutive repeats there are per peak. For each stretch of consecutive repeats, only report the final count and remove all other numbers. Consecutive repeats are defined by me as either directly after one another, e.g. TTAGGGTTAGGG (i.e. 6 bp apart), or with one base in between, e.g. TTAGGGNTTAGGG (i.e. 7 bp apart).
# awk: removes all fields with empty characters.
# sed: replace the first space it finds with tab
# sed: replace all spaces it finds with nothing (i.e. remove them).
# awk: sorts the second column (containing counts of grouped consecutive repeats) in descending order so the count of the largest group comes first, then second largest, etc. (first creates an array 'a' containing all its values, then sorts).
########## start: array details ##########
# print column 1 (peak names)
# then split column 2 (nr of consecutive repeats) on delimiter ',' into array called 'a'.
# then sort array 'a' numerically by descending values and re-assign it to 'a'.
# then for each value in the sorted array, print it, with a tab preceeding it, then pring a newline at the end of the line so the next row starts on another line.
########## end: array details ##########
# sort: sort on second column in descending order using -V, --version-sort (i.e. sort on highest count of grouped consecutive repeats, then on 2nd highest, etc.).
# awk: count and report the number of groups of 4+ consecutive repeats and report it in the last column, then duplicate the last column before the first column.
# awk: remove the last column, since there is a copy of it in the first column).
# awk: swap the 1st and 2nd columns together, so the peak ID/name is the 1st column, and the number of grouped repeats in the 2nd column.
# awk: begin by printing the header line, then only select the first 3 columns: peakID, nr of consecutive (4+) repeat groups per peak, and max nr of consecutive repeats per peak.




# awk: starting from the 4th column onwards, for each line, remove any field where the value is less than 3. This is because there are a lot of 1's and 2's in some fields which means some rows have too many fields, and this causes issues when reading files in R using eg. fread or data.table or read.csv.


########## start: optional add at the end before less -S ##########
#     sort -k2 -rV | awk '$2>1{print}' |\
# sort: sort descending on number of groups of consecutive repeats.
# awk: only keep peaks with more than 1 group of consecutive repeats.
########## end: optional add at the end before less -S ##########


# example.
#awk '{sub(/regexp/, replacement, target)}'
#awk '{sub(/_[1-9]+/," ", "\t")}1'
#sed 's/regexp/replacement/flags'
