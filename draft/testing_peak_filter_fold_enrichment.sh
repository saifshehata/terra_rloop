#! /bin/bash
#$ -N test
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 2

# # testing_peak_filter_fold_enrichment

# cat peaks/terra_sense_peaks.narrowPeak |  wc -l
# 27707
# cat peaks/terra_sense_peaks.narrowPeak | awk '$7 >= 10{print $0}' | wc -l
# 3758


# cat overlaps/terra_intersect.narrowPeak | wc -l
# 689
# cat overlaps/terra_intersect.narrowPeak | awk '$7>=10{print$0}'| wc -l
# 294


# cat overlaps/terra_no_intersect.narrowPeak | wc -l
# 26936
# cat overlaps/terra_no_intersect.narrowPeak | awk '$7>=10{print$0}'| wc -l
# 3441


# cat overlaps/terra_intersect_4tandem.narrowPeak | wc -l
# 176
# cat overlaps/terra_intersect_4tandem.narrowPeak | awk '$7>=10{print$0}'| wc -l
# 160


# cat overlaps/terra_intersect_no4tandem.narrowPeak | wc -l
# 513
# cat overlaps/terra_intersect_no4tandem.narrowPeak | awk '$7>=10{print$0}'| wc -l
# 134


cat fimo.tsv | cut -f3- | egrep '(TTAGGGTTAGGGTTAGGGTTAGGG|ttagggttagggttagggttaggg)' | cut -f1-3 | sort -k1V |\
awk '{a[$1]=a[$1] FS $2} END{for(i in a) print i a[i]}' |\
awk 'BEGIN {if($3-$2!=6) {$2=$2" "$2+23} } {
    for (i=3;i<NF;i++) {
        if($(i+1)-$i==6) {$i=""} 
        else if($(i-1)=="") {$i=$i+23 " " $(i+1)}
        else {$i=$i " " $(i+1)}
        } 
    }1' |\
    awk '{for(i=1;i<=NF;i++) printf $i""FS ; print ""}' |\
    awk '{ for (i=2;i<NF;i=i+2) {print$1, $i, $(i+1)} }' \
> fimo.bed


# | awk '{$NF+=23}1'
