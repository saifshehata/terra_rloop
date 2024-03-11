# save cluster_1 and cluster_2 after deeptools heatmap of terra read coverage on a) terra intersecting peaks, b) terra intersecting peaks with 4tandem repeats, and c) terra intersecting peaks without 4tandem repeats
paste <(cat results/profile/terra_reads_on_intersecting_peaks.bed | egrep '(cluster_1|cluster_2)' | cut -f4,13 | sort -k2,2V -k1,1V) <(cat results/profile/terra_reads_on_intersecting_peaks_4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f4,13 | sort -k2,2V -k1,1V) <(cat results/profile/terra_reads_on_intersecting_peaks_no4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f4,13 | sort -k2,2V -k1,1V) > terra_clusters.txt

# get the peaks from clusters 1 and 2 in b) that are present in a)
grep -wf <(cat results/profile/terra_reads_on_intersecting_peaks_4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f4,13 | sort -k2,2V -k1,1V | cut -f1) terra_clusters.txt

# get the peaks from clusters 1 and 2 in c) that are present in a)
grep -wf <(cat results/profile/terra_reads_on_intersecting_peaks_no4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f4,13 | sort -k2,2V -k1,1V | cut -f1) terra_clusters.txt


# get the nr of peaks from clusters 1 and 2 in b) that are present in a)
grep -wf <(cat results/profile/terra_reads_on_intersecting_peaks_4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f4,13 | sort -k2,2V -k1,1V | cut -f1) <(cut -f1,2 terra_clusters.txt) | wc -l

# get the nr of peaks from clusters 1 and 2 in c) that are present in a)
grep -wf <(cat results/profile/terra_reads_on_intersecting_peaks_no4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f4,13 | sort -k2,2V -k1,1V | cut -f1) <(cut -f1,2 terra_clusters.txt) | wc -l

# get chr start end and peak name for clusters 1 and 2 to later check if the enrichment of TERRA peaks at these peak regions is artificial
cat results/profile/terra_reads_on_intersecting_peaks_no4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f1-4,13 | sort -k5,5V -k1,1V
# terra_sense_peak_3203

# get gaps/centromere regions in mm10 reference genome
curl https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz | zcat | column -t | less
curl https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/gap.txt.gz | zcat | column -t | less

# get blacklisted regions
curl http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz | zcat | column -t | less
https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz?raw=true
https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz
https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz


# get the high coverage intersecting peaks without 4 tandem repeats that intersect with blacklisted regions
intersectBed -wo -a <(cat results/profile/terra_reads_on_intersecting_peaks_no4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f1-4) -b <(curl http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz | zcat) | column -t
intersectBed -wo -a <(cat results/profile/terra_reads_on_intersecting_peaks_no4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f1-4) -b <(curl https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz | zcat) | column -t



# get the high coverage intersecting peaks with 4 tandem repeats that intersect with blacklisted regions
intersectBed -wo -a <(cat results/profile/terra_reads_on_intersecting_peaks_4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f1-4) -b <(curl http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz | zcat) | column -t


# get the high coverage intersecting peaks with 4 tandem repeats that intersect with gapped regions
intersectBed -wo -a <(cat results/profile/terra_reads_on_intersecting_peaks_no4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f1-4) -b <(curl https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/gap.txt.gz | zcat | cut -f2-) | column -t

# get the high coverage intersecting peaks with 4 tandem repeats that intersect with cytoband regions
intersectBed -wo -a <(cat results/profile/terra_reads_on_intersecting_peaks_no4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f1-4) -b <(curl https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz | zcat) | column -t


#grep -wf <(peak names in cluster_1 and cluster_2 of bed file in /profile) <(peak annotation file)
grep -wf <(cat results/profile/terra_reads_on_intersecting_peaks_4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f4) results/annotate/terra_peaks_annotate_motif_counts.txt | cut -f16


paste <(cat results/profile/terra_reads_on_intersecting_peaks.bed | egrep '(cluster_1|cluster_2)' | cut -f4,13 | sort -k2,2V -k1,1V) \
<(cat results/profile/terra_reads_on_intersecting_peaks_4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f4,13 | sort -k2,2V -k1,1V) \
<(cat results/profile/rloop_reads_on_intersecting_peaks.bed | egrep '(cluster_1|cluster_2)' | cut -f4,13 | sort -k2,2V -k1,1V) \
<(cat results/profile/rloop_reads_on_intersecting_peaks_4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f4,13 | sort -k2,2V -k1,1V) \
<(cat results/profile/atrx_reads_on_intersecting_peaks.bed | egrep '(cluster_1|cluster_2)' | cut -f4,13 | sort -k2,2V -k1,1V) \
<(cat results/profile/atrx_reads_on_intersecting_peaks_4tandem.bed | egrep '(cluster_1|cluster_2)' | cut -f4,13 | sort -k2,2V -k1,1V)


paste <(cat results/profile/terra_reads_on_intersecting_peaks.bed | egrep '(cluster_1)' | cut -f4,13 | sort -k2,2V -k1,1V) \
<(cat results/profile/rloop_reads_on_intersecting_peaks.bed | egrep '(cluster_1)' | cut -f4,13 | sort -k2,2V -k1,1V) \
<(cat results/profile/atrx_reads_on_intersecting_peaks.bed | egrep '(cluster_1)' | cut -f4,13 | sort -k2,2V -k1,1V)

paste <(cat results/profile/terra_reads_on_intersecting_peaks_4tandem.bed | egrep '(cluster_1)' | cut -f4,13 | sort -k2,2V -k1,1V) \
<(cat results/profile/rloop_reads_on_intersecting_peaks_4tandem.bed | egrep '(cluster_1)' | cut -f4,13 | sort -k2,2V -k1,1V) \
<(cat results/profile/atrx_reads_on_intersecting_peaks_4tandem.bed | egrep '(cluster_1)' | cut -f4,13 | sort -k2,2V -k1,1V)


# get ensemble gene names of genes closest to/associated with peaks with highest terra coverage from cluster_1 of heatmap 
grep -wf <(cat terra_reads_on_intersecting_peaks.bed | grep cluster_1 | cut -f4) $annotate/terra_intersect_annotate_motif_counts.txt | cut -f15