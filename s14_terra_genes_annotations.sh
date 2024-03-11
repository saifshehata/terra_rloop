#! /bin/bash
#$ -N terra_genes
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
meme=$project_dir/results/meme
genes=$project_dir/results/genes

mkdir -p $genes
cd $genes

################################################################################################ 
# download excel .xlsx table from publication containing the names of genes differentially 
# expressed upon TERRA knockdown
################################################################################################ 

# get xls file from Jeanie lee's paper. choose one comand.
# wget https://www.cell.com/cms/10.1016/j.cell.2017.06.017/attachment/7c712ac5-cef9-41a4-99ac-c6d5bc53a69a/mmc1.xlsx
wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5552367/bin/NIHMS887216-supplement-8.xlsx -O NIHMS887216-supplement-8.xlsx

# convert xlsx to csv using csvkit
in2csv NIHMS887216-supplement-8.xlsx > NIHMS887216-supplement-8.csv

# extract gene names into a file.
cat NIHMS887216-supplement-8.csv | cut -f1 -d, | sed 1d > terra_knock_down_genes.txt

################################################################################################ 
# get annotation of those peaks with or without 4 tandem telomeric repeats
################################################################################################ 

# define variables for files with names of intersecting peaks having 4 tandem telomeric repeats
TERRA_4t_PEAK_NAMES=$overlaps/terra_intersect_4tandem_peak_names.txt
RLOOP_4t_PEAK_NAMES=$overlaps/rloop_intersect_4tandem_peak_names.txt

# define variables for annotated intersecting peaks
TERRA_INTERSECT_ANN_COUNTS=$annotate/terra_intersect_annotate_motif_counts.txt
RLOOP_INTERSECT_ANN_COUNTS=$annotate/rloop_intersect_annotate_motif_counts.txt

# filter annotation file for only those peaks that have, or lack, 4 tandem repeats i.e. no need to re-annotate
### containing 4 tandem repeats ###
grep -wf $TERRA_4t_PEAK_NAMES $TERRA_INTERSECT_ANN_COUNTS > $annotate/terra_intersect_4tandem_annotate_motif_counts.txt
grep -wf $RLOOP_4t_PEAK_NAMES $RLOOP_INTERSECT_ANN_COUNTS > $annotate/rloop_intersect_4tandem_annotate_motif_counts.txt

### lacking 4 tandem repeats ###
grep -v -wf $TERRA_4t_PEAK_NAMES $TERRA_INTERSECT_ANN_COUNTS > $annotate/terra_intersect_no4tandem_annotate_motif_counts.txt
grep -v -wf $RLOOP_4t_PEAK_NAMES $RLOOP_INTERSECT_ANN_COUNTS > $annotate/rloop_intersect_no4tandem_annotate_motif_counts.txt


#########################################################################################################################################
# from the differentially expressed genes upon TERRA knock-down by Jeanie Lee, get those genes that contain TERRA-R-loop overlapping peaks as well as min 4 tandem repeats 
#########################################################################################################################################

# first, get the names of genes near intersecting peaks from the annotation files
cat $TERRA_INTERSECT_ANN_COUNTS | cut -f16 | sed 1d | sort | uniq | sort > $genes/terra_intersect_annotate_gene_names.txt
cat $RLOOP_INTERSECT_ANN_COUNTS | cut -f16 | sed 1d | sort | uniq | sort > $genes/rloop_intersect_annotate_gene_names.txt

cat $annotate/terra_intersect_4tandem_annotate_motif_counts.txt | cut -f16 | sed 1d | sort | uniq | sort > $genes/terra_intersect_4tandem_annotate_gene_names.txt
cat $annotate/rloop_intersect_4tandem_annotate_motif_counts.txt | cut -f16 | sed 1d | sort | uniq | sort > $genes/rloop_intersect_4tandem_annotate_gene_names.txt

# then, get the genes present in both lists, the 1st being genes near TERRA-R-loop intersecting peaks, and the 2nd being the list of genes differentially expressed upon TERRA knockdown
grep -xf <(cat terra_knock_down_genes.txt) terra_intersect_annotate_gene_names.txt > $genes/terra_intersect_knock_down_genes.txt
grep -xf <(cat terra_knock_down_genes.txt) rloop_intersect_annotate_gene_names.txt > $genes/rloop_intersect_knock_down_genes.txt

grep -xf <(cat terra_knock_down_genes.txt) terra_intersect_4tandem_annotate_gene_names.txt > $genes/terra_intersect_4tandem_knock_down_genes.txt
grep -xf <(cat terra_knock_down_genes.txt) rloop_intersect_4tandem_annotate_gene_names.txt > $genes/rloop_intersect_4tandem_knock_down_genes.txt

################################################################################################ 
# now get the frequency of exons, introns etc. within the intersecting peak regions
################################################################################################ 

# get number of intron, intergenic, exon, UTR, etc. regions that are found at peak regions (all, intersecting, and non-intersecting)
### all peaks ###
cat $annotate/terra_peaks_annotate_motif_counts.txt | cut -f8 | sed 1d | sed -r 's/(\(.+\).*)//' | sort | uniq -c | sort -rn  > $genes/terra_peaks_gene_annotations.txt
cat $annotate/rloop_peaks_annotate_motif_counts.txt | cut -f8 | sed 1d | sed -r 's/(\(.+\).*)//' | sort | uniq -c | sort -rn > $genes/rloop_peaks_gene_annotations.txt

### intersecting peaks ###
cat $TERRA_INTERSECT_ANN_COUNTS | cut -f8 | sed 1d | sed -r 's/(\(.+\).*)//' | sort | uniq -c | sort -rn > $genes/terra_intersect_gene_annotations.txt
cat $RLOOP_INTERSECT_ANN_COUNTS | cut -f8 | sed 1d | sed -r 's/(\(.+\).*)//' | sort | uniq -c | sort -rn > $genes/rloop_intersect_gene_annotations.txt

### intersecting peaks with min 4 tandem repeats ###
cat $annotate/terra_intersect_4tandem_annotate_motif_counts.txt | cut -f8 | sed 1d | sed -r 's/(\(.+\).*)//' | sort | uniq -c | sort -rn > $genes/terra_intersect_4tandem_gene_annotations.txt
cat $annotate/rloop_intersect_4tandem_annotate_motif_counts.txt | cut -f8 | sed 1d | sed -r 's/(\(.+\).*)//' | sort | uniq -c | sort -rn > $genes/rloop_intersect_4tandem_gene_annotations.txt

### intersecting peaks without 4 tandem repeats ###
cat $annotate/terra_intersect_no4tandem_annotate_motif_counts.txt | cut -f8 | sed 1d | sed -r 's/(\(.+\).*)//' | sort | uniq -c | sort -rn > $genes/terra_intersect_no4tandem_gene_annotations.txt
cat $annotate/rloop_intersect_no4tandem_annotate_motif_counts.txt | cut -f8 | sed 1d | sed -r 's/(\(.+\).*)//' | sort | uniq -c | sort -rn > $genes/rloop_intersect_no4tandem_gene_annotations.txt

### non-intersecting peaks ###
cat $annotate/terra_no_intersect_annotate_motif_counts.txt | cut -f8 | sed 1d | sed -r 's/(\(.+\).*)//' | sort | uniq -c | sort -rn > $genes/terra_no_intersect_gene_annotations.txt
cat $annotate/rloop_no_intersect_annotate_motif_counts.txt | cut -f8 | sed 1d | sed -r 's/(\(.+\).*)//' | sort | uniq -c | sort -rn > $genes/rloop_no_intersect_gene_annotations.txt

################################################################################################ 
# now filter the peak annotation files (only intersecting peaks with/without 4 tandem repeats) to get bed files with only those peaks that lie within introns or that are intergenic
# will use these bed files when creating profiles/heatmaps in following scripts
################################################################################################ 
### intersecting peaks with min 4 tandem repeats that are intergenic###
cat $annotate/terra_intersect_4tandem_annotate_motif_counts.txt | sed 1d | grep Intergenic | awk 'OFS="\t" {print $2,$3,$4,$1,$8}' > $genes/terra_intersect_4tandem_intergenic.bed
cat $annotate/rloop_intersect_4tandem_annotate_motif_counts.txt | sed 1d | grep Intergenic | awk 'OFS="\t" {print $2,$3,$4,$1,$8}' > $genes/rloop_intersect_4tandem_intergenic.bed

### intersecting peaks with min 4 tandem repeats that are intronic###
cat $annotate/terra_intersect_4tandem_annotate_motif_counts.txt | sed 1d | grep intron | awk 'OFS="\t" {print $2,$3,$4,$1,$8}' > $genes/terra_intersect_4tandem_intron.bed
cat $annotate/rloop_intersect_4tandem_annotate_motif_counts.txt | sed 1d | grep intron | awk 'OFS="\t" {print $2,$3,$4,$1,$8}' > $genes/rloop_intersect_4tandem_intron.bed

### intersecting peaks without 4 tandem repeats that are intergenic###
cat $annotate/terra_intersect_no4tandem_annotate_motif_counts.txt | sed 1d | grep Intergenic | awk 'OFS="\t" {print $2,$3,$4,$1,$8}' > $genes/terra_intersect_no4tandem_intergenic.bed
cat $annotate/rloop_intersect_no4tandem_annotate_motif_counts.txt | sed 1d | grep Intergenic | awk 'OFS="\t" {print $2,$3,$4,$1,$8}' > $genes/rloop_intersect_no4tandem_intergenic.bed

### intersecting peaks without 4 tandem repeats that are intronic###
cat $annotate/terra_intersect_no4tandem_annotate_motif_counts.txt | sed 1d | grep intron | awk 'OFS="\t" {print $2,$3,$4,$1,$8}' > $genes/terra_intersect_no4tandem_intron.bed
cat $annotate/rloop_intersect_no4tandem_annotate_motif_counts.txt | sed 1d | grep intron | awk 'OFS="\t" {print $2,$3,$4,$1,$8}' > $genes/rloop_intersect_no4tandem_intron.bed

################################################################################################ 
# now do the same filtering for the peak annotation files (only intersecting peaks) to get bed files with only those peaks that lie within introns or that are intergenic
# will use these bed files when creating profiles/heatmaps in following scripts
################################################################################################ 
cat $RLOOP_INTERSECT_ANN_COUNTS | sed 1d | grep Intergenic | awk 'OFS="\t" {print $2,$3,$4,$1,$8}' > $genes/terra_intersect_intergenic.bed
cat $RLOOP_INTERSECT_ANN_COUNTS | sed 1d | grep Intergenic | awk 'OFS="\t" {print $2,$3,$4,$1,$8}' > $genes/rloop_intersect_intergenic.bed

cat $RLOOP_INTERSECT_ANN_COUNTS | sed 1d | grep intron | awk 'OFS="\t" {print $2,$3,$4,$1,$8}' > $genes/terra_intersect_intron.bed
cat $RLOOP_INTERSECT_ANN_COUNTS | sed 1d | grep intron | awk 'OFS="\t" {print $2,$3,$4,$1,$8}' > $genes/rloop_intersect_intron.bed

#############################################


# # get gene names from terra overlapping peak annotation file.
# # cat terra_overlap_narrowPeak_ann_nrRepeats.txt | cut -f16 | sort | uniq | sort > terra_overlap_narrowPeak_ann_geneNames.txt

# ###
# ### intersect with RNA-Seq differentially expressed genes upon Terra-KD from terra paper Jeanie Lee ###
# ###


# # convert from the binary xlsx file to csv. or do this by hand!
# libreoffice --headless --convert-to csv NIHMS887216-supplement-8.xlsx

# # extract gene names into a file.
# cat NIHMS887216-supplement-8.csv | cut -f1 -d, > terra_knockDown_diffExpGenes.txt

# # get the genes from the paper that are in the overlapping peak genes.
# grep -xf <(cat terra_knockDown_diffExpGenes.txt) terra_overlap_narrowPeak_ann_geneNames.txt > terra_overlap_diffExpGenes.txt

# ###
# ### now do the same for the subset of overlapping genes with at least 4 consecutive telomeric repeats.
# ###

# # go to directory
# cd /proj/nb_storage/private/terra_rloop_project/results/peaks/annotation/

# # define variable for overlapping peaks with minimum 4 consecutive repeats.
# TERRA_OVERLAP_4CR=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_min4repeats.narrowPeak

# # annotate these peaks.
# annotatePeaks.pl $TERRA_OVERLAP_4CR mm10 -size given -cpu 10 > terra_overlap_min4repeats_ann.txt

# # get gene names the annotation file.
# cat terra_overlap_min4repeats_ann.txt | cut -f16 | sed 1d | sort | uniq | sort > terra_overlap_min4repeats_ann_geneNames.txt

# # get the genes from the paper that are in the overlapping_min4cR peak genes.
# grep -xf <(cat terra_knockDown_diffExpGenes.txt) terra_overlap_min4repeats_ann_geneNames.txt > terra_overlap_min4cR_diffExpGenes.txt


# cat terra_overlap_diffExpGenes.txt

# #Cnot2
# #Ctnnbl1
# #Dpy19l1
# #Il6
# #Nedd4l
# #Pcnx
# #Pde7b
# #Pola1
# #Sh3bgrl2
# #Tnks

# cat terra_overlap_min4cR_diffExpGenes.txt

# #Nedd4l
# #Pde7b
# #Sh3bgrl2
# #Tnks

# # get info about the diffGenes intersecting with overlapping regions.
# #cat NIHMS887216-supplement-8.csv | grep -wE "(padj|Cnot2|Ctnnbl1|Dpy19l1|Il6|Nedd4l|Pcnx|Pde7b|Pola1|Sh3bgrl2|Tnks)" | sed '1 s/,/geneName,/' | tr "," "\t" | awk 'NR=1{print $0;next}{print $0| "sort -r -k7,7 -V"}' | sort -r -k7,7 -V
# #grep -wf <(cat terra_overlap_diffExpGenes.txt) NIHMS887216-supplement-8.csv | sed '1 s/,/geneName,/' | tr "," "\t" | (sed -u 1q; sort -k7,7 -V) | column -t 

# cat NIHMS887216-supplement-8.csv | grep -wE "(padj|Cnot2|Ctnnbl1|Dpy19l1|Il6|Nedd4l|Pcnx|Pde7b|Pola1|Sh3bgrl2|Tnks)" | sed '1 s/,/geneName,/' | tr "," "\t" | (sed -u 1q; sort -k7,7 -V) | cut -f1,3,7 | tr "\t" "," > terra_overlap_diffExpGenes_info.csv

# # get annotation info (intron/intergenic, etc.) about these genes, peak distance to nearest gene, CpG, GC% and nr of telomeric repeats in peak.
# cat terra_overlap_narrowPeak_ann_nrRepeats.txt | cut -f1,8,9,16 | grep -wE "(padj|Cnot2|Ctnnbl1|Dpy19l1|Il6|Nedd4l|Pcnx|Pde7b|Pola1|Sh3bgrl2|Tnks)"
# cat terra_overlap_narrowPeak_ann_nrRepeats.txt | cut -f1,8,10,16,20-22 | grep -wE "(padj|Cnot2|Ctnnbl1|Dpy19l1|Il6|Nedd4l|Pcnx|Pde7b|Pola1|Sh3bgrl2|Tnks)"

# # extract diffExp genes containing min4repeats from the diff_ExpGenes_info.txt file. the head -1; part is to first print the header, then the output of the grep.
# head -1 terra_overlap_diffExpGenes_info.csv > terra_overlap_min4repeats_diffExpGenes_info.csv
# grep -wf<(cat terra_overlap_min4cR_diffExpGenes.txt) terra_overlap_diffExpGenes_info.csv >> terra_overlap_min4repeats_diffExpGenes_info.csv


# ########################################################################################################################################################
# # INCREASE PEAK LENGTH TO  +/- 500KB AND REPEAT
# ########################################################################################################################################################

# ### now get all genomic sequences +/- 500kb from all overlapping peaks (with 4 consec. repeats), and annotate to get genes, then check for overlap with terraKD differentially expressed genes again. this is to check for genes that terra potentially regulates at enhancer sites.


# # define variables.
# TERRA_OVERLAP_4CR=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_min4repeats.narrowPeak
# TERRA_OVERLAP_NO4CR=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_noConsRepeats.narrowPeak
# TERRA_OVERLAP=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap.narrowPeak

# CHROM_SIZES=/proj/nb_storage/private/terra_rloop_project/refs/mouse/GRCm38/mm10.chrom.sizes

# # go to directory
# cd /proj/nb_storage/private/terra_rloop_project/results/peaks/

# # extend peaks by 500kb in both directions.
# bedtools slop -b 500000 -i $TERRA_OVERLAP -g $CHROM_SIZES | cut -f1-3 > overlap/terra_overlap_slop500kb.bed

# # go to directory
# cd /proj/nb_storage/private/terra_rloop_project/results/peaks/

# # annotate these peaks and get the gene names.
# annotatePeaks.pl overlap/terra_overlap_slop500kb.bed mm10 | cut -f16 | sed 1d | sort | uniq | sort > annotation/terra_overlap_slop500kb_ann_geneNames.txt

# # go to directory
# cd /proj/nb_storage/private/terra_rloop_project/results/peaks/annotation/

# # get the genes from the paper that are in the overlapping_min4cR peak genes.
# grep -xf <(cat terra_knockDown_diffExpGenes.txt) terra_overlap_slop500kb_ann_geneNames.txt > terra_overlap_slop500kb_diffExpGenes.txt


# ###RESULT: no difference in the intersecting genes with the diffExp list from Jeanie Lee's paper.

