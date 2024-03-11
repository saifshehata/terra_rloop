#! /bin/bash
#$ -N genome_track
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 4

######################################################################################################################################################
# create genome track images
######################################################################################################################################################

# go to directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
mapped_reads=$project_dir/results/mapped_reads
annotate=$project_dir/results/annotate
genome=$project_dir/raw_data/genome
profile=$project_dir/results/profile
genome_tracks=$project_dir/results/figures/genome_tracks

mkdir -p $genome_tracks
cd $genome_tracks

# define variables
TERRA_BAM=$mapped_reads/SRR2062968_pe_sort_rmdup.bam
RLOOP_BAM=$mapped_reads/SRR2075686_se_sort_rmdup.bam
ATRX_BAM=$mapped_reads/SRR057567_se_sort_rmdup.bam

GFF=$genome/GCA_000001635.5_GRCm38.p3_full_analysis_set.refseq_annotation.gff.gz
#GTF=/home/saif/terra_rloop_project/genome_browser/mm10.ncbiRefSeq.gtf

################# using svist4get to visualize genome tracks. must create bedgraph files as it does not like bigwig files #################

# create bedgraph file.
# bedtools genomecov -ibam $TERRA_BAM -bg > $mapped_reads/SRR2062968_pe_sort_rmdup.bg
# bedtools genomecov -ibam $RLOOP_BAM -bg > $mapped_reads/SRR2075686_se_sort_rmdup.bg
# bedtools genomecov -ibam $ATRX_BAM -bg > $mapped_reads/SRR057567_se_sort_rmdup.bg

# doanload gtf files
wget -P $genome --timestamping http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz
wget -P $genome --timestamping http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz

# First filter the GTF file by removing XM_ and XR_ accessions to avoid all of them being incorporated in the final figure when using genomic regions (-w) instead of transcript IDs (-t).
zcat $GFF | egrep -v '(XM_|XR_)' > $genome/GCA_000001635.5_GRCm38.p3_full_analysis_set.refseq_annotation_filtered.gff
zcat $genome/mm10.ncbiRefSeq.gtf.gz | egrep -v '(XM_|XR_)' > $genome/mm10.ncbiRefSeq_filtered.gtf
zcat $genome/mm10.refGene.gtf.gz | egrep -v '(XM_|XR_)' > $genome/mm10.refGene_filtered.gtf


# define variables
TERRA_BG=$mapped_reads/SRR2062968_pe_sort_rmdup.bg
RLOOP_BG=$mapped_reads/SRR2075686_se_sort_rmdup.bg
ATRX_BG=$mapped_reads/SRR057567_se_sort_rmdup.bg

REPEATS=$annotate/terra_intersect_motif_positions.bed
GFF_FILTERED=$genome/GCA_000001635.5_GRCm38.p3_full_analysis_set.refseq_annotation_filtered.gff
GTF_FILTERED=$genome/mm10.ncbiRefSeq_filtered.gtf
GTF2_FILTERED=$genome/mm10.refGene_filtered.gtf
REF=$genome/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna

# Creates the 'svist4get_data' folder in the current working directory. The folder will contain sample data sets and adjustable configuration file templates.
svist4get --sampledata

# nedd4l whole gene
# svist4get -gtf $GTF -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -t NM_001114386.1 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'Nedd4l gene' -gi 65099500-65109500 -gil '' -hf 65099500 65109500 'Peak' -bgb max -o nedd4l_gene.pdf

svist4get -gtf $GTF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr18 64887755 65217826 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'Nedd4l gene' -gi 65099500-65109500 -gil '' -hf 65099500 65109500 'Peak' -bgb max -c svist4get_data/A4_p1.cfg -o nedd4l_gene_gtf1.pdf

svist4get -gtf $GTF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr18 65103965 65104162 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'Nedd4l peak' -gi 65104021-65104087 -gi 65103976-65103982 -gi 65104003-65104009 -gi 65104113-65104119 -gi 65104145-65104151 -gil '' '' '' '' '' -hf 65104021 65104087 'Telomeric Repeats' -hf 65103976 65103982 '' -hf 65104003 65104009 '' -hf 65104113 65104119 '' -hf 65104145 65104151 '' -bgb max -nts -c svist4get_data/A4_p1.cfg -o nedd4l_peak_gtf1.pdf

# # get telomeric repeat positions from the bed file you created from HOMER annotatePeaks.pl.
# cat $REPEATS | grep chr18 | grep 651

# Tnks gene track zoomOut.
svist4get -gtf $GTF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr8 34826000 35148000 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'Tnks gene' -gi 35130000-35139000 -gil '' -hf 35130000 35139000 'Peak' -bgb max -c svist4get_data/A4_p1.cfg -o tnks_gene.pdf

# tnks rloop peak region.
svist4get -gtf $GTF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr8 35134877 35135032 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'Tnks proximal peak' -gi 35134903-35134993 -gil '' -hf 35134903 35134993 'Telomeric Repeats' -bgb max -nts -c svist4get_data/A4_p1.cfg -o tnks_peak.pdf

# # get telomeric repeat positions from the bed file you created from HOMER annotatePeaks.pl.
# cat $REPEATS | grep chr8 | grep 351


# # removing the auto visualization of amino acod track.

# SVIST_PATH=$(which svist4get) # /cluster/khiom/sshehata001/miniconda3/envs/ngs/bin/svist4get
# sudo nano $SVIST_PATH
# # change  "parameters.config['show_aa_seq_track'] == 'auto' to "parameters.config['show_aa_seq_track'] = 0"
# # save changes and rerun code


##############################################################################
# min 4 repeats intron
##############################################################################

# map2k6
# chr11:110422188-110422411

cat $REPEATS | grep chr11 | grep 110

# map2k6 rloop peak region.
svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr11 110422240 110422390 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'MAP2K6 peak' -gi 110422264-110422330 -gi 110422367-110422373 -gil '' '' -hf 110422264 110422330 'Telomeric Repeats' -hf 110422367 110422373 '' -bgb max -nts -c svist4get_data/A4_p1.cfg -o map2k6_peak.pdf

# bard1
# SRR2075686_rloop_peak_637
# chr1:71079268-71079495

cat $REPEATS | grep chr1 | grep 71079

svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr1 71079316 71079435 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'BARD1 peak' -gi 71079344-71079410 -gil '' -hf 71079344 71079410 'Telomeric Repeats' -hf 71079344 71079410 '' -bgb max -nts -c svist4get_data/A4_p1.cfg -o bard1_peak.pdf

# pde7b
# SRR2075686_rloop_peak_1952
# chr10:20459730-20459984
# chr10:20,459,781-20,459,927

cat $REPEATS | grep chr10 | grep 20459

svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr10 20459771 20459937 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'PDE7B peak' -gi 20459783-20459789 -gi 20459795-20459801 -gi 20459813-20459843 -gi 20459855-20459891 -gi 20459903-20459921 -gil '' '' '' '' '' -hf 20459783 20459789 '' -hf 20459795 20459801 '' -hf 20459813 20459843 '' -hf 20459855 20459891 'Telomeric Repeats' -hf 20459903 20459921 '' -bgb max -nts -c svist4get_data/A4_p1.cfg -o pde7b_peak.pdf

chr10	20459915	20459921	telomeric_repeat	8.299739	-
chr10	20459909	20459915	telomeric_repeat	8.299739	-
chr10	20459903	20459909	telomeric_repeat	8.299739	-

chr10	20459885	20459891	telomeric_repeat	8.299739	-
chr10	20459879	20459885	telomeric_repeat	8.299739	-
chr10	20459873	20459879	telomeric_repeat	8.299739	-
chr10	20459867	20459873	telomeric_repeat	8.299739	-
chr10	20459861	20459867	telomeric_repeat	8.299739	-
chr10	20459855	20459861	telomeric_repeat	8.299739	-

chr10	20459837	20459843	telomeric_repeat	8.299739	-
chr10	20459831	20459837	telomeric_repeat	8.299739	-
chr10	20459825	20459831	telomeric_repeat	8.299739	-
chr10	20459819	20459825	telomeric_repeat	8.299739	-
chr10	20459813	20459819	telomeric_repeat	8.299739	-

chr10	20459795	20459801	telomeric_repeat	8.299739	-

chr10	20459783	20459789	telomeric_repeat	8.299739	-

# pde7b whole gene
svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr10 20396004 20727068 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'PDE7B gene' -gi 20452966-20467404 -gil '' -hf 20452966 20467404 'Peak' -bgb max -c svist4get_data/A4_p1.cfg -o pde7b_gene.pdf

##############################################################################
# min 4 repeats intergenic
##############################################################################

##########################################################################
cd results

# check annotation of cluster_1 intersecting peaks. Same result if checked in intersecting peaks with 4 tandem repeats
grep -wf <(cat profile/terra_reads_on_rloop_intersecting_peaks.bed | grep cluster_1 | cut -f4) annotate/terra_intersect_annotate_motif_counts.txt | sort -V | cut -f1,2,3,4,8

# check annotation of cluster_2 intersecting peaks with 4 tandem repeats
grep -wf <(cat profile/terra_reads_on_intersecting_peaks_4tandem.bed | grep cluster_2 | cut -f4) annotate/terra_intersect_annotate_motif_counts.txt | sort -V | cut -f1,2,3,4,8

# check annotation of cluster_2 intersecting peaks 
grep -wf <(cat profile/terra_reads_on_intersecting_peaks.bed | grep cluster_2 | cut -f4) annotate/terra_intersect_annotate_motif_counts.txt | sort -V | cut -f1,2,3,4,8

# check number of introns, intergenic, exonx, etc. in certain clusters
grep -wf <(cat profile/terra_reads_on_intersecting_peaks.bed | grep cluster_1 | cut -f4) annotate/terra_intersect_annotate_motif_counts.txt | sort -V | cut -f8 | sed -r 's|(\(.+\))||g' | sort | uniq -c | sort -rn
grep -wf <(cat profile/terra_reads_on_intersecting_peaks.bed | grep cluster_2 | cut -f4) annotate/terra_intersect_annotate_motif_counts.txt | sort -V | cut -f8 | sed -r 's|(\(.+\))||g' | sort | uniq -c | sort -rn
grep -wf <(cat profile/terra_reads_on_intersecting_peaks.bed | grep cluster_3 | cut -f4) annotate/terra_intersect_annotate_motif_counts.txt | sort -V | cut -f8 | sed -r 's|(\(.+\))||g' | sort | uniq -c | sort -rn
grep -wf <(cat profile/terra_reads_on_intersecting_peaks.bed | cut -f4) annotate/terra_intersect_annotate_motif_counts.txt | sort -V | cut -f8 | sed -r 's|(\(.+\)).*||g' | sort | uniq -c | sort -rn

# check number of introns, intergenic, exonx, etc. in intersecting peaks
cat annotate/terra_intersect_annotate_motif_counts.txt | cut -f8 | sed -r 's|(\(.+\)).*||g' | sort | uniq -c | sort -rn

# difference between intersecting peaks before and after excluding any overlaps with blacklisted regions and telomeric regions that I defined
diff <(cat annotate/terra_intersect_annotate_motif_counts.txt | cut -f1,2,3,4,8 | sed -r 's|(\(.+\)).*||g' | grep TTS | sort -V) <(grep -wf <(cat profile/terra_reads_on_intersecting_peaks.bed | cut -f4) annotate/terra_intersect_annotate_motif_counts.txt | sort -V | cut -f1,2,3,4,8 | sed -r 's|(\(.+\)).*||g' | grep TTS | sort -V)
# cluster_1 intergenic peaks that are not at telomeres
chr2:57629491-5762925
chr3:149052455-149052936    chr3:149052400-149053001
chr9:96,852,861-96,853,270

# cluster_2 intergenic peaks that are not at telomeres
# near telomeres
chr3:6115900-6116289        chr3:6,116,043-6,116,248
chr4:153,152,933-153,153,502
# far from telomeres
chr8:83756756-83757133
chr9:27213041-27213460      chr9:27,213,148-27,213,307
chr12:93576804-93577316     chr12:93,576,960-93,577,268
chr13:110603046-110603212   chr13:110,603,037-110,603,142
chr14:86198993-86199141
chr14:88479083-88479622     chr14:88,479,244-88,479,451

# cluster_2 intron peaks that are not at telomeres
chr6:4923775-4924292
chr13:113171841-113172524   chr13:113,172,030-113,172,270   chr13:113,172,103-113,172,221
chr17:62,737,496-62,737,687     chr17:62,737,536-62,737,647
chr17:67741361-67741541
chrX:129056463-129056911    chrX:129,056,527-129,056,736
chrX:141,276,555-141,276,706    chrX:141,276,577-141,276,668

# very interesting terra peak with 2 smaller rloop peaks inside, each with a smaller atrx peak inside just above the telomeric repeats
terra_sense_peak_7647
chr14:69161841-69162690


map2k6
chr11:110422188-110422411





# get peaks from different intron/4tR groups
cd results
cat profile/terra_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed | cut -f4,13 | column -t | less -S

# rloop_intersect_4tandem_intergenic.bed
chr2:124003615-124003831 
chr2:57629423-57629789
chr13:110743970-110744557 # 2 atrx peaks 
chr13:13860789-13861072 
chr5:106176708-106176931
# rloop_intersect_4tandem_intron.bed
chr10:20459715-20459996 
chr16:93789122-93789456
chr12:73772214-73772432
chr17:67741289-67741599
# rloop_intersect_no4tandem_intergenic.bed
chr13:74,530,478-74,530,632 # ACAG repeats
chr6:82909752-82909962 
chr4:64497807-64498124 # 2 atrx peaks
chr19:44660333-44660627 # no atrx peak
# rloop_intersect_no4tandem_intron.bed
chr17:13550672-13551660  # 2 atrx peaks, repeats
chr1:124037229-124037403 # no atrx peak, repeats
chr2:73709121-73709293 # no atrx peak, repeats
##############################################################################
# no consensus repeats intron
##############################################################################

# mink1
# SRR2062968_terra_peak_4135
# chr11:70609756-70609970
svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr11 70609796 70609940 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'MINK1 peak' -bgb max -nts -c svist4get_data/A4_p1.cfg -o mink1_peak.pdf

# arnt
# SRR2062968_terra_peak_16249
# chr3:95479625-95479897
svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr11 95479625 95479897 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'ARNT peak' -bgb max -nts -c svist4get_data/A4_p1.cfg -o arnt_peak.pdf

# vegfc
# SRR2062968_terra_peak_24919
# chr8:54095563-54095883

# svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr8 54095563 54095883 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'VEGFC peak region' -bgb max -nts -c svist4get_data/A4_p1.cfg -o vegfc_peak.pdf

# SRR2075686_rloop_peak_25347
# chr8:54095616-54095823
# chr8:54095647-54095792 zoomed in
svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr8 54095647 54095792 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'VEGFC peak' -bgb max -nts -c svist4get_data/A4_p1.cfg -o vegfc_peak.pdf

# ctnnbl1
# SRR2075686_rloop_peak_15912
# chr2:157811536-157811692
# chr2:157,811,535-157,811,690
svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr2 157811535 157811690 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'CTNNBL1 peak' -bgb max -nts -c svist4get_data/A4_p1.cfg -o ctnnbl1_peak.pdf

# pola1
# chrX:93,562,792-93,563,527
svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chrX 93562700 93563600 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'POLA1 peak' -bgb max -nts -c svist4get_data/A4_p1.cfg -o pola1_peak.pdf

# sh3bgrl2
# chr9:83,559,186-83,559,285
cat $REPEATS | grep chr9 | grep 83559

svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr9 83559186 83559285 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'SH3BGRL2 peak' -gi 83559212-83559254 -gil '' -hf 83559212 83559254 'Telomeric Repeats' -bgb max -nts -c svist4get_data/A4_p1.cfg -o sh3bgrl2_peak.pdf

chr9	83559212	83559218	telomeric_repeat	8.299739	+
chr9	83559218	83559224	telomeric_repeat	8.299739	+
chr9	83559224	83559230	telomeric_repeat	8.299739	+
chr9	83559230	83559236	telomeric_repeat	8.299739	+
chr9	83559236	83559242	telomeric_repeat	8.299739	+
chr9	83559242	83559248	telomeric_repeat	8.299739	+
chr9	83559248	83559254	telomeric_repeat	8.299739	+
chr9	83559079	83559085	telomeric_repeat	8.299739	-


##############################################################################
# no consensus repeats intergenic
##############################################################################

# pcnx peak
# SRR2075686_rloop_peak_5944
# chr12:82013784-82013983

svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr12 82013808 82013958 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'PCNX proximal peak' -bgb max -nts -c svist4get_data/A4_p1.cfg -o pcnx_peak.pdf

# pcnx gene
# chr12:81,858,030-82,025,000
# chr12:82011417-82014339
svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr12 81858030 82025000 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'PCNX gene' -gi 82011417-82014339 -gil '' -hf 82011417 82014339 'Peak' -bgb max -c svist4get_data/A4_p1.cfg -o pcnx_gene.pdf


# pcnx
# chr12:82,044,642-82,044,791
# chr12:82,044,602-82,044,808

# dpy19l1 peak
# chr9:24,542,029-24,542,170
# chr9:24,541,977-24,542,170
svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr9 24542001 24542170 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'DPY19L1 proximal peak' -bgb max -nts -c svist4get_data/A4_p1.cfg -o dpy19l1_peak.pdf

# dpy19l gene
# chr9:24,411,396-24,699,243 - dpy19l1 and 2
# chr9:24,409,627-24,546,424 - dpy19l1 only
# chr9:24537242-24546975 - peak (genes 1 & 2)
# chr9:24540175-24543977 - peak (gene 1)
svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr9 24411396 24699243 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'DPY19L1 and 2 genes' -gi 24537242-24546975 -gil '' -hf 24537242 24546975 'Peak' -bgb max -c svist4get_data/A4_p1.cfg -o dpy19l1_dpy19l2_genes.pdf

svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr9 24409627 24546424 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'DPY19L1 gene' -gi 24540175-24543977 -gil '' -hf 24540175 24543977 'Peak' -bgb max -c svist4get_data/A4_p1.cfg -o dpy19l1_gene.pdf

