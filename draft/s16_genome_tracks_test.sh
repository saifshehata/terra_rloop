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


# intergenic 4 tandem repeats
svist4get -gtf $GTF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr5 93289561 93289789 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'Intergenic peak with >= 4 tandem repeats' -gi 93289638-93289710 -gil '' -hf 93289638 93289710 'Telomeric Repeats' -bgb max -nts -c svist4get_data/A4_p1.cfg -o chr5_93289561_932897891_4tR_intergenic.pdf

# intergenic no4 tandem repeats
svist4get -gtf $GTF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr13 119759111 119759329 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'Intergenic peak with < 4 tandem repeats' -o chr13_119759111_119759329_no4tR_intergenic.pdf


# intron 4 tandem repeats
svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr1 71079316 71079435 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'Intron (BARD1) peak with >= 4 tandem repeats' -gi 71079344-71079410 -gil '' -hf 71079344 71079410 'Telomeric Repeats' -hf 71079344 71079410 '' -bgb max -nts -c svist4get_data/A4_p1.cfg -o bard1_peak.pdf

svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr11 110422240 110422390 -bl 'TERRA' 'R-Loop' 'ATRX' -it ' Intron (MAP2K6) peak with >= 4 tandem repeats' -gi 110422264-110422330 -gi 110422367-110422373 -gil '' '' -hf 110422264 110422330 'Telomeric Repeats' -hf 110422367 110422373 '' -bgb max -nts -c svist4get_data/A4_p1.cfg -o map2k6_peak.pdf


# intron no 4 tandem repeats
svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr8 54095647 54095792 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'Intron (VEGFC) peak with < 4 tandem repeats' -bgb max -nts -c svist4get_data/A4_p1.cfg -o vegfc_peak.pdf

svist4get -gtf $GFF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr11 70609796 70609940 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'Intron (MINK1) peak with < 4 tandem repeats' -bgb max -nts -c svist4get_data/A4_p1.cfg -o mink1_peak.pdf





# create genome tracks
# svist4get -gtf $GTF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr2 124003615 124003831 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'Intergenic >= 4 tandem repeats' -o 1.png

# rloop_intersect_4tandem_intergenic.bed
# chr2:124003615-124003831 
# chr2:57629423-57629789
# chr13:110743970-110744557 # 2 atrx peaks 
# chr13:13860789-13861072 
# chr5:106176708-106176931



# svist4get -gtf $GTF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr10 20459715 20459996 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'Intron >= 4 tandem repeats' -o 2.png

# rloop_intersect_4tandem_intron.bed
# chr10:20459715-20459996 
# chr16:93789122-93789456
# chr12:73772214-73772432
# chr17:67741289-67741599



# svist4get -gtf $GTF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr13 74530478 74530632 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'Intergenic < 4 tandem repeats' -o 3.png

# rloop_intersect_no4tandem_intergenic.bed
# chr13:74530478-74530632 # ACAG repeats
# chr6:82909752-82909962 
# chr4:64497807-64498124 # 2 atrx peaks
# chr19:44660333-44660627 # no atrx peak



# svist4get -gtf $GTF_FILTERED -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w chr17 13550672 13551660 -bl 'TERRA' 'R-Loop' 'ATRX' -it 'Intron < 4 tandem repeats' -o 4.png

# rloop_intersect_no4tandem_intron.bed
# chr17:13550672-13551660  # 2 atrx peaks, repeats
# chr1:124037229-124037403 # no atrx peak, repeats
# chr2:73709121-73709293 # no atrx peak, repeats

