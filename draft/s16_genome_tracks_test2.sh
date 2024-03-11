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
wget -P $genome --timestamping http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz
gunzip -c $genome/mm10.refGene.gtf.gz > $genome/mm10.refGene.gtf

# define variables
TERRA_BG=$mapped_reads/SRR2062968_pe_sort_rmdup.bg
RLOOP_BG=$mapped_reads/SRR2075686_se_sort_rmdup.bg
ATRX_BG=$mapped_reads/SRR057567_se_sort_rmdup.bg

REPEATS=$annotate/terra_intersect_motif_positions.bed
GTF=$genome/mm10.refGene.gtf
REF=$genome/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna



# get 5 genomic ranges from the each group and create genomic track for each 
cat $profile/terra_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed | cut -f4,13 | sed 1d | grep _4tandem_intergenic | head -40 | tail -5 | awk '{print $1, $2, $1}' | sed 's|:| |' | sed 's|-| |' | sed 's|:|_|' | sed 's|-|_|' | parallel --colsep ' ' "svist4get -gtf $GTF -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w {1} {2} {3} -bl 'TERRA' 'R-Loop' 'ATRX' -it '{4.}' -bgb max -nts -c $genome_tracks/svist4get_data/A4_p1.cfg -o $genome_tracks/{5}.pdf"

echo 'done _4tandem_intergenic'

cat $profile/terra_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed | cut -f4,13 | sed 1d | grep _4tandem_intron | head -30 | tail -5 | awk '{print $1, $2, $1}' | sed 's|:| |' | sed 's|-| |' | sed 's|:|_|' | sed 's|-|_|' | parallel --colsep ' ' "svist4get -gtf $GTF -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w {1} {2} {3} -bl 'TERRA' 'R-Loop' 'ATRX' -it '{4.}' -bgb max -nts -c $genome_tracks/svist4get_data/A4_p1.cfg -o $genome_tracks/{5}.pdf"

echo 'done _4tandem_intron'

cat $profile/terra_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed | cut -f4,13 | sed 1d | grep _no4tandem_intergenic | head -60 | tail -5 | awk '{print $1, $2, $1}' | sed 's|:| |' | sed 's|-| |' | sed 's|:|_|' | sed 's|-|_|' | parallel --colsep ' ' "svist4get -gtf $GTF -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w {1} {2} {3} -bl 'TERRA' 'R-Loop' 'ATRX' -it '{4.}' -bgb max -nts -c $genome_tracks/svist4get_data/A4_p1.cfg -o $genome_tracks/{5}.pdf"

echo 'done _no4tandem_intergenic'

cat $profile/terra_reads_on_rloop_intersecting_peaks_4tR_no4tR_intergenic_intron.bed | cut -f4,13 | sed 1d | grep _no4tandem_intron | head -60 | tail -5 | awk '{print $1, $2, $1}' | sed 's|:| |' | sed 's|-| |' | sed 's|:|_|' | sed 's|-|_|' | parallel --colsep ' ' "svist4get -gtf $GTF -fa $REF -bg $TERRA_BG $RLOOP_BG $ATRX_BG -w {1} {2} {3} -bl 'TERRA' 'R-Loop' 'ATRX' -it '{4.}' -bgb max -nts -c $genome_tracks/svist4get_data/A4_p1.cfg -o $genome_tracks/{5}.pdf"

echo 'done _no4tandem_intron'
