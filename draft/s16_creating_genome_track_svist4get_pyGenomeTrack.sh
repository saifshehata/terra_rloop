# pyGenomeTrack

cd /home/saif/terra_rloop_project/genome_browser

# define variables.
TERRA=/home/saif/terra_rloop_project/genome_browser/SRR2062968_terra_AS.bwa_mm10.sorted.rg.rmdup.bw
RLOOP=/home/saif/terra_rloop_project/genome_browser/SRR2075686_rloop.bwa_mm10.sorted.rg.rmdup.bw
ATRX=/home/saif/terra_rloop_project/genome_browser/SRR057566_Atrx.bwa_mm10.sorted.rg.rmdup.bw
REPEATS=/home/saif/terra_rloop_project/genome_browser/terra_overlap_telRepPos.bed
GTF=/home/saif/terra_rloop_project/genome_browser/mm10.ncbiRefSeq.gtf

# make track file.
make_tracks_file --trackFiles $TERRA $RLOOP $ATRX $REPEATS $GTF -o tracks.ini
make_tracks_file --trackFiles SRR2062968_terra_AS.bwa_mm10.sorted.rg.rmdup.bw SRR2075686_rloop.bwa_mm10.sorted.rg.rmdup.bw SRR057566_Atrx.bwa_mm10.sorted.rg.rmdup.bw terra_overlap_telRepPos.bed -o tracks.ini

# create region image.
pyGenomeTracks --tracks tracks.ini --region chr18:65,103,612-65,104,547 --outFileName nice_image.pdf

chr8:34,802,117-35,148,004



################# using svist4get to visualize genome tracks. must create bedgraph files as it does not like bigwig files #################

# load modules.
ml bioinfo-tools 
ml ucsc-utilities
ml BEDTools

# go to file and create
cd /proj/nb_storage/private/terra_rloop_project/results/alignments/

# create bedgraph file.
bedtools genomecov -ibam terra/SRR2062968_terra_AS.bwa_mm10.sorted.rg.rmdup.bam -bg > terra/SRR2062968_terra_AS.bwa_mm10.sorted.rg.rmdup.bg
bedtools genomecov -ibam rloop/SRR2075686_rloop.bwa_mm10.sorted.rg.rmdup.bam -bg > rloop/SRR2075686_rloop.bwa_mm10.sorted.rg.rmdup.bg
bedtools genomecov -ibam law_paper/SRR057566_Atrx.bwa_mm10.sorted.rg.rmdup.bam -bg > law_paper/SRR057566_Atrx.bwa_mm10.sorted.rg.rmdup.bg


### go to local linux file (ubuntu virtual machine). ###
cd /home/saif/terra_rloop_project/genome_browser/genome_tracks

# copy file from server to local linux.
rsync -rvP shehata@rackham.uppmax.uu.se:/crex/proj/nb_storage/private/terra_rloop_project/results/alignments/terra/SRR2062968_terra_AS.bwa_mm10.sorted.rg.rmdup.bg ./
rsync -rvP shehata@rackham.uppmax.uu.se:/crex/proj/nb_storage/private/terra_rloop_project/results/alignments/rloop/SRR2075686_rloop.bwa_mm10.sorted.rg.rmdup.bg ./
rsync -rvP shehata@rackham.uppmax.uu.se:/crex/proj/nb_storage/private/terra_rloop_project/results/alignments/law_paper/SRR057566_Atrx.bwa_mm10.sorted.rg.rmdup.bg ./
# First filter the GTF file by removing XM_ and XR_ accessions to avoid all of them being incorporated in the final figure when using genomic regions (-w) instead of transcript IDs (-t).
cat ./../mm10.ncbiRefSeq.gtf | grep -v XM_ | grep -v XR_ > mm10.ncbiRefSeq_no_XM_XR_accessions.gtf



# define variables and create genome track images.
TERRA=/home/saif/terra_rloop_project/genome_browser/genome_tracks/SRR2062968_terra_AS.bwa_mm10.sorted.rg.rmdup.bg
RLOOP=/home/saif/terra_rloop_project/genome_browser/genome_tracks/SRR2075686_rloop.bwa_mm10.sorted.rg.rmdup.bg
ATRX=/home/saif/terra_rloop_project/genome_browser/genome_tracks/SRR057566_Atrx.bwa_mm10.sorted.rg.rmdup.bg
REPEATS=/home/saif/terra_rloop_project/genome_browser/genome_tracks/terra_overlap_telRepPos.bed
GTF=/home/saif/terra_rloop_project/genome_browser/mm10.ncbiRefSeq.gtf
GTF_FILTERED=/home/saif/terra_rloop_project/genome_browser/genome_tracks/mm10.ncbiRefSeq_no_XM_XR_accessions.gtf

# nedd4l whole gene
# svist4get -gtf $GTF -fa mm10.fa -bg $TERRA $RLOOP $ATRX -t NM_001114386.1 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'Nedd4l gene' -gi 65099500-65109500 -gil '' -hf 65099500 65109500 'Peak' -bgb max -o nedd4l_gene.pdf

svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr18 64887755 65217826 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'Nedd4l gene' -gi 65099500-65109500 -gil '' -hf 65099500 65109500 'Peak' -bgb max -c svist4get_data/A4_p1.cfg -o nedd4l_gene.pdf

# # nedd4l transcript.
# svist4get -gtf $GTF -fa mm10.fa -bg $TERRA $RLOOP $ATRX -t NM_031881.2 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'Nedd4l intron region' -gi 6510100-65108000 -gil '' -hf 65101000 65108000 'Peak' -bgb max -o nedd4l_transcript.pdf

# svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr18 65023526 65217826 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'Nedd4l intron region' -gi 65102000-65107000 -gil '' -hf 65102000 65107000 'Peak' -bgb max -o nedd4l_transcript.pdf


# # nedd4l intron 4 of 29.
# svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr18 65079900 65144000 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'Nedd4l intron region' -gi 65108000-65100000 -gil '' -hf 65102000 65107000 'Peak' -bgb max -o nedd4l_intron.pdf

# nedd4l peak region (only the inner rloop peak region is shown as terra peak is larger). 
# svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr18 65103939 65104162 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'Nedd4l peak region' -gi 65104021-65104087 -gi 65103976-65103982 -gi 65104003-65104009 -gi 65104113-65104119 -gi 65104145-65104151 -gil '' '' '' '' '' -hf 65104021 65104087 'Telomeric Repeats' -hf 65103976 65103982 '' -hf 65104003 65104009 '' -hf 65104113 65104119 '' -hf 65104145 65104151 '' -bgb max -nts -c svist4get_data/A4_p1.cfg -o nedd4l_peak.pdf

svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr18 65103965 65104162 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'Nedd4l peak' -gi 65104021-65104087 -gi 65103976-65103982 -gi 65104003-65104009 -gi 65104113-65104119 -gi 65104145-65104151 -gil '' '' '' '' '' -hf 65104021 65104087 'Telomeric Repeats' -hf 65103976 65103982 '' -hf 65104003 65104009 '' -hf 65104113 65104119 '' -hf 65104145 65104151 '' -bgb max -nts -c svist4get_data/A4_p1.cfg -o nedd4l_peak.pdf

# get telomeric repeat positions from the bed file you created from HOMER annotatePeaks.pl.
cat terra_overlap_telRepPos.bed | grep chr18 | grep 651

chr18   65104145        65104151        telomeric_repeat        8.299739        -
chr18   65104113        65104119        telomeric_repeat        8.299739        -

chr18   65104081        65104087        telomeric_repeat        8.299739        -
chr18   65104075        65104081        telomeric_repeat        8.299739        -
chr18   65104069        65104075        telomeric_repeat        8.299739        -
chr18   65104063        65104069        telomeric_repeat        8.299739        -
chr18   65104057        65104063        telomeric_repeat        8.299739        -
chr18   65104051        65104057        telomeric_repeat        8.299739        -
chr18   65104045        65104051        telomeric_repeat        8.299739        -
chr18   65104039        65104045        telomeric_repeat        8.299739        -
chr18   65104033        65104039        telomeric_repeat        8.299739        -
chr18   65104027        65104033        telomeric_repeat        8.299739        -
chr18   65104021        65104027        telomeric_repeat        8.299739        -

chr18   65104003        65104009        telomeric_repeat        8.299739        -
chr18   65103976        65103982        telomeric_repeat        8.299739        -


Nedd4l
chr18:64887756-65217826
id = NM_001114386.1

Nedd4l
chr18:65023527-65217826
id = NM_031881.2
65102997-65104998

SRR2062968_terra_peak_12499
chr18:65103898-65104242
Score = 427.0
Signal Value: 16.91357
pValue (-log10): 46.87128
qValue (-log10): 42.74085

SRR2075686_rloop_peak_12969
chr18:65103939-65104162
Score = 677.0
Signal Value: 20.07798
pValue (-log10): 72.58802
qValue (-log10): 67.71745


# Tnks gene track zoomOut.
svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr8 34826000 35148000 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'Tnks gene' -gi 35130000-35139000 -gil '' -hf 35130000 35139000 'Peak' -bgb max -c svist4get_data/A4_p1.cfg -o tnks_gene.pdf

# tnks rloop peak region.
svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr8 35134877 35135032 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'Tnks proximal peak' -gi 35134903-35134993 -gil '' -hf 35134903 35134993 'Telomeric Repeats' -bgb max -nts -c svist4get_data/A4_p1.cfg -o tnks_peak.pdf


SRR2075686_rloop_peak_25209
chr8:35134827-35135072
Score = 1382.0
Signal Value: 40.27322
pValue (-log10): 143.72543
qValue (-log10): 138.26538

# get telomeric repeat positions from the bed file you created from HOMER annotatePeaks.pl.
cat terra_overlap_telRepPos.bed | grep chr8 | grep 351

chr8	35134903	35134909	telomeric_repeat	8.299739	+
chr8	35134909	35134915	telomeric_repeat	8.299739	+
chr8	35134915	35134921	telomeric_repeat	8.299739	+
chr8	35134921	35134927	telomeric_repeat	8.299739	+
chr8	35134927	35134933	telomeric_repeat	8.299739	+
chr8	35134933	35134939	telomeric_repeat	8.299739	+
chr8	35134939	35134945	telomeric_repeat	8.299739	+
chr8	35134945	35134951	telomeric_repeat	8.299739	+
chr8	35134951	35134957	telomeric_repeat	8.299739	+
chr8	35134957	35134963	telomeric_repeat	8.299739	+
chr8	35134963	35134969	telomeric_repeat	8.299739	+
chr8	35134969	35134975	telomeric_repeat	8.299739	+
chr8	35134975	35134981	telomeric_repeat	8.299739	+
chr8	35134981	35134987	telomeric_repeat	8.299739	+
chr8	35134987	35134993	telomeric_repeat	8.299739	+

# removing the auto visualization of amino acod track.

sudo nano /home/saif/miniconda3/envs/bioinfo/bin/svist4get

# change  "parameters.config['show_aa_seq_track'] = auto" to "parameters.config['show_aa_seq_track'] = 0"
# save changes and rerun code


##############################################################################
# min 4 repeats intron
##############################################################################

# map2k6
# chr11:110422188-110422411

cat terra_overlap_telRepPos.bed | grep chr11 | grep 110

# map2k6 rloop peak region.
svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr11 110422240 110422390 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'MAP2K6 peak' -gi 110422264-110422330 -gi 110422367-110422373 -gil '' '' -hf 110422264 110422330 'Telomeric Repeats' -hf 110422367 110422373 '' -bgb max -nts -c svist4get_data/A4_p1.cfg -o map2k6_peak.pdf

# bard1
# SRR2075686_rloop_peak_637
# chr1:71079268-71079495

cat terra_overlap_telRepPos.bed | grep chr1 | grep 71079

svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr1 71079316 71079435 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'BARD1 peak' -gi 71079344-71079410 -gil '' -hf 71079344 71079410 'Telomeric Repeats' -hf 71079344 71079410 '' -bgb max -nts -c svist4get_data/A4_p1.cfg -o bard1_peak.pdf

# pde7b
# SRR2075686_rloop_peak_1952
# chr10:20459730-20459984
# chr10:20,459,781-20,459,927

cat terra_overlap_telRepPos.bed | grep chr10 | grep 20459

svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr10 20459771 20459937 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'PDE7B peak' -gi 20459783-20459789 -gi 20459795-20459801 -gi 20459813-20459843 -gi 20459855-20459891 -gi 20459903-20459921 -gil '' '' '' '' '' -hf 20459783 20459789 '' -hf 20459795 20459801 '' -hf 20459813 20459843 '' -hf 20459855 20459891 'Telomeric Repeats' -hf 20459903 20459921 '' -bgb max -nts -c svist4get_data/A4_p1.cfg -o pde7b_peak.pdf

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
svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr10 20396004 20727068 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'PDE7B gene' -gi 20452966-20467404 -gil '' -hf 20452966 20467404 'Peak' -bgb max -c svist4get_data/A4_p1.cfg -o pde7b_gene.pdf

##############################################################################
# min 4 repeats intergenic
##############################################################################

map2k6
chr11:110422188-110422411


##############################################################################
# no consensus repeats intron
##############################################################################

# mink1
# SRR2062968_terra_peak_4135
# chr11:70609756-70609970
svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr11 70609796 70609940 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'MINK1 peak' -bgb max -nts -c svist4get_data/A4_p1.cfg -o mink1_peak.pdf

# arnt
# SRR2062968_terra_peak_16249
# chr3:95479625-95479897
svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr11 95479625 95479897 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'ARNT peak' -bgb max -nts -c svist4get_data/A4_p1.cfg -o arnt_peak.pdf

# vegfc
# SRR2062968_terra_peak_24919
# chr8:54095563-54095883

# svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr8 54095563 54095883 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'VEGFC peak region' -bgb max -nts -c svist4get_data/A4_p1.cfg -o vegfc_peak.pdf

# SRR2075686_rloop_peak_25347
# chr8:54095616-54095823
# chr8:54095647-54095792 zoomed in
svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr8 54095647 54095792 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'VEGFC peak' -bgb max -nts -c svist4get_data/A4_p1.cfg -o vegfc_peak.pdf

# ctnnbl1
# SRR2075686_rloop_peak_15912
# chr2:157811536-157811692
# chr2:157,811,535-157,811,690
svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr2 157811535 157811690 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'CTNNBL1 peak' -bgb max -nts -c svist4get_data/A4_p1.cfg -o ctnnbl1_peak.pdf

# pola1
# chrX:93,562,792-93,563,527
svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chrX 93562700 93563600 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'POLA1 peak' -bgb max -nts -c svist4get_data/A4_p1.cfg -o pola1_peak.pdf

# sh3bgrl2
# chr9:83,559,186-83,559,285
cat terra_overlap_telRepPos.bed | grep chr9 | grep 83559

svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr9 83559186 83559285 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'SH3BGRL2 peak' -gi 83559212-83559254 -gil '' -hf 83559212 83559254 'Telomeric Repeats' -bgb max -nts -c svist4get_data/A4_p1.cfg -o sh3bgrl2_peak.pdf

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

svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr12 82013808 82013958 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'PCNX proximal peak' -bgb max -nts -c svist4get_data/A4_p1.cfg -o pcnx_peak.pdf

# pcnx gene
# chr12:81,858,030-82,025,000
# chr12:82011417-82014339
svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr12 81858030 82025000 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'PCNX gene' -gi 82011417-82014339 -gil '' -hf 82011417 82014339 'Peak' -bgb max -c svist4get_data/A4_p1.cfg -o pcnx_gene.pdf


# pcnx
# chr12:82,044,642-82,044,791
# chr12:82,044,602-82,044,808

# dpy19l1 peak
# chr9:24,542,029-24,542,170
# chr9:24,541,977-24,542,170
svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr9 24542001 24542170 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'DPY19L1 proximal peak' -bgb max -nts -c svist4get_data/A4_p1.cfg -o dpy19l1_peak.pdf

# dpy19l gene
# chr9:24,411,396-24,699,243 - dpy19l1 and 2
# chr9:24,409,627-24,546,424 - dpy19l1 only
# chr9:24537242-24546975 - peak (genes 1 & 2)
# chr9:24540175-24543977 - peak (gene 1)
svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr9 24411396 24699243 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'DPY19L1 and 2 genes' -gi 24537242-24546975 -gil '' -hf 24537242 24546975 'Peak' -bgb max -c svist4get_data/A4_p1.cfg -o dpy19l1_dpy19l2_genes.pdf

svist4get -gtf $GTF_FILTERED -fa mm10.fa -bg $TERRA $RLOOP $ATRX -w chr9 24409627 24546424 -bl 'TERRA CHIRT-Seq' 'R-Loop DRIP-Seq' 'ATRX ChIP-Seq' -it 'DPY19L1 gene' -gi 24540175-24543977 -gil '' -hf 24540175 24543977 'Peak' -bgb max -c svist4get_data/A4_p1.cfg -o dpy19l1_gene.pdf

