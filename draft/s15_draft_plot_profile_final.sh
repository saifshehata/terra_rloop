

# plot_profile_final.sh

# load modules.
module load bioinfo-tools
#module load HOMER
module load deepTools

######################################################################################################################################################
# plot profile for terra reads on terra OVERLAPPING_4CR (INTERGENIC & INTRON) VS OVERLAPPING_NO4CR (INTERGENIC & INTRON)
######################################################################################################################################################

# go to directory.
cd /proj/nb_storage/private/terra_rloop_project/results/figures/plot_profile/

# define variables. copied from above.
TERRA_AS_BW=/proj/nb_storage/private/terra_rloop_project/results/alignments/terra/SRR2062968_terra_AS.bwa_mm10.sorted.rg.rmdup.bw
TERRA_OVERLAP_4CR_INTERGENIC=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_min4repeats_intergenic.narrowPeak
TERRA_OVERLAP_4CR_INTRON=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_min4repeats_intron.narrowPeak
TERRA_OVERLAP_NO4CR_INTERGENIC=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_noConsRepeats_intergenic.narrowPeak
TERRA_OVERLAP_NO4CR_INTRON=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_noConsRepeats_intron.narrowPeak

#--regionsFileName $TERRA_OVERLAP_4CR_INTERGENIC $TERRA_OVERLAP_4CR_INTRON $TERRA_OVERLAP_NO4CR_INTERGENIC $TERRA_OVERLAP_NO4CR_INTRON \
computeMatrix scale-regions --scoreFileName $TERRA_AS_BW \
--regionsFileName $TERRA_OVERLAP_4CR_INTRON $TERRA_OVERLAP_4CR_INTERGENIC $TERRA_OVERLAP_NO4CR_INTRON $TERRA_OVERLAP_NO4CR_INTERGENIC \
--startLabel "peak start" \
--endLabel "peak end" \
--samplesLabel "Terra reads" \
--upstream 1000 \
--regionBodyLength 1000 \
--downstream 1000 \
--binSize 10 \
--outFileName terraOnPeaks_Overlap_4crVSno4cr_intergenicVSintron.matrix.gz \
--outFileNameMatrix terraOnPeaks_Overlap_4crVSno4cr_intergenicVSintron.tab \
--numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile terraOnPeaks_Overlap_4crVSno4cr_intergenicVSintron.matrix.gz \
--outFileName terraOnPeaks_Overlap_4crVSno4cr_intergenicVSintron.pdf


######################################################################################################################################################
# plot profile for ATRX reads on terra OVERLAPPING_4CR (INTERGENIC & INTRON) VS OVERLAPPING_NO4CR (INTERGENIC & INTRON)
######################################################################################################################################################

# go to directory.
cd /proj/nb_storage/private/terra_rloop_project/results/figures/plot_profile/

# define variables. copied from above.
ATRX_LAW1=/proj/nb_storage/private/terra_rloop_project/results/alignments/law_paper/SRR057566_Atrx.bwa_mm10.sorted.rg.rmdup.bw
ATRX_LAW2=/proj/nb_storage/private/terra_rloop_project/results/alignments/law_paper/SRR057567_Atrx.bwa_mm10.sorted.rg.rmdup.bw
TERRA_OVERLAP_4CR_INTERGENIC=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_min4repeats_intergenic.narrowPeak
TERRA_OVERLAP_4CR_INTRON=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_min4repeats_intron.narrowPeak
TERRA_OVERLAP_NO4CR_INTERGENIC=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_noConsRepeats_intergenic.narrowPeak
TERRA_OVERLAP_NO4CR_INTRON=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_noConsRepeats_intron.narrowPeak

#--regionsFileName $TERRA_OVERLAP_4CR_INTERGENIC $TERRA_OVERLAP_4CR_INTRON $TERRA_OVERLAP_NO4CR_INTERGENIC $TERRA_OVERLAP_NO4CR_INTRON \
computeMatrix scale-regions --scoreFileName $ATRX_LAW2 \
--regionsFileName $TERRA_OVERLAP_4CR_INTRON $TERRA_OVERLAP_4CR_INTERGENIC $TERRA_OVERLAP_NO4CR_INTRON $TERRA_OVERLAP_NO4CR_INTERGENIC \
--startLabel "peak start" \
--endLabel "peak end" \
--samplesLabel "ATRX-2 reads" \
--upstream 1000 \
--regionBodyLength 1000 \
--downstream 1000 \
--binSize 10 \
--outFileName atrxOnPeaks_Overlap_4crVSno4cr_intergenicVSintron.matrix.gz \
--outFileNameMatrix atrxOnPeaks_Overlap_4crVSno4cr_intergenicVSintron.tab \
--numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile atrxOnPeaks_Overlap_4crVSno4cr_intergenicVSintron.matrix.gz \
--outFileName atrxOnPeaks_Overlap_4crVSno4cr_intergenicVSintron.pdf


######################################################################################################################################################
# plot profile for RLOOP reads on terra OVERLAPPING_4CR (INTERGENIC & INTRON) VS OVERLAPPING_NO4CR (INTERGENIC & INTRON)
######################################################################################################################################################

# go to directory.
cd /proj/nb_storage/private/terra_rloop_project/results/figures/plot_profile/

# define variables. copied from above.
RLOOP_BW=/proj/nb_storage/private/terra_rloop_project/results/alignments/rloop/SRR2075686_rloop.bwa_mm10.sorted.rg.rmdup.bw
TERRA_OVERLAP_4CR_INTERGENIC=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_min4repeats_intergenic.narrowPeak
TERRA_OVERLAP_4CR_INTRON=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_min4repeats_intron.narrowPeak
TERRA_OVERLAP_NO4CR_INTERGENIC=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_noConsRepeats_intergenic.narrowPeak
TERRA_OVERLAP_NO4CR_INTRON=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_noConsRepeats_intron.narrowPeak

#--regionsFileName $TERRA_OVERLAP_4CR_INTERGENIC $TERRA_OVERLAP_4CR_INTRON $TERRA_OVERLAP_NO4CR_INTERGENIC $TERRA_OVERLAP_NO4CR_INTRON \
computeMatrix scale-regions --scoreFileName $RLOOP_BW \
--regionsFileName $TERRA_OVERLAP_4CR_INTRON $TERRA_OVERLAP_4CR_INTERGENIC $TERRA_OVERLAP_NO4CR_INTRON $TERRA_OVERLAP_NO4CR_INTERGENIC \
--startLabel "peak start" \
--endLabel "peak end" \
--samplesLabel "R-loop reads" \
--upstream 1000 \
--regionBodyLength 1000 \
--downstream 1000 \
--binSize 10 \
--outFileName rloopOnPeaks_Overlap_4crVSno4cr_intergenicVSintron.matrix.gz \
--outFileNameMatrix rloopOnPeaks_Overlap_4crVSno4cr_intergenicVSintron.tab \
--numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile rloopOnPeaks_Overlap_4crVSno4cr_intergenicVSintron.matrix.gz \
--outFileName rloopOnPeaks_Overlap_4crVSno4cr_intergenicVSintron.pdf


######################################################################################################################################################
# plot profile for RLOOP reads on terra OVERLAPPING_INTRON peaks (4CR VS NO4CR)
######################################################################################################################################################

# go to directory.
cd /proj/nb_storage/private/terra_rloop_project/results/figures/plot_profile/

# define variables. copied from above.
RLOOP_BW=/proj/nb_storage/private/terra_rloop_project/results/alignments/rloop/SRR2075686_rloop.bwa_mm10.sorted.rg.rmdup.bw
TERRA_OVERLAP_4CR_INTRON=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_min4repeats_intron.narrowPeak
TERRA_OVERLAP_NO4CR_INTRON=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_noConsRepeats_intron.narrowPeak

computeMatrix scale-regions --scoreFileName $RLOOP_BW \
--regionsFileName $TERRA_OVERLAP_4CR_INTRON $TERRA_OVERLAP_NO4CR_INTRON \
--startLabel "peak start" \
--endLabel "peak end" \
--samplesLabel "R-loop reads" \
--upstream 1000 \
--regionBodyLength 1000 \
--downstream 1000 \
--binSize 10 \
--outFileName rloopOnPeaks_Overlap_intron_4crVSno4cr.matrix.gz \
--outFileNameMatrix rloopOnPeaks_Overlap_intron_4crVSno4cr.tab \
--numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile rloopOnPeaks_Overlap_intron_4crVSno4cr.matrix.gz \
--outFileName rloopOnPeaks_Overlap_intron_4crVSno4cr.pdf


###################### previous profiles to edit and combine ################################

### test ###
      2 Isl1
      2 Cr2
      2 Asmt
      2 Adcyap1



######################################################################################################################################################
# plot profile for terra reads on terra OVERLAPPING_4CR_INTERGENIC VS OVERLAPPING_4CR_INTRON
######################################################################################################################################################

# go to directory.
cd /proj/nb_storage/private/terra_rloop_project/results/figures/plot_profile/

# define variables. copied from above.
TERRA_AS_BW=/proj/nb_storage/private/terra_rloop_project/results/alignments/terra/SRR2062968_terra_AS.bwa_mm10.sorted.rg.rmdup.bw
TERRA_OVERLAP_4CR_INTERGENIC=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_min4repeats_intergenic.narrowPeak
TERRA_OVERLAP_4CR_INTRON=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_min4repeats_intron.narrowPeak
TERRA_OVERLAP_NO4CR_INTRON=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_noConsRepeats_intron.narrowPeak

computeMatrix scale-regions --scoreFileName $TERRA_AS_BW \
--regionsFileName $TERRA_OVERLAP_4CR_INTERGENIC $TERRA_OVERLAP_4CR_INTRON \
--startLabel "peak start" \
--endLabel "peak end" \
--samplesLabel "Terra reads" \
--upstream 1000 \
--regionBodyLength 1000 \
--downstream 1000 \
--binSize 10 \
--outFileName terraOnPeaks_Overlap4crIntergenic_vs_Overlap4crIntron.matrix.gz \
--outFileNameMatrix terraOnPeaks_Overlap4crIntergenic_vs_Overlap4crIntron.tab \
--numberOfProcessors max

# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile terraOnPeaks_Overlap4crIntergenic_vs_Overlap4crIntron.matrix.gz \
--outFileName terraOnPeaks_Overlap4crIntergenic_vs_Overlap4crIntron.pdf


######################################################################################################################################################
# plot profile for terra reads on terra OVERLAPPING_NO4CR_INTERGENIC VS OVERLAPPING_NO4CR_INTRON
######################################################################################################################################################

# go to directory.
cd /proj/nb_storage/private/terra_rloop_project/results/figures/plot_profile/

# define variables. copied from above.
TERRA_AS_BW=/proj/nb_storage/private/terra_rloop_project/results/alignments/terra/SRR2062968_terra_AS.bwa_mm10.sorted.rg.rmdup.bw
TERRA_OVERLAP_NO4CR_INTERGENIC=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_noConsRepeats_intergenic.narrowPeak
TERRA_OVERLAP_NO4CR_INTRON=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_noConsRepeats_intron.narrowPeak

computeMatrix scale-regions --scoreFileName $TERRA_AS_BW \
--regionsFileName $TERRA_OVERLAP_NO4CR_INTERGENIC $TERRA_OVERLAP_NO4CR_INTRON \
--startLabel "peak start" \
--endLabel "peak end" \
--samplesLabel "Terra reads" \
--upstream 1000 \
--regionBodyLength 1000 \
--downstream 1000 \
--binSize 10 \
--outFileName terraOnPeaks_OverlapNo4crIntergenic_vs_OverlapNo4crIntron.matrix.gz \
--outFileNameMatrix terraOnPeaks_OverlapNo4crIntergenic_vs_OverlapNo4crIntron.tab \
--numberOfProcessors max


# create a profile plot for scores over sets of genomic regions.
plotProfile --matrixFile terraOnPeaks_OverlapNo4crIntergenic_vs_OverlapNo4crIntron.matrix.gz \
--outFileName terraOnPeaks_OverlapNo4crIntergenic_vs_OverlapNo4crIntron.pdf


###
