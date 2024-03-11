#! /bin/bash
#$ -N meme_common_motifs
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 6

# go to directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
peaks=$project_dir/results/peaks
overlaps=$project_dir/results/overlaps
annotate=$project_dir/results/annotate
meme=$project_dir/results/meme

mkdir -p $meme
cd $meme


################################################################################################ 
# Filter intersecting peaks by the presence of >= 4 tandem telomeric repeats 
################################################################################################ 

# define variables for intersecting peaks
TERRA_INTERSECT=$overlaps/terra_intersect.narrowPeak
RLOOP_INTERSECT=$overlaps/rloop_intersect.narrowPeak

# define variables for files to be created
TERRA_4t_PEAK_NAMES=$overlaps/terra_intersect_4tandem_peak_names.txt
RLOOP_4t_PEAK_NAMES=$overlaps/rloop_intersect_4tandem_peak_names.txt

# get names of peaks that contain at least 4 tandem repeats, for both terra and rloop
cat $annotate/terra_intersect_annotate_motif_positions.txt | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | cut -d" " -f1 > $TERRA_4t_PEAK_NAMES
cat $annotate/rloop_intersect_annotate_motif_positions.txt | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | cut -d" " -f1 > $RLOOP_4t_PEAK_NAMES

# for both terra and rloop, filter intersecting peaks for those that contain, or lack, at least 4 tandem repeats
### containing 4 tandem repeats ###
grep -wf $TERRA_4t_PEAK_NAMES $TERRA_INTERSECT > $overlaps/terra_intersect_4tandem.narrowPeak
grep -wf $RLOOP_4t_PEAK_NAMES $RLOOP_INTERSECT > $overlaps/rloop_intersect_4tandem.narrowPeak

### lacking 4 tandem repeats ###
grep -v -wf $TERRA_4t_PEAK_NAMES $TERRA_INTERSECT > $overlaps/terra_intersect_no4tandem.narrowPeak
grep -v -wf $RLOOP_4t_PEAK_NAMES $RLOOP_INTERSECT > $overlaps/rloop_intersect_no4tandem.narrowPeak

# # example of another way to get the same result using join (faster than grep?), although the final sorting is slightly different due to sort -V
# join -1 1 -2 4 -o 2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10 <(cat $annotate/terra4tandem_peak_names.txt | sort -k1,1) <(cat $TERRA_INTERSECT | sort -k4,4) | sort -k1V -k2n | sed 's/ /\t/g'  > $overlaps/terra_intersect_4tandem.narrowPeak

################################################################################################ 
# now find common motifs within the intersecting peaks with or without 4 tandem telomeric repeats
# using the meme suite
################################################################################################ 

# first get the genomic sequences from the peak coordinates.
#bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF>
REF=$project_dir/raw_data/genome/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna
bedtools getfasta -fi $REF -bed $overlaps/terra_intersect_4tandem.narrowPeak -name -fo $meme/terra_intersect_4tandem.fasta
bedtools getfasta -fi $REF -bed $overlaps/rloop_intersect_4tandem.narrowPeak -name -fo $meme/rloop_intersect_4tandem.fasta

bedtools getfasta -fi $REF -bed $overlaps/terra_intersect_no4tandem.narrowPeak -name -fo $meme/terra_intersect_no4tandem.fasta
bedtools getfasta -fi $REF -bed $overlaps/rloop_intersect_no4tandem.narrowPeak -name -fo $meme/rloop_intersect_no4tandem.fasta

# then find the common motifs in those sequences.
#meme  <dataset> [optional arguments] -oc <output dir[will replace existing directory]> -revcomp -nmotifs <nmotifs> -evt <ev>
meme $meme/terra_intersect_4tandem.fasta -dna -revcomp -nmotifs 15 -oc $meme/meme_terra_intersect_4tandem
meme $meme/rloop_intersect_4tandem.fasta -dna -revcomp -nmotifs 15 -oc $meme/meme_rloop_intersect_4tandem

meme $meme/terra_intersect_no4tandem.fasta -dna -revcomp -nmotifs 15 -oc $meme/meme_terra_intersect_no4tandem
meme $meme/rloop_intersect_no4tandem.fasta -dna -revcomp -nmotifs 15 -oc $meme/meme_rloop_intersect_no4tandem



#################################

# ########################################################################################################################################################
# # EXTRA USEFUL COMMANDS TO GET INTRON/INTERGENIC COUNT AND GENE COUNT FROM ANNOTATION AND GTF FILES, RESPECTIVELY
# ########################################################################################################################################################

# # go to directory
# cd /proj/nb_storage/private/terra_rloop_project/results/peaks/annotation/

# # get nr of introns/exons
# cat terra_overlap_narrowPeak_ann_nrRepeats.txt | cut -f8 | sed 1d | sed -r 's/(\(.+\).*)//' | sort | uniq -c | sort 

# # 2 non-coding 
# # 4 3' UTR 
# # 6 exon 
# # 12 NA
# # 14 promoter-TSS 
# # 15 TTS 
# # 309 intron 
# # 448 Intergenic


# cat terra_overlap_min4repeats_ann.txt | cut -f8 | sed 1d | sed -r 's/(\(.+\).*)//' | sort | uniq -c | sort 

# # 1 TTS 
# # 2 3' UTR 
# # 63 intron 
# # 108 intergenic

# cat terra_overlap_noConsRepeats_ann.txt | cut -f8 | sed 1d | sed -r 's/(\(.+\).*)//' | sort | uniq -c | sort 

# # 2 3' UTR 
# # 2 non-coding 
# # 6 exon 
# # 12 NA
# # 14 TTS 
# # 14 promoter-TSS 
# # 246 intron 
# # 340 Intergenic
# # 
# # get nr genes
# cat mm10.ncbiRefSeq.gtf | cut -f9 | cut -f10 -d" "| sort | uniq | sort | grep -v -e "1" | sed -r 's/"|;//g' | head
# ########################################################################################################################################################
# # GET GO TERMS AND GENOME ONTOLOGY FOR OVERLAPPING PEAKS +/- 4 CONSECUTIVE REPEATS
# ########################################################################################################################################################

# # go to directory
# cd /proj/nb_storage/private/terra_rloop_project/results/peaks/annotation/

# # define variables.
# TERRA_OVERLAP_4CR=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_min4repeats.narrowPeak
# TERRA_OVERLAP_NO4CR=/proj/nb_storage/private/terra_rloop_project/results/peaks/overlap/terra_overlap_noConsRepeats.narrowPeak

# # for terra_overlap_min4repeats.narrowPeak
# annotatePeaks.pl $TERRA_OVERLAP_4CR mm10 -size given -go geneOnt_terra_overlap_min4repeats -genomeOntology genomeOnt_terra_overlap_min4repeats -cpu 10 > terra_overlap_min4repeats_ann.txt

# # for terra_overlap_noConsRepeats.narrowPeak
# annotatePeaks.pl $TERRA_OVERLAP_NO4CR mm10 -size given -go geneOnt_terra_overlap_noConsRepeats -genomeOntology genomeOnt_terra_overlap_noConsRepeats -cpu 10 > terra_overlap_noConsRepeats_ann.txt

