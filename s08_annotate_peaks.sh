#! /bin/bash
#$ -N annotate_peaks
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 5

# go to directory
main_dir=/cluster/khiom/sshehata001/
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
peaks=$project_dir/results/peaks
overlaps=$project_dir/results/overlaps
annotate=$project_dir/results/annotate
mkdir -p $annotate
cd $annotate


# add the mm10 genome to be used for annotating, otherwise this error will appear
# !!!!Genome mm10 not found in /cluster/khiom/sshehata001/miniconda3/envs/ngs/share/homer/.//config.txt

#         To check if is available, run "perl /cluster/khiom/sshehata001/miniconda3/envs/ngs/share/homer/.//configureHomer.pl -list"
#         If so, add it by typing "perl /cluster/khiom/sshehata001/miniconda3/envs/ngs/share/homer/.//configureHomer.pl -install mm10"
# usage: perl <main_dir>/miniconda3/envs/<env_name>/share/homer/.//configureHomer.pl -install mm10
perl $main_dir/miniconda3/envs/ngs/share/homer/.//configureHomer.pl -install mm10

# define some variables.
TERRA_PEAKS=$peaks/terra_sense_peaks_filtered.narrowPeak
TERRA_INTERSECT=$overlaps/terra_intersect.narrowPeak
TERRA_NO_INTERSECT=$overlaps/terra_no_intersect.narrowPeak

RLOOP_PEAKS=$peaks/rloop_peaks_filtered.narrowPeak
RLOOP_INTERSECT=$overlaps/rloop_intersect.narrowPeak
RLOOP_NO_INTERSECT=$overlaps/rloop_no_intersect.narrowPeak

# create a telomeric repeat (TTAGGG) motif matrix file using HOMER's seq2profile.pl script. This takes into consideration both TTAGGG and the reverse complement CCCTAA
# usage: seq2profile.pl <consensus> [# mismatches] [name] > output.motif
seq2profile.pl TTAGGG 0 telomeric_repeat > $annotate/telomeric.repeat

################################################################################################ 
# ANNOTATE PEAKS TWO TIMES: ONCE BY REPORTING TOTAL REPEAT COUNT, AND ONCE BY REPORTING ALL REPEAT POSITIONS
################################################################################################ 

# annotate peaks using HOMER's annotatePeaks.pl script. only report number of motifs/repeats per peak using -nmotifs.
annotatePeaks.pl $TERRA_PEAKS mm10 -size given -m telomeric.repeat -nmotifs -cpu 10 > terra_peaks_annotate_motif_counts.txt
annotatePeaks.pl $RLOOP_PEAKS mm10 -size given -m telomeric.repeat -nmotifs -cpu 10 > rloop_peaks_annotate_motif_counts.txt
annotatePeaks.pl $TERRA_INTERSECT mm10 -size given -m telomeric.repeat -nmotifs -cpu 10 > terra_intersect_annotate_motif_counts.txt
annotatePeaks.pl $RLOOP_INTERSECT mm10 -size given -m telomeric.repeat -nmotifs -cpu 10 > rloop_intersect_annotate_motif_counts.txt
annotatePeaks.pl $TERRA_NO_INTERSECT mm10 -size given -m telomeric.repeat -nmotifs -cpu 10 > terra_no_intersect_annotate_motif_counts.txt
annotatePeaks.pl $RLOOP_NO_INTERSECT mm10 -size given -m telomeric.repeat -nmotifs -cpu 10 > rloop_no_intersect_annotate_motif_counts.txt

# annotate peaks. show all motifs/repeats per peak, with distance of each repeat to peak centre by NOT using -nmotifs, and extract telomeric repeat positions in .bed files.
annotatePeaks.pl $TERRA_PEAKS mm10 -size given -m telomeric.repeat -mbed terra_peaks_motif_positions.bed -cpu 10 > terra_peaks_annotate_motif_positions.txt
annotatePeaks.pl $RLOOP_PEAKS mm10 -size given -m telomeric.repeat -mbed rloop_peaks_motif_positions.bed -cpu 10 > rloop_peaks_annotate_motif_positions.txt
annotatePeaks.pl $TERRA_INTERSECT mm10 -size given -m telomeric.repeat -mbed terra_intersect_motif_positions.bed -cpu 10 > terra_intersect_annotate_motif_positions.txt
annotatePeaks.pl $RLOOP_INTERSECT mm10 -size given -m telomeric.repeat -mbed rloop_intersect_motif_positions.bed -cpu 10 > rloop_intersect_annotate_motif_positions.txt
annotatePeaks.pl $TERRA_NO_INTERSECT mm10 -size given -m telomeric.repeat -mbed terra_no_intersect_motif_positions.bed -cpu 10 > terra_no_intersect_annotate_motif_positions.txt
annotatePeaks.pl $RLOOP_NO_INTERSECT mm10 -size given -m telomeric.repeat -mbed rloop_no_intersect_motif_positions.bed -cpu 10 > rloop_no_intersect_annotate_motif_positions.txt


# example of how t0 annotate and also add Gene/Genome Ontology.
#annotatePeaks.pl $RLOOP_INTERSECT mm10 -size given -m telomeric.repeat -nmotifs -mbed rloop_overlap_repPos.bed -go rloop_overlap_geneOnt -genomeOntology rloop_overlap_genomeOnt -cpu 10 > rloop_intersect_ann.txt


# ################################################################################################ 
# # Filter intersecting peaks by the presence of >= 4 tandem telomeric repeats 
# ################################################################################################ 

# # define variables
# ### annotated intersecting peaks ###
# TERRA_INTERSECT_ANN_POSITIONS=$annotate/terra_intersect_annotate_motif_positions.txt
# RLOOP_INTERSECT_ANN_POSITIONS=$annotate/rloop_intersect_annotate_motif_positions.txt

# # get names of peaks that contain at least 4 tandem repeats, for both terra and rloop
# cat $TERRA_INTERSECT_ANN_POSITIONS | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | cut -d" " -f1 > $annotate/terra4tandem_peak_names.txt
# cat $RLOOP_INTERSECT_ANN_POSITIONS | cut -f1,22 | grep -E "(TTAGGG|CCCTAA)" | sed -r 's|(\(.{13}\))||g'| tr "," "\t" | awk '{for (i=2;i<NF;i++) $i=$(i+1)-$i}1' | awk 'NF{NF--};1' | grep -E "( (6|7)){3,}" | cut -d" " -f1 > $annotate/rloop4tandem_peak_names.txt

# # for both terra and rloop, filter intersecting peaks for those that contain at least 4 tandem repeats
# grep -wf $annotate/terra4tandem_peak_names.txt $TERRA_INTERSECT > $overlaps/terra_intersect_4tandem.narrowPeak
# grep -wf $annotate/rloop4tandem_peak_names.txt $RLOOP_INTERSECT > $overlaps/rloop_intersect_4tandem.narrowPeak

# # # example of another way to get the same result using join (faster than grep?), although the final sorting is slightly different due to sort -V
# # join -1 1 -2 4 -o 2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10 <(cat $annotate/terra4tandem_peak_names.txt | sort -k1,1) <(cat $TERRA_INTERSECT | sort -k4,4) | sort -k1V -k2n | sed 's/ /\t/g'  > $overlaps/terra_intersect_4tandem.narrowPeak

# # for both terra and rloop, get peaks lacking 4 tandem repeats
# grep -v -wf $annotate/terra4tandem_peak_names.txt $TERRA_INTERSECT > $overlaps/terra_intersect_no4tandem.narrowPeak
# grep -v -wf $annotate/rloop4tandem_peak_names.txt $RLOOP_INTERSECT > $overlaps/rloop_intersect_no4tandem.narrowPeak

# # remove not-needed files
# rm $annotate/terra4tandem_peak_names.txt
# rm $annotate/rloop4tandem_peak_names.txt



