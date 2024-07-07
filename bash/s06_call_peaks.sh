#! /bin/bash
#$ -N call_peaks
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 3

# this script uses the mapped reads to call peaks using MACS2

# # load some modules
# module load bioinfo-tools
# module load MACS

# # go to directory.
# project_dir=/proj/nb_storage/private/terra_rloop_project
# mkdir -p $project_dir/results/peaks/
# cd $project_dir/results/peaks/
# basenames=$project_dir/data/runids/basenames.txt

# go to directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
mkdir -p $project_dir/results/peaks
cd $project_dir/results/peaks
mapped_reads=$project_dir/results/mapped_reads

# define variables for peak calling
TERRA_T=$mapped_reads/SRR2062968_pe_sort_rmdup.bam
TERRA_C1=$mapped_reads/SRR2062969_pe_sort_rmdup.bam
TERRA_C2=$mapped_reads/SRR2062971_pe_sort_rmdup.bam
ATRX_T=$mapped_reads/SRR057567_se_sort_rmdup.bam
ATRX_C=$mapped_reads/SRR080724_se_sort_rmdup.bam
RLOOP_T=$mapped_reads/SRR2075686_se_sort_rmdup.bam

# call peaks
macs2 callpeak -t $TERRA_T -c $TERRA_C1 -f BAMPE -g mm -B -n terra_sense
macs2 callpeak -t $TERRA_T -c $TERRA_C2 -f BAMPE -g mm -B -n terra_input
macs2 callpeak -t $ATRX_T -c $ATRX_C -f BAM -g mm -B -n atrx_input
macs2 callpeak -t $RLOOP_T -f BAM -g mm -B -n rloop




# # a seemingly more complex but more automatic way of selecting the files for peak calling

# # define variable
# basenames=$project_dir/raw_data/runinfo/basenames.txt

# # call peaks for TERRA vs Sense (paired-end) mapped reads. -g: effective genome size or organism (mm for mus musculus). -s: read length. -B/--bdg: store the fragment pileup, control lambda, -log10pvalue and -log10qvalue scores in bedGraph files. -n: The prefix string for output files. -q: q-value (minimum FDR) cutoff, default 0.05.
# cat $basenames | sort -n | egrep '(TERRA_ChIRT|Sense_ChIRT)' | cut -d_ -f1 | parallel --dryrun --max-args 2 "macs2 callpeak -t $mapped_reads/{1}_pe_sort_rmdup.bam -c $mapped_reads/{2}_pe_sort_rmdup.bam -f BAMPE -g mm -B -n {1}_terraVSsense"

# # call peaks for TERRA vs input (paired-end) mapped reads
# cat $basenames | sort -n | egrep '(TERRA_ChIRT|input_ChIRT)' | cut -d_ -f1 | parallel --dryrun --max-args 2 "macs2 callpeak -t $mapped_reads/{1}_pe_sort_rmdup.bam -c $mapped_reads/{2}_pe_sort_rmdup.bam -f BAMPE -g mm -B -n {1}_terraVSinput"

# # call peaks for ATRX (single-end) mapped reads
# cat $basenames | sort -n | egrep '(ATRX_Mouse_ES_ChIPseq|ATRX_Mouse_ES_Input)' | cut -d_ -f1 | parallel --dryrun --max-args 2 "macs2 callpeak -t $mapped_reads/{1}_se_sort_rmdup.bam -c $mapped_reads/{2}_se_sort_rmdup.bam -f BAM -g mm -B -n {1}_atrxVSinput"

# # call peaks for R-loop (single-end) mapped reads
# cat $basenames | sort -n | grep E14_DRIP-seq | cut -d_ -f1 | parallel --dryrun "macs2 callpeak -t $mapped_reads/{}_se_sort_rmdup.bam -f BAM -g mm -B -n {}_rloop"



# # testing
# TERRA_T=$mapped_reads/$(cat $basenames | sort -n | grep TERRA_ChIRT | cut -d_ -f1)_pe_sort_rmdup.bam
