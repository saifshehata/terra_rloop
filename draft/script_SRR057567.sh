#! /bin/bash
#$ -N SRR057567
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 8

project_dir=/cluster/khiom/sshehata001/proj/terra_rloop

# get reads
cd $project_dir/raw_data/reads
fastq-dump --gzip --split-3 SRR057567

# trim reads
cd $project_dir/raw_data/trimmed_reads
reads_dir=$project_dir/raw_data/reads
trimmomatic SE $reads_dir/SRR057567.fastq.gz SRR057567_trim.fq.gz ILLUMINACLIP:$ADAPTERS_SE:2:30:10 SLIDINGWINDOW:4:30 MINLEN:50

# map reads
cd $project_dir/results/mapped_reads
reads_dir=$project_dir/raw_data/trimmed_reads
REF=$project_dir/raw_data/genome/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna
bwa mem -t 8 $REF $reads_dir/SRR057567_trim.fq.gz | samtools sort -@ 5 -O bam -o SRR057567_se_sort.bam

# remove duplicates
samtools markdup -@ 5 -r -s -f SRR057567_se_sort_rmdup_stats SRR057567_se_sort.bam SRR057567_se_sort_rmdup.bam

# index
samtools index SRR057567_se_sort_rmdup.bam SRR057567_se_sort_rmdup.bai

# convert bam to bigwig
bamCoverage --bam SRR057567_se_sort_rmdup.bam --outFileName SRR057567_se_sort_rmdup.bw --outFileFormat bigwig --binSize 1 --numberOfProcessors 6 --normalizeUsing RPKM --effectiveGenomeSize 2652783500 --extendReads 130 --centerReads --ignoreDuplicates

