#! /bin/bash
#$ -N map_reads
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 8

# this script maps the trimmed fastq reads to the reference genome

# # load some modules
# module load bioinfo-tools
# module load bwa
# module load deepTools # for bamCoverage
# module load samtools

# # go directory
# project_dir=/proj/nb_storage/private/terra_rloop_project
# mkdir -p $project_dir/results/mapped_reads
# cd $project_dir/results/mapped_reads
# trimmed_reads=$project_dir/data/fastq/trimmed

# go to directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
mkdir -p $project_dir/results/mapped_reads
cd $project_dir/results/mapped_reads
trimmed_reads=$project_dir/raw_data/trimmed_reads

# define reference genome base/index variable.
# REF=$project_dir/refs/mouse/GRCm38/mm10.fa
REF=$project_dir/raw_data/genome/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna

################################################################################################ 
# MAP SINGLE-END READS TO REFERENCE MM10
################################################################################################ 

# map single-end reads | sort by query position.
# ls $trimmed_reads | grep trim | grep -v unpaired | egrep -v '(_1|_2)' | cut -d_ -f1 | egrep 'SRR[0-9]+$' | parallel "echo 'bwa mem -t 8 $REF {}_trim.fq.gz | samtools view -@ 5 -bSu -T $REF | samtools sort -@ 5 > {}.sorted.bam'"
ls $trimmed_reads | grep -v unpaired | egrep 'SRR[0-9]+_trim.fq.gz' | cut -d_ -f1 | parallel "bwa mem -t 8 $REF $trimmed_reads/{}_trim.fq.gz | samtools sort -@ 5 -O bam -o {}_se_sort.bam"

# remove duplicates from single-end reads
ls | egrep 'SRR[0-9]+_se_sort.bam$' | parallel "samtools markdup -@ 5 -r -s -f {.}_rmdup_stats {} {.}_rmdup.bam"

# index bam file.
ls | egrep 'SRR[0-9]+_se_sort_rmdup.bam$' | parallel "samtools index {} {.}.bai"

# create .bw files for IGV
ls | egrep 'SRR[0-9]+_se_sort_rmdup.bam$' | parallel "bamCoverage --bam {} --outFileName {.}.bw --outFileFormat bigwig --binSize 1 --numberOfProcessors 6 --normalizeUsing RPKM --effectiveGenomeSize 2652783500 --extendReads 130 --centerReads --ignoreDuplicates"

###--effectiveGenomeSize 2652783500 according to https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html.
###--extendReads 130 \ # automatically determined for paired-end data.

################################################################################################ 
# MAP PAIRED-END READS TO REFERENCE MM10
################################################################################################ 

# map paired-end reads | sort by query name for fixmate
ls $trimmed_reads | grep -v unpaired | egrep 'SRR[0-9]+_[1-2]_trim.fq.gz' | cut -d_ -f1 | uniq | parallel "bwa mem -t 8 $REF $trimmed_reads/{}_1_trim.fq.gz $trimmed_reads/{}_2_trim.fq.gz | samtools sort -@ 5 -O bam -n -o {}_pe_nsort.bam"

# Add ms and MC tags for markdup to use later.
ls | egrep 'SRR[0-9]+_pe_nsort.bam$' | parallel "samtools fixmate -@ 5 -m {} {.}_fixmate.bam"

# resort by position. samtools markdup needs position order
ls | egrep 'SRR[0-9]+_pe_nsort_fixmate.bam$' | cut -d_ -f1,2 | parallel "samtools sort -@ 5 -o {}_sort.bam {}_nsort_fixmate.bam"

# remove duplicates from paired-end reads
ls | egrep 'SRR[0-9]+_pe_sort.bam$' | parallel "samtools markdup -@ 5 -r -s -f {.}_rmdup_stats {} {.}_rmdup.bam"

# index bam file.
ls | egrep 'SRR[0-9]+_pe_sort_rmdup.bam$' | parallel "samtools index {} {.}.bai"

# create .bw files for IGV
ls | egrep 'SRR[0-9]+_pe_sort_rmdup.bam$' | parallel "bamCoverage --bam {} --outFileName {.}.bw --outFileFormat bigwig --binSize 1 --numberOfProcessors 6 --normalizeUsing RPKM --effectiveGenomeSize 2652783500 --extendReads --centerReads --ignoreDuplicates"

# remove unnecessary files.
ls | grep nsort | parallel rm {}

