#! /bin/bash
#$ -N trim_reads
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 2


# this script trims the fastq reads

# # load some modules
# module load bioinfo-tools
# module load sratools
# module load gnuparallel
# module load trimmomatic
# module load java

# # go to fastq directory
# project_dir=/proj/nb_storage/private/terra_rloop_project
# mkdir -p $project_dir/data/fastq/trimmed/
# cd $project_dir/data/fastq/

# go to directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
mkdir -p $project_dir/raw_data/trimmed_reads
cd $project_dir/raw_data/trimmed_reads
reads_dir=$project_dir/raw_data/reads

# list the adapters in trimmomatic to decide which of them to use for trimming. I chose TruSeq3-SE.fa.
# ls /sw/apps/bioinfo/trimmomatic/0.36/rackham/adapters/ 
ls -l /cluster/khiom/sshehata001/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic/adapters 

# define variable for chosen adapter.
# ADAPTERS_SE=/sw/apps/bioinfo/trimmomatic/0.36/rackham/adapters/TruSeq3-SE.fa
# ADAPTERS_PE=/sw/apps/bioinfo/trimmomatic/0.36/rackham/adapters/TruSeq3-PE.fa
# define adapter to use..
ADAPTERS_SE=/cluster/khiom/sshehata001/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic/adapters/TruSeq3-SE.fa
    #alt# adapter=/cluster/khiom/sshehata001/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa
ADAPTERS_PE=/cluster/khiom/sshehata001/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic/adapters/TruSeq3-PE.fa

# trim single-end reads
# ls | egrep -v '(_1|_2)' | cut -d. -f1 | parallel 'trimmomatic SE {}.fastq.gz {.}_trim.fq.gz SLIDINGWINDOW:4:30 MINLEN:50'
ls $reads_dir | egrep -v '(_1|_2)' | cut -d. -f1 | parallel "trimmomatic SE $reads_dir/{}.fastq.gz {}_trim.fq.gz ILLUMINACLIP:$ADAPTERS_SE:2:30:10 SLIDINGWINDOW:4:30 MINLEN:50"

# trim paired-end reads
# ls $reads_dir | egrep '(_1|_2)' | cut -d. -f1 | parallel --max-args 2 'trimmomatic PE $reads_dir/{1}.fastq.gz $reads_dir/{2}.fastq.gz {1}_trim.fq.gz {1}_unpaired_trim.fq.gz {2}_trim.fq.gz {2}_unpaired_trim.fq.gz SLIDINGWINDOW:4:30 MINLEN:30'
ls $reads_dir | egrep '(_1|_2)' | cut -d. -f1 | parallel --max-args 2 "trimmomatic PE $reads_dir/{1}.fastq.gz $reads_dir/{2}.fastq.gz {1}_trim.fq.gz {1}_unpaired_trim.fq.gz {2}_trim.fq.gz {2}_unpaired_trim.fq.gz ILLUMINACLIP:$ADAPTERS_PE:2:30:10 SLIDINGWINDOW:4:30 MINLEN:30"


# move all trimmed files to the trimmed reads directory
# mv *trim.fq.gz $project_dir/raw_data/trimmed_reads