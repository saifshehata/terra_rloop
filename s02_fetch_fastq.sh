#! /bin/bash
#$ -N fetch_fastq
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 2


# This script downloads fastq files using SRR run ID accessions

# # load some modules
# module load bioinfo-tools
# module load sratools

# # go to directory
# project_dir=/proj/nb_storage/private/terra_rloop_project
# mkdir -p $project_dir/data/fastq/
# cd $project_dir/data/fastq/

# go to directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
mkdir -p $project_dir/raw_data/reads
cd $project_dir/raw_data/reads

# define variable for runids file with SRA accessions
runids=$project_dir/raw_data/runinfo/runids_to_fetch.txt

# download fastq files
cat $runids | parallel fastq-dump --gzip --split-3 {}

# compress files to save space.
# ls | grep .fastq | parallel gzip --verbose --best {}
    #alt# gzip -r



