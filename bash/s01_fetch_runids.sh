#! /bin/bash
#$ -N fetch_runids
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 1

# This script fetches the SRR run ID accessions for the fastq files to downloads

# # go to directory.
# project_dir=/proj/nb_storage/private/terra_rloop_project
# mkdir -p $project_dir/data/runinfo/
# cd $project_dir/data/runinfo/

# go directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
mkdir -p $project_dir/raw_data/runinfo
cd $project_dir/raw_data/runinfo

# fetch SRR IDs for mouse TERRA binding datasets in mouse from publication using GEO accession. SRA Bioproject nr: PRJNA315253
esearch -db gds -query GSE79180 | elink -target sra | efetch -format docsum | xtract -pattern DocumentSummary -element Run@acc,Bioproject,Title | grep ChIRT | egrep '(input|Sense|TERRA)' | awk '{print $1,$5,$6,$7,$8}' | sed 's/;//' > terra_runinfo.txt

# fetch SRR IDs for mouse R-loop binding datasets using GEO accession. SRA Bioproject nr: PRJNA287823
esearch -db gds -query GSE70189 | elink -target sra | efetch -format docsum | xtract -pattern DocumentSummary -element Run@acc,Bioproject,Title | grep "Mus musculus" | grep E14 | awk '{print $1,$4,$5}' | sed 's/;//' > rloop_runinfo.txt

# fetch SRR IDs for mouse ATRX binding datasets using GEO accession. SRA Bioproject nr: PRJNA127699. Works fine even with individual GSM IDs (e.g. GSM551138)
esearch -db gds -query GSE22162 | elink -target sra | efetch -format docsum | xtract -pattern DocumentSummary -element Run@acc,Bioproject,Title | grep ES | awk '{print $1,$4,$5}{print $2,$4,$5}' | grep -v PRJN | sed 's/GSM551138: //' | sed 's/ //2' > atrx_runinfo.txt

# extract SRR IDs from runinfo files into a separate text file.
cat terra_runinfo.txt | awk '{print $1}' > runids_to_fetch.txt
cat rloop_runinfo.txt | awk '{print $1}' >> runids_to_fetch.txt
cat atrx_runinfo.txt | awk '{print $1}' >> runids_to_fetch.txt

# create basenames
cat terra_runinfo.txt rloop_runinfo.txt atrx_runinfo.txt | tr ' ' '_' > basenames.txt