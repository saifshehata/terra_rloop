#! /bin/bash
#$ -N fetch_ref_genome
#$ -cwd
#$ -V
#$ -b n
#$ -j y
#$ -pe smp 2



# This script downloads the mouse mm10 reference genome

################################################################################################ 
# FETCH REFERENCE GENOME MM10 - THE ANALYSIS SET THAT IS MADE FOR SEQUENCE ANALYSIS PIPELINES
# for details, see https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
################################################################################################ 

# # load some modules
# module load bioinfo-tools
# module load bwa
# module load samtools

# # go directory
# project_dir=/proj/nb_storage/private/terra_rloop_project
# mkdir -p $project_dir/raw_data/ref_genome
# cd $project_dir/raw_data/ref_genome

# go directory
project_dir=/cluster/khiom/sshehata001/proj/terra_rloop
mkdir -p $project_dir/raw_data/genome
cd $project_dir/raw_data/genome

# define variable for the mm10 analysis set genome with no alt contigs
fasta=https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/seqs_for_alignment_pipelines/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.gz

# define variable for the ready-prepared index file for the reference genome (i.e. no need to use samtools faidx).
samidx=https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/seqs_for_alignment_pipelines/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.fai

# define variable for the ready-prepared bwa index file for the reference genome (i.e. no need to use bwa index).
bwaidx=https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/seqs_for_alignment_pipelines/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.bwa_index.tar.gz

# define variables for the GFF3 and GTF annotation files for that genome. NCBI Homo sapiens Updated Annotation Release 109.20190125 from 25 January 2019
gff=https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/seqs_for_alignment_pipelines/GCA_000001635.5_GRCm38.p3_full_analysis_set.refseq_annotation.gff.gz

# define variable for the chromosome sizes file.
chromSizes=https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

# # download the files
# echo -e "$fasta\n""$samidx\n""$bwaidx\n""$gff\n""$chromSizes" | parallel wget --timestamping {}

# download the files, but first check if the reference genome exists in unzipped format before re-downloading
if [ -f GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna ]; then
echo 'Unzipped reference genome already exists. Downloading other files...'
echo -e "$samidx\n""$bwaidx\n""$gff\n""$chromSizes" | parallel wget --timestamping {}
else 
echo -e "$fasta\n""$samidx\n""$bwaidx\n""$gff\n""$chromSizes" | parallel wget --timestamping {}
fi

# unzip the reference genome (if it is not already unzipped) for bwa mapping
if [ ! -f GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna ] && [ -f GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.gz ]; then
gunzip GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.gz
echo 'Reference genome unzipped'
fi

# unpack the bwa index .tar file if the files to be unpacked do not exist
if [ -f GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.sa ] && [ -f GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.amb ] && [ -f GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.ann ] && [ -f GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.pac ] && [ -f GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.bwt ]
then echo 'ALL tar files present. tar xvzf <file.bwa_index.tar.gz> will NOT be performed!'
else echo 'At least one tar file is missing. tar xvzf <file.bwa_index.tar.gz> is being performed...'
tar xvzf GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.bwa_index.tar.gz
fi


# alternative simpler commands without the if statements. The if statements are designed to avoid command repetition if the files are already there, e.g. the script is run several times

# # unpack the bwa index .tar file
# tar xvzf GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.bwa_index.tar.gz









# # download the mouse reference genome and its .gtf features file from ucsc. downloaded 2021-04-06.
# wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
# wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz # there are other versions as well on the website!
# wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

# # unzip downloaded genome files.
# gunzip mm10.fa.gz
# gunzip mm10.ncbiRefSeq.gtf.gz

# # index reference genome.
# bwa index mm10.fa # for bwa-mem aligner.
# samtools faidx mm10.fa # for IGV visualization.


