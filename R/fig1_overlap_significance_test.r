##########################################################################################################
# Statistical test to check if the overlap between TERRA and R-loop peaks (~700 peaks) would be expected to occur by chance.
##########################################################################################################

# set working directory.
setwd("/cluster/khiom/sshehata001/proj/terra_rloop/results/figures/Rplots")

# load rtracklayer to use import() in order to import bed files directly as GRanges objects.
library('rtracklayer')

# import TERRA and R-loop peaks as GRanges objects.
terra_all <- import("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/peak_lengths/terra_peaks_peak_lengths.bed", format = 'bed')
rloop_all <- import("/cluster/khiom/sshehata001/proj/terra_rloop/results/repeat_analysis/peak_lengths/rloop_peaks_peak_lengths.bed", format = 'bed')


#install regioneR package to use its overlapPermTest() for the statistical test.
#BiocManager::install('regioneR', update=FALSE)

# load regioneR package
library('regioneR')

#################################################################################################
# check installed genomes
installed.genomes()
# install mouse mm10 genome BioStrings file if not already installed, otherwise it will not work!
# BiocManager::install("BSgenome.Mmusculus.UCSC.mm10", update=FALSE)
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10.masked", update = FALSE)
#################################################################################################
library(BSgenome.Mmusculus.UCSC.mm10)

# run permutation test to check if overlap is by chance (for 100 permutations, p-value = 0.0099).Running with 1000 permutations may take ~30-60min.
pt <- overlapPermTest(A=terra_all, B=rloop_all, ntimes=1000, genome="mm10" , alternative = "auto")

# display the output of the test as text.
pt

# plot permutation test figure and save it as pdf.
pdf("fig1_overlap_significance_test.pdf", width = 4.5, height = 3.5)
plot(pt)
dev.off()


##########################################################################################################
# End 
##########################################################################################################
