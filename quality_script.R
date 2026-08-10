# Load library
library(dada2)

# Defined path to the samples folder
path <- "/mnt/d/IGIB_Data/Metagenomics_analysis/raw_fastqs_26_samples"

#Describing the forward and reverse reads separately as a sorted list
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Picking only the basename of the samples
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Plotting the quality of the samples(as of only 2 out of the 26)
pdf("/mnt/d/IGIB_Data/Metagenomics_analysis/analysis_output/quality_plots.pdf", width=10, height=7)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

dev.off()
