library(dada2)

# 1. Paths to your TRIMMED files (the output from your previous script)
filt_path <- "/mnt/d/IGIB_Data/Metagenomics_analysis/analysis_output/filtered_fastq_output"
output_base <- "/mnt/d/IGIB_Data/Metagenomics_analysis/analysis_output"

# 2. Grab all trimmed files safely
filtFs <- sort(list.files(filt_path, pattern="_R1_001.fastq.gz_trimmed.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(filt_path, pattern="_R2_001.fastq.gz_trimmed.fastq.gz", full.names = TRUE))

# Extract clean sample names
sample.names <- sapply(strsplit(basename(filtFs), "_R1"), `[`, 1)
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# 3. Learn Error Rates (The Machine Learning Step)
cat("\nStep 1: Learning error rates (this may take some time)...\n")
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Optional: Save the error plots to check later
pdf(file.path(output_base, "error_plots.pdf"))
plotErrors(errF, nominalQ=TRUE)
dev.off()

# 4. Sample Inference (Denoising)
cat("\nStep 2: Denoising sequences...\n")
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# 5. Merge Paired End Reads (Using your calculated 30bp overlap buffer)
cat("\nStep 3: Merging forward and reverse reads...\n")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# 6. Construct Sequence Table (ASV Count Matrix)
cat("\nStep 4: Building the Amplicon Sequence Variant (ASV) table...\n")
seqtab <- makeSequenceTable(mergers)

# 7. Remove Chimeras (PCR artifacts)
cat("\nStep 5: Removing chimeric sequences...\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# 8. Save the final processed matrix
saveRDS(seqtab.nochim, file.path(output_base, "seqtab_nochim.rds"))
cat("\nCore pipeline finished successfully! Saved final ASV table to 'seqtab_nochim.rds'\n")

