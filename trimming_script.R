library(dada2)

# 1. Path to your raw files
path <- "/mnt/d/IGIB_Data/Metagenomics_analysis/raw_fastqs_26_samples"

# 2. Grab all fastq.gz files safely
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Check if any files were found
if(length(fnFs) == 0) stop("No Forward (_R1) files found in raw_fastqs/")
if(length(fnRs) == 0) stop("No Reverse (_R2) files found in raw_fastqs/")

# 3. Create the output directories using absolute paths (Safe!)
output_base <- "/mnt/d/IGIB_Data/Metagenomics_analysis/analysis_output"
if(!dir.exists(output_base)) dir.create(output_base, recursive = TRUE)

filt_path <- file.path(output_base, "filtered_fastq_output")
if(!dir.exists(filt_path)) dir.create(filt_path, recursive = TRUE)

# 4. Generate guaranteed distinct output names inside the correct folder
filtFs <- file.path(filt_path, paste0(basename(fnFs), "_trimmed.fastq.gz"))
filtRs <- file.path(filt_path, paste0(basename(fnRs), "_trimmed.fastq.gz"))

# 5. Run the trimming step
cat("Starting filtering and trimming...\n")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     trimLeft=c(17, 21),        # Trims primers
                     truncLen=c(275, 215),      # Truncates lengths
                     maxN=0, 
                     maxEE=c(2, 6),             # Drops poor reads
                     truncQ=2, 
                     rm.phix=TRUE,
                     compress=TRUE, 
                     multithread=TRUE)

# 6. Save the tracking table and print results
summary_path <- file.path(output_base, "trimming_summary.txt")
write.table(out, summary_path, sep="\t", quote=FALSE)

cat("\nTrimming complete! Summary of the first few files:\n")
print(head(out))


