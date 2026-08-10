library(dada2)

# 1. Load the binary RDS file
output_base <- "/mnt/d/IGIB_Data/Metagenomics_analysis/analysis_output"
seqtab.nochim <- readRDS(file.path(output_base, "seqtab_nochim.rds"))

# 2. Extract the actual DNA sequences
asv_sequences <- colnames(seqtab.nochim)

# 3. Create clean, short IDs for your ASVs (e.g., ASV_1, ASV_2...)
asv_headers <- paste0("ASV_", seq_len(length(asv_sequences)))

# 4. Save the FASTA file (Contains the actual DNA sequences for alignment)
asv_fasta <- c(rbind(paste0(">", asv_headers), asv_sequences))
write(asv_fasta, file.path(output_base, "asvs.fasta"))

# 5. Save the ASV Count Matrix (Rows = Samples, Columns = ASV IDs)
# Rename columns from long DNA strings to short ASV IDs for readability
asv_tab <- seqtab.nochim
colnames(asv_tab) <- asv_headers
write.csv(asv_tab, file.path(output_base, "asv_table.csv"), quote=FALSE)

cat("\nExport complete! Created 'asvs.fasta' and 'asv_table.csv' in your output folder.\n")

