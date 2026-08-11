library(phyloseq)

# 1. Define paths
output_base <- "/mnt/d/IGIB_Data/Metagenomics_analysis/analysis_output"

# 2. Load DADA2 files
seqtab <- readRDS(file.path(output_base, "seqtab_nochim.rds"))
taxa <- readRDS(file.path(output_base, "taxa_final.rds"))

# 3. Build Phyloseq object and aggregate raw counts
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), tax_table(taxa))
df <- psmelt(ps)

# 4. Clean up NA values across levels safely
df[is.na(df)] <- "Unclassified"

# 5. Build Galaxy-style hierarchy: Count followed by taxonomic levels
krona_data <- data.frame(
  Abundance = df$Abundance,
  Kingdom   = df$Kingdom,
  Phylum    = df$Phylum,
  Class     = df$Class,
  Order     = df$Order,
  Family    = df$Family,
  Genus     = df$Genus,
  Species   = df$Species
)

# 6. Save text file
write.table(krona_data, file.path(output_base, "krona_input.txt"), 
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

cat("\nKrona text structure saved to 'krona_input.txt'\n")

