library(phyloseq)

# 1. Define paths
output_base <- "/mnt/d/IGIB_Data/Metagenomics_analysis/analysis_output"

# 2. Load DADA2 files
seqtab <- readRDS(file.path(output_base, "seqtab_nochim.rds"))
taxa <- readRDS(file.path(output_base, "taxa_final.rds"))

# 3. Build Phyloseq object and melt
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), tax_table(taxa))
df <- psmelt(ps)

# 4. Clean up NA values safely
df[is.na(df)] <- "Unclassified"

# ==============================================================================
# OPTION A: ONE PLOT WITH SAMPLES NESTED INSIDE
# Use this if you want the sample names to be the very first ring of the wheel.
# ==============================================================================
krona_nested <- data.frame(
  Abundance = df$Abundance,
  Sample    = df$Sample,  # Sample name becomes the root level ring
  Kingdom   = df$Kingdom,
  Phylum    = df$Phylum,
  Class     = df$Class,
  Order     = df$Order,
  Family    = df$Family,
  Genus     = df$Genus,
  Species   = df$Species
)

write.table(krona_nested, file.path(output_base, "krona_nested_samples.txt"), 
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
