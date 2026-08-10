library(dada2)

# 1. Paths to data and reference files
output_base <- "/mnt/d/IGIB_Data/Metagenomics_analysis/analysis_output"
db_base <- "/mnt/d/IGIB_Data/Metagenomics_analysis/SILVA_database"

seqtab.nochim <- readRDS(file.path(output_base, "seqtab_nochim.rds"))
train_set <- file.path(db_base, "silva_nr99_v138.2_toGenus_trainset.fa.gz")
species_set <- file.path(db_base, "silva_v138.2_assignSpecies.fa.gz")

# 2. Assign Taxonomy down to Genus level (Kingdom, Phylum, Class, Order, Family, Genus)
cat("\nStep 1: Assigning higher-level taxonomy (this utilizes significant RAM/CPU)...\n")
taxa <- assignTaxonomy(seqtab.nochim, train_set, multithread=TRUE)

# 3. Add Species-level annotations (Examines 100% exact string matching)
cat("\nStep 2: Refining exact species assignments...\n")
taxa <- addSpecies(taxa, species_set)

# 4. Save the structural binary file for downstream R packages (like Phyloseq)
saveRDS(taxa, file.path(output_base, "taxa_final.rds"))

# 5. Export a readable CSV layout matching your previously generated ASV IDs
asv_sequences <- colnames(seqtab.nochim)
asv_headers <- paste0("ASV_", seq_len(length(asv_sequences)))

taxa_exported <- taxa
rownames(taxa_exported) <- asv_headers
write.csv(taxa_exported, file.path(output_base, "asv_taxonomy.csv"), quote=FALSE)

cat("\nTaxonomic classification complete! Output saved to 'asv_taxonomy.csv'\n")
print(head(taxa_exported))
