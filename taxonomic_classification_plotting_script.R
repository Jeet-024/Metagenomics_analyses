library(phyloseq)
library(ggplot2)

# 1. Define paths
output_base <- "/mnt/d/IGIB_Data/Metagenomics_analysis/analysis_output"

# 2. Load the binary DADA2 output files
seqtab <- readRDS(file.path(output_base, "seqtab_nochim.rds"))
taxa <- readRDS(file.path(output_base, "taxa_final.rds"))

# 3. Create a unified Phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               tax_table(taxa))

# Clean up sample names for the plot axis labels
sample_names(ps) <- sapply(strsplit(sample_names(ps), "_R1"), `[`, 1)

# 4. Transform raw counts to Relative Abundance (Percentages)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# 5. Filter out rare taxa (Keep Phyla with more than 1 percent average abundance)
ps_top <- filter_taxa(ps_rel, function(x) mean(x) > 0.01, TRUE)

# 6. Generate the Phylum-level Stacked Bar Plot
p_phylum <- plot_bar(ps_top, fill="Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  theme_bw() +
  labs(title="Microbiome Composition at Phylum Level (Above 1 Percent Abundance)", 
       x="Samples", y="Relative Abundance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# 7. Save the plot to a high-resolution PDF
pdf(file.path(output_base, "taxonomic_composition_phylum.pdf"), width=12, height=7)
print(p_phylum)
dev.off()

cat("\nPlotting complete! Saved composition chart to 'taxonomic_composition_phylum.pdf'\n")

