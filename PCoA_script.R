library(tidyverse)
library(vegan)

# ==============================================================================
# 1. LOAD YOUR RAW DATA
# ==============================================================================
abundance <- read_csv("/mnt/d/IGIB_Data/Metagenomics_analysis/analysis_output/comparative_study/asv_table_03_updated.csv")
metadata  <- read_csv("/mnt/d/IGIB_Data/Metagenomics_analysis/analysis_output/comparative_study/sample_metadata.csv") %>%
  mutate(Site = as.factor(Site))

# ==============================================================================
# 2. PREPARE STRICTOR NUMERIC MATRIX
# ==============================================================================
asv_df <- as.data.frame(abundance)
rownames(asv_df) <- asv_df$Sample_ID
asv_matrix <- as.matrix(asv_df[, -1])

# ==============================================================================
# 3. CALCULATE BETA DIVERSITY DISTANCE MATRIX
# ==============================================================================
# Calculate Bray-Curtis dissimilarity (standard for microbiome community structures)
# For presence/absence only, you could swap to method = "jaccard", binary = TRUE
bray_dist <- vegdist(asv_matrix, method = "bray")

# ==============================================================================
# 4. RUN PRINCIPAL COORDINATE ANALYSIS (PCoA)
# ==============================================================================
# cmdscale executes Classical Multidimensional Scaling (PCoA)
pcoa_output <- cmdscale(bray_dist, k = 2, eig = TRUE)

# Extract coordinates for Axis 1 and Axis 2
pcoa_data <- as.data.frame(pcoa_output$points) %>%
  rename(PCoA1 = V1, PCoA2 = V2) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(metadata, by = "Sample_ID")

# Calculate percent variation explained by each coordinate axis safely
eigenvalues <- pcoa_output$eig
# Filter out negative eigenvalues if any exist to maintain exact denominator bounds
eigenvalues_pos <- eigenvalues[eigenvalues > 0]
var_explained <- (eigenvalues_pos / sum(eigenvalues_pos)) * 100

# ==============================================================================
# 5. GENERATE AND SAVE THE CORRECT PCoA PLOT WITH LABELS
# ==============================================================================
library(ggrepel) # Load the package to manage text spacing safely

pcoa_plot <- ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2, color = Site)) +
  geom_point(size = 4, alpha = 0.8) +
  
  # --- ADDED: Auto-adjusting sample labels ---
  # 'label = Sample_ID' assigns text names from your data frame column
  geom_text_repel(
    aes(label = Sample_ID),
    size = 3,                   # Size of the font
    color = "black",            # Keeps text legible regardless of point color
    max.overlaps = Inf,         # Ensures all 24 labels are displayed
    box.padding = 0.3,          # Space around the text box
    point.padding = 0.3,        # Distance between the label and data point
    segment.color = 'grey50',   # Color of the small connector line
    segment.alpha = 0.6         # Transparency of the connector line
  ) +
  
  stat_ellipse(type = "t", linetype = 2, alpha = 0.4) + 
  theme_bw() + 
  labs(
    x = paste0("PCoA 1 (", round(var_explained[1], 1), "%)"),  
    y = paste0("PCoA 2 (", round(var_explained[2], 1), "%)"),  
    title = "PCoA of Microbial Beta Diversity (Bray-Curtis)",
    color = "Geographic Site"
  ) +
  theme(panel.grid.minor = element_blank())

print(pcoa_plot)
ggsave("/mnt/d/IGIB_Data/Metagenomics_analysis/analysis_output/comparative_study/pcoa_beta_diversity_labeled.png", 
       plot = pcoa_plot, width = 8, height = 6, dpi = 300) # Slightly widened to allow room for text

# ==============================================================================
# 6. STATISTICAL BETA DIVERSITY SIGNIFICANCE (PERMANOVA)
# ==============================================================================
cat("\n====================================================================\n")
cat(" PERMANOVA TEST STATISTICAL OUTPUT\n")
cat("====================================================================\n")

# Align metadata explicitly with distance vector matrices
metadata_aligned <- metadata %>%
  filter(Sample_ID %in% rownames(asv_matrix)) %>%
  arrange(match(Sample_ID, rownames(asv_matrix)))

# Run the adonis2 test to find out if "Site" clusters are truly different
permanova_res <- adonis2(bray_dist ~ Site, data = metadata_aligned, permutations = 999)
print(permanova_res)

