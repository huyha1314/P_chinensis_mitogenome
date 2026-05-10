#!/usr/bin/env Rscript

# ---------------------------------------------------------
# 1. Load Required Libraries Silently
# ---------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork) # For combining the plots
  library(tidyverse)
  library(ggtree)
  library(treeio)
  library(ggrepel)
  library(ape)      
})

# ---------------------------------------------------------
# 2. Setup File Paths & Shared Colors
# ---------------------------------------------------------
# Paths for RNA Editing Bar Plot
rna_file <- "data/VN.merge.tsv"

# Paths for MEME Tree Plot
csv_file  <- "data/Table_S1_Positive_Selection_Summary.csv"
tree_file <- "data/hyphy_annotated_tree.nwk"

# Output Files
output_pdf <- "results/Figure_7.pdf"
output_png <- "results/Figure_7.png"   

# Master Color Palette (Synchronizes both plots)
custom_colors <- c(
  "ccmFC" = "#F8766D",  
  "ccmFN" = "#A3A500",  
  "nad4"  = "#00BF7D",  
  "rps3"  = "#00B0F6",  
  "rps4"  = "#E76BF3",  
  "sdh4"  = "#E68613",  
  "Other" = "#9ECAE1"   
)

# ---------------------------------------------------------
# 3. Build TOP PLOT (RNA Editing Bar Chart)
# ---------------------------------------------------------
cat("Processing RNA Editing data...\n")
df_rna <- read.table(rna_file, header = FALSE, sep = "\t", quote = "", 
                 comment.char = "", fill = TRUE, stringsAsFactors = FALSE)

# Extract Gene and Chromosome
df_rna$Gene <- sub(".*chr[0-9]+_(.*)_CDS.*", "\\1", df_rna$V1)
df_rna <- df_rna[!is.na(df_rna$Gene) & df_rna$Gene != "" & df_rna$Gene != df_rna$V1, ]
df_rna$Chromosome <- ifelse(grepl("chr1", df_rna$V1), "Chromosome 1", 
                     ifelse(grepl("chr2", df_rna$V1), "Chromosome 2", "Unknown"))

# Count and assign colors
gene_counts <- aggregate(V1 ~ Gene + Chromosome, data = df_rna, FUN = length)
colnames(gene_counts) <- c("Gene", "Chromosome", "Count")
gene_counts$ColorGroup <- ifelse(gene_counts$Gene %in% names(custom_colors), gene_counts$Gene, "Other")

# Generate RNA Plot
p_rna <- ggplot(gene_counts, aes(x = Gene, y = Count, fill = ColorGroup)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.2) +
  geom_text(aes(label = Count), vjust = -0.8, size = 3.5, color = "black") +
  scale_fill_manual(values = custom_colors) +
  labs(y = "Number of RNA editing sites", x = NULL) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  facet_grid(~Chromosome, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 14, margin = margin(r = 10)), 
    axis.text.x = element_text(angle = 60, hjust = 1, size = 11, color = "black"), 
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "none", 
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 12, face = "bold", color = "black")
  )

# ---------------------------------------------------------
# 4. Build BOTTOM PLOT (Phylogenetic Tree)
# ---------------------------------------------------------
cat("Loading tree and MEME data...\n")
gene_map <- c(
  "OG0000009_cleaned" = "ccmFC", "OG0000010_cleaned" = "ccmFN",
  "OG0000018_cleaned" = "nad4",  "OG0000029_cleaned" = "rps3",
  "OG0000030_cleaned" = "rps4",  "OG0000031_cleaned" = "sdh4"
)

meme_data <- read_csv(csv_file, show_col_types = FALSE)
tree <- read.tree(tree_file)
tree <- ape::drop.tip(tree, "cacao_mt") 

selection_clean <- meme_data %>%
  mutate(Gene_Name = gene_map[Orthogroup]) %>%
  mutate(Node_Label = strsplit(as.character(`Affected Lineages (Tips/Nodes)`), ",\\s*")) %>%
  unnest(Node_Label) %>%
  mutate(Node_Label = str_replace_all(Node_Label, " ", "_")) %>%
  mutate(Mutation_Label = paste0(Gene_Name, " (Site ", `Site (Codon)`, ")")) %>%
  select(Node_Label, everything())

# Generate Tree Plot
p_tree <- ggtree(tree) %<+% selection_clean +  
  geom_tree(linewidth = 0.8) +
  geom_tiplab(aes(label = str_replace_all(label, "_", " ")), 
              size = 3.5, offset = 0.0001, fontface = "italic") +
  geom_point(aes(color = Gene_Name, size = `p-value`), alpha = 0.85, na.rm = TRUE) +
  geom_label_repel(aes(label = Mutation_Label, fill = Gene_Name), 
                   color = "white", size = 2.5, fontface = "bold",
                   box.padding = 0.8,      
                   point.padding = 0.5, 
                   nudge_x = -0.0005,      
                   force = 5,              
                   direction = "y", 
                   segment.color = 'grey50', segment.size = 0.5,
                   na.rm = TRUE, show.legend = FALSE) +
  theme_tree2() + 
  
  # Ensure custom colors map correctly
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  
  # FIX: Range added here to prevent the giant brown bubble
  scale_size_continuous(trans = "reverse", range = c(2, 8)) + 
  
  labs(subtitle = "Phylogenetic Distribution of Episodic Positive Selection (HyPhy MEME)",
       color = "Target Gene", size = "p-value") +
  theme(
    plot.subtitle = element_text(face = "italic", size = 12, hjust = 0.5),
    # FIX: Legend and axis theme items moved safely inside the theme() block
    legend.position = "right",
    legend.box.margin = margin(l = -20), 
    legend.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10, color = "black")
  ) +
  # FIX: xlim updated to 2.0 to give species names room to dodge the labels
  xlim(-0.002, max(ape::node.depth.edgelength(tree)) * 2.0)

# ---------------------------------------------------------
# 5. Combine and Export Graphic
# ---------------------------------------------------------
cat("Stitching plots together...\n")

master_plot <- p_rna / p_tree +
  plot_layout(heights = c(1, 1.5)) + 
  plot_annotation(
    title = "Mitochondrial RNA Editing and Positive Selection Landscape",
    tag_levels = 'A',
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  ) & 
  theme(plot.tag = element_text(size = 18, face = "bold"))
  
cat("Saving high-resolution plot to", output_pdf, "...\n")
ggsave(output_pdf, plot = master_plot, width = 12, height = 14, dpi = 600, bg="white")
ggsave(output_png, plot = master_plot, width = 12, height = 14, dpi = 600, bg="white")

cat("Process complete! Check the output files.\n")