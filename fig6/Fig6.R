# Load required libraries
suppressMessages(library(treeio))
suppressMessages(library(ggtree))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(ape))
suppressMessages(library(dplyr))

# ==========================================
# 1. AUTO-FIX HPD TO TEXT FORMAT
# ==========================================
raw_text <- readLines("data/FigTree.tre")
fixed_text <- gsub("\\[&95%HPD=\\[(.*?)\\]\\]", "[&HPD_95=\"[\\1]\"]", raw_text)
fixed_text <- gsub("\\[&95%HPD=\\{(.*?)\\}\\]", "[&HPD_95=\"[\\1]\"]", fixed_text)
writeLines(fixed_text, "data/FigTree_fixed.tre")
tree <- read.beast("data/FigTree_fixed.tre")
tree@phylo$root.edge <- NULL

# ==========================================
# 2. COMPLETE NOMENCLATURE CLEANUP
# ==========================================
labels <- tree@phylo$tip.label
labels <- sub("_cp", "", labels)

# Update Genus Names to Shorea
labels <- sub("Anthoshorea_henryana", "Shorea_henryana", labels)
labels <- sub("Rubroshorea_leprosula", "Shorea_leprosula", labels)
labels <- sub("Anthoshorea_roxburghii", "Shorea_roxburghii", labels)

# Format Accessions: COMPLETELY REMOVE ACCESSION NUMBERS
labels <- sub("_MZ[A-Za-z0-9.]+", "", labels)
labels <- sub("_NC_[A-Za-z0-9.]+", "", labels)
labels <- gsub("_", " ", labels)
labels[labels == "cacao"] <- "Theobroma cacao"

# Format Focal Species
focal_spp <- c("Parashorea chinensis (VN1)", "Parashorea chinensis (CH1)", "Parashorea chinensis (CH2)")
labels[labels == "VN1"] <- focal_spp[1]
labels[labels == "CH1"] <- focal_spp[2]
labels[labels == "CH2"] <- focal_spp[3]
tree@phylo$tip.label <- labels

# ==========================================
# PANEL A: FULL DEEP-TIME PHYLOGENY
# ==========================================
pA <- ggtree(tree, linewidth = 1, color = "black")

# Background Taxonomic Shading (Perfectly mapped to monophyletic blocks)
rect_layers <- list(
  annotate("rect", xmin = -1.3, xmax = 0.02, ymin = -0.5, ymax = 1.5, fill = "#E8F5E9", alpha = 1), # Malvaceae (soft green)
  annotate("rect", xmin = -1.3, xmax = 0.02, ymin = 1.5, ymax = 9.5, fill = "#FCE4EC", alpha = 1), # Tribe Shoreae (soft pink)
  annotate("rect", xmin = -1.3, xmax = 0.02, ymin = 9.5, ymax = 16.5, fill = "#F3E5F5", alpha = 1) # Tribe Dipterocarpeae (soft purple)
)
pA$layers <- c(rect_layers, pA$layers)

# Calculate Node Ages and 95% HPD in Ma
d <- pA$data
max_x <- max(d$x)

d <- d %>%
  mutate(
    age_ma = (max_x - x) * 100,
    hpd_lower = sapply(HPD_95, function(val) {
      if (is.null(val) || any(is.na(val)) || length(val) < 1) return(NA)
      clean <- gsub("\\[|\\]|\\s", "", val[1])
      as.numeric(clean) * 100
    }),
    hpd_upper = sapply(HPD_95, function(val) {
      if (is.null(val) || any(is.na(val)) || length(val) < 2) return(NA)
      clean <- gsub("\\[|\\]|\\s", "", val[2])
      as.numeric(clean) * 100
    }),
    node_label = ifelse(!is.na(hpd_lower) & !is.na(hpd_upper) & age_ma >= 1.0,
                        sprintf("%.1f\n(%.1f-%.1f)", age_ma, hpd_lower, hpd_upper),
                        ifelse(age_ma >= 1.0, sprintf("%.1f", age_ma), ""))
  )
pA$data <- d

pA <- pA +
  # FIX: Prevent text from being hidden behind branches/nodes by nudging left/up and using vjust=0 (bottom-aligned)
  geom_nodelab(data = function(df) {
                 # Nudge Node 28 left to prevent overlap with Node 29, keeping it hovering cleanly above its branch line
                 df$x[df$node == 28] <- df$x[df$node == 28] - 0.06
                 df
               },
               aes(label = node_label), color = "steelblue4", size = 3.86, 
               nudge_x = -0.015, nudge_y = 0.25, vjust = 0, hjust = 1, fontface = "bold") +
  
  # Calibration / Support Dots
  geom_nodepoint(color = "black", size = 3) +
  
  # Plot Standard Tip Labels in Black
  geom_tiplab(data = function(x) x[!(x$label %in% focal_spp), ], 
              size = 6.2, fontface = "italic", offset = 0.02, color = "black") +
  
  # Plot Focal Species in RED
  geom_tiplab(data = function(x) x[x$label %in% focal_spp, ], 
              color = "darkred", size = 6.2, fontface = "bold.italic", offset = 0.02) +
  
  # FIX: Annotate text directly onto the background colors
  annotate("text", x = -1.28, y = 0.5, label = "Malvaceae", fontface = "bold", size = 6.9, hjust = 0, color = "darkgreen", alpha = 0.7) +
  annotate("text", x = -1.28, y = 5.5, label = "Tribe Shoreae", fontface = "bold", size = 6.9, hjust = 0, color = "indianred", alpha = 0.7) +
  annotate("text", x = -1.28, y = 13.0, label = "Tribe Dipterocarpeae", fontface = "bold", size = 6.9, hjust = 0, color = "purple4", alpha = 0.7) +
  
  # FIX: Bounding box around the Parashorea clade (source material for Panel B)
  annotate("rect", xmin = -0.02, xmax = 0.018, ymin = 4.5, ymax = 7.5, 
           color = "darkred", fill = NA, linetype = "dashed", linewidth = 0.8) +
           
  # TAXONOMIC TRIBE BRACKETS (Names updated to match clean tip labels)
  geom_strip(taxa1 = "Dipterocarpus zeylanicus", taxa2 = "Vatica rassak", 
             label = "Tribe Dipterocarpeae", barsize = 1.5, offset = 0.25, offset.text = 0.03, 
             color = "black", textcolor = "black", fontsize = 5.5, fontface = "bold") +
             
  geom_strip(taxa1 = "Hopea odorata", taxa2 = "Shorea roxburghii", 
             label = "Tribe Shoreae", barsize = 1.5, offset = 0.25, offset.text = 0.03, 
             color = "black", textcolor = "black", fontsize = 5.5, fontface = "bold") +
             
  geom_strip(taxa1 = "Theobroma cacao", taxa2 = "Theobroma cacao", 
             label = "Malvaceae", barsize = 1.5, offset = 0.25, offset.text = 0.03, 
             color = "black", textcolor = "black", fontsize = 5.5, fontface = "bold") +
             
  theme_tree2() +
  labs(
    title = "A. Phylogenomics and Evolutionary History of Dipterocarpaceae",
    x = "Millions of Years Ago (Ma)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 22.1, hjust = 0, margin = margin(b = 10)),
    axis.title.x = element_text(size = 19.3, face = "bold", margin = margin(t = 25)), 
    axis.text.x = element_text(size = 16.6, color = "black"),
    axis.line.x = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.8, color = "black"),
    # Adjusted right margin back to 386 to prevent long tip labels from truncating
    plot.margin = margin(t=10, r=386, b=0, l=10) 
  )

# Strict cutoff at 0.02 using expand = c(0,0) and coord_cartesian to widen right boundary
pA <- revts(pA) + 
  scale_x_continuous(
    breaks = seq(-1.2, 0, by = 0.2),
    labels = function(x) abs(x) * 100,
    expand = c(0, 0)
  ) + 
  scale_y_continuous(limits = c(-1.5, 17.5), expand = c(0, 0)) + 
  coord_cartesian(xlim = c(-1.3, 0.02), ylim = c(-1.5, 17.5), expand = FALSE, clip = "off") 

# FIX: Close the gap with the X-axis by using ymin = -1.5 for the geological boxes to exactly match the plot limits
pA <- pA + 
  annotate("rect", xmin = -1.3, xmax = -0.66, ymin = -1.5, ymax = -0.5, fill = "#7FC64E", alpha = 0.8) +
  annotate("text", x = -0.98, y = -1.0, label = "Cretaceous", size = 5.5, fontface = "bold", color = "white") +
  
  annotate("rect", xmin = -0.66, xmax = -0.2303, ymin = -1.5, ymax = -0.5, fill = "#FDA75F", alpha = 0.8) +
  annotate("text", x = -0.445, y = -1.0, label = "Paleogene", size = 5.5, fontface = "bold", color = "white") +
  
  annotate("rect", xmin = -0.2303, xmax = -0.0258, ymin = -1.5, ymax = -0.5, fill = "#FFE619", alpha = 0.8) +
  annotate("text", x = -0.14, y = -1.0, label = "Neogene", size = 5.0, fontface = "bold", color = "black") +
  
  annotate("rect", xmin = -0.0258, xmax = 0.02, ymin = -1.5, ymax = -0.5, fill = "#F9F97F", alpha = 0.8) +
  annotate("text", x = -0.001, y = -1.0, label = "Quaternary", size = 4.0, fontface = "bold", color = "black")

# ==========================================
# PANEL B: ZOOMED PLEISTOCENE ISOLATION
# ==========================================
non_focal <- tree@phylo$tip.label[!tree@phylo$tip.label %in% focal_spp]
tree_B <- treeio::drop.tip(tree, non_focal)

pB <- ggtree(tree_B, linewidth = 1.2, color = "black")

# Calculate Node Ages and 95% HPD in kya for Panel B
d_B <- pB$data
max_x_B <- max(d_B$x)

d_B <- d_B %>%
  mutate(
    # Multiply raw branch length by 100,000 to convert to kya
    age_kya = (max_x_B - x) * 100000,
    
    hpd_lower = sapply(HPD_95, function(val) {
      if (is.null(val) || any(is.na(val)) || length(val) < 1) return(NA)
      clean <- gsub("\\[|\\]|\\s", "", val[1])
      as.numeric(clean) * 100000
    }),
    
    hpd_upper = sapply(HPD_95, function(val) {
      if (is.null(val) || any(is.na(val)) || length(val) < 2) return(NA)
      clean <- gsub("\\[|\\]|\\s", "", val[2])
      as.numeric(clean) * 100000
    }),
    
    # Format the labels as Age (Lower - Upper)
    node_label = ifelse(!is.na(hpd_lower) & !is.na(hpd_upper),
                        sprintf("%.1f\n(%.1f-%.1f)", age_kya, hpd_lower, hpd_upper),
                        sprintf("%.1f", age_kya))
  )
pB$data <- d_B

pB <- pB +
  # Add formatted divergence times and HPD intervals for the zoomed panel
  geom_nodelab(aes(label = node_label), color = "darkred", size = 4.8, 
               nudge_x = -0.00001, nudge_y = 0.1, vjust = 0, hjust = 1, fontface = "bold") +
  
  geom_nodepoint(color = "darkred", size = 3) +
  
  geom_tiplab(color = "darkred", size = 8.3, fontface = "bold.italic", offset = 0.00001) +
  
  theme_tree2() +
  labs(
    title = "B. Late Pleistocene Isolation of Focal Species (Hainan Population CH2)",
    x = "Thousands of Years Ago (kya)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 22.1, hjust = 0, margin = margin(b = 10)),
    # FIX: Emphasize the difference between the two X-axes by making Panel B's axis distinct (dark red)
    axis.title.x = element_text(size = 19.3, face = "bold", color = "darkred", margin = margin(t = 10)),
    axis.text.x = element_text(size = 16.6, face = "bold", color = "darkred"),
    axis.line.x = element_line(linewidth = 0.8, color = "darkred"),
    axis.ticks.x = element_line(linewidth = 0.8, color = "darkred"),
    plot.margin = margin(t=0, r=386, b=10, l=10) # Restored margin to match Panel A
  )

# Strict cutoff at 0 using expand = c(0,0)
pB <- revts(pB) + 
  scale_x_continuous(
    breaks = c(-0.001, -0.00075, -0.0005, -0.00025, 0),
    labels = function(x) abs(x) * 100000,
    expand = c(0, 0) # STOPS AXIS STRICTLY AT 0
  ) +
  coord_cartesian(xlim = c(-0.001, 0), clip = "off")

pB <- pB + 
  geom_vline(xintercept = -0.000336, color = "darkred", linetype = "dashed", linewidth = 1.2) +
  # FIX: Added uncertainty (95% HPD confidence interval) to the split
  annotate("text", x = -0.000350, y = 1.5, label = "Hainan Isolation", 
           color = "darkred", hjust = 1, fontface = "bold", size = 7.6)

# ==========================================
# COMBINE WITH PATCHWORK & EXPORT
# ==========================================
composite_plot <- pA / pB + plot_layout(heights = c(1.3, 1))

pdf_out <- "results/Figure_6.pdf"
ggsave(pdf_out, plot = composite_plot, width = 14, height = 16, device = "pdf")

png_out <- "results/Figure_6.png"
ggsave(png_out, plot = composite_plot, width = 14, height = 16, dpi = 300, bg = "white", device = png)

cat(paste0("Success! Saved to ", pdf_out, " AND ", png_out, "\n"))
