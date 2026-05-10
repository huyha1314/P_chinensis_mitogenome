# 1. Load required libraries
library(phytools)
library(ape)

# 2. Load the three tree files
tree_nc   <- read.tree("data/nc.treefile")
tree_cp   <- read.tree("data/tree_cp_3.treefile")
tree_mito <- read.tree("data/pro_mito.treefile")

# 3. BULLETPROOF LABEL CLEANER
clean_labels <- function(labels) {
  labels <- gsub("_", "", labels)
  labels <- gsub(" ", "", labels)

  # Strip Accessions
  labels <- gsub("MZ[0-9]+\\.?[0-9]*", "", labels)
  labels <- gsub("NC[0-9]+\\.?[0-9]*", "", labels)
  
  # Re-insert formatting underscores
  labels <- gsub("Shorea", "Shorea_", labels)
  labels <- gsub("Dipterocarpus", "Dipterocarpus_", labels)
  labels <- gsub("Vatica", "Vatica_", labels)
  labels <- gsub("Hopea", "Hopea_", labels)
  labels <- gsub("OutgroupTheobromacacao", "Outgroup_Theobroma_cacao", labels)
  
  # Handle Parashorea group explicitly
  labels <- gsub("CH1Parashoreachinensis", "CH1_Parashorea_chinensis", labels)
  labels <- gsub("CH2Parashoreachinensis", "CH2_Parashorea_chinensis", labels)
  labels <- gsub("VN1Parashoreachinensis", "VN1_Parashorea_chinensis", labels)
  
  # Catch typos
  labels <- gsub("CH1Parashorechinensis", "CH1_Parashorea_chinensis", labels)
  labels <- gsub("CH2Parashorechinensis", "CH2_Parashorea_chinensis", labels)
  labels <- gsub("VN1Parashorechinensis", "VN1_Parashorea_chinensis", labels)
  
  return(labels)
}

# Apply the cleaner
tree_nc$tip.label   <- clean_labels(tree_nc$tip.label)
tree_mito$tip.label <- clean_labels(tree_mito$tip.label)
tree_cp$tip.label   <- clean_labels(tree_cp$tip.label)

# 4. ROOT THE TREES
out_name <- "Outgroup_Theobroma_cacao"
tree_nc   <- root(tree_nc, outgroup = out_name, resolve.root = TRUE)
tree_mito <- root(tree_mito, outgroup = out_name, resolve.root = TRUE)
tree_cp   <- root(tree_cp, outgroup = out_name, resolve.root = TRUE)

# 5. CREATE COPHYLOGENY OBJECTS
print("Computing Nuclear vs Chloroplast tanglegram...")
cophy_nc_cp   <- cophylo(tree_nc, tree_cp, rotate = FALSE)
print("Computing Mitochondrial vs Chloroplast tanglegram...")
cophy_mito_cp <- cophylo(tree_mito, tree_cp, rotate = FALSE)

# 6. BEAUTIFY LABELS FOR THE FINAL PLOT
clean_spaces <- function(cophy_obj) {
  cophy_obj$trees[[1]]$tip.label <- gsub("_", " ", cophy_obj$trees[[1]]$tip.label)
  cophy_obj$trees[[2]]$tip.label <- gsub("_", " ", cophy_obj$trees[[2]]$tip.label)
  cophy_obj$assoc[,1] <- gsub("_", " ", cophy_obj$assoc[,1])
  cophy_obj$assoc[,2] <- gsub("_", " ", cophy_obj$assoc[,2])
  return(cophy_obj)
}
cophy_nc_cp   <- clean_spaces(cophy_nc_cp)
cophy_mito_cp <- clean_spaces(cophy_mito_cp)

# 7. DEFINE HIGHLIGHT COLORS (Green shade matching iTOL)
cols_nc_cp   <- ifelse(grepl("Parashorea", cophy_nc_cp$assoc[,1]), "magenta", "mediumseagreen")
cols_mito_cp <- ifelse(grepl("Parashorea", cophy_mito_cp$assoc[,1]), "magenta", "mediumseagreen")

# Helper function for bootstraps
clean_bootstraps <- function(labels) {
  labels[is.na(labels)] <- ""
  return(labels)
}

# ==========================================
# --- PLOT A: NUCLEUS VS CHLOROPLAST ---
# ==========================================
pdf("results/Panel_A_Nuclear_vs_CP.pdf", width = 14, height = 8)
par(mar = c(1, 1, 3, 1)) # Add top margin for titles

plot(cophy_nc_cp, 
     link.type = "curved",    # Curved lines are better for crossed branches
     link.col = cols_nc_cp,   
     link.lwd = 2.5,          
     fsize = 0.9, 
     pts = FALSE)

title(main = "A. Nuclear vs. Chloroplast Phylogeny", cex.main = 1.8, font.main = 2, line = 1.5)
mtext("Nuclear Tree", side = 3, line = 0, adj = 0.15, cex = 1.2, font = 2)
mtext("Chloroplast Tree", side = 3, line = 0, adj = 0.85, cex = 1.2, font = 2)

# Labels are called AFTER plot() and set to black
nodelabels.cophylo(text = clean_bootstraps(cophy_nc_cp$trees[[1]]$node.label), 
                   which = "left", frame = "none", cex = 0.65, adj = c(1.2, -0.4), col = "black")
nodelabels.cophylo(text = clean_bootstraps(cophy_nc_cp$trees[[2]]$node.label), 
                   which = "right", frame = "none", cex = 0.65, adj = c(-0.2, -0.4), col = "black")
dev.off()


# ==========================================
# --- PLOT B: MITO VS CHLOROPLAST ---
# ==========================================
pdf("results/Panel_B_Mito_vs_CP.pdf", width = 14, height = 8)
par(mar = c(1, 1, 3, 1)) 

plot(cophy_mito_cp, 
     link.type = "curved", 
     link.col = cols_mito_cp, 
     link.lwd = 2.5, 
     fsize = 0.9, 
     pts = FALSE)

title(main = "B. Mitochondrial vs. Chloroplast Phylogeny", cex.main = 1.8, font.main = 2, line = 1.5)
mtext("Mitochondrial Tree", side = 3, line = 0, adj = 0.15, cex = 1.2, font = 2)
mtext("Chloroplast Tree", side = 3, line = 0, adj = 0.85, cex = 1.2, font = 2)

# Labels are called AFTER plot() and set to black
nodelabels.cophylo(text = clean_bootstraps(cophy_mito_cp$trees[[1]]$node.label), 
                   which = "left", frame = "none", cex = 0.65, adj = c(1.2, -0.4), col = "black")
nodelabels.cophylo(text = clean_bootstraps(cophy_mito_cp$trees[[2]]$node.label), 
                   which = "right", frame = "none", cex = 0.65, adj = c(-0.2, -0.4), col = "black")
dev.off()

print("Done! Check your folder for Panel_A_Nuclear_vs_CP.pdf and Panel_B_Mito_vs_CP.pdf")

# ==========================================
# --- FIGURE 5: COMBINED TANGLEGRAMS ---
# ==========================================
graphics.off()

# Helper to draw one panel
draw_combined_panel <- function(cophy_obj, link_cols, panel_letter, label_left, label_right) {
  if (panel_letter == "A") {
    title_text <- "A. Nuclear vs. Chloroplast Phylogeny"
  } else {
    title_text <- "B. Mitochondrial vs. Chloroplast Phylogeny"
  }
  
  # Plot the tanglegram
  plot(cophy_obj,
       link.type = "curved",
       link.col  = link_cols,
       link.lwd  = 2.5,
       fsize     = 0.9,
       pts       = FALSE)
  
  # phytools::plot.cophylo resets par() internally, so we force the title
  # into the current viewport using xpd=NA and figure coordinates
  usr <- par("usr")
  par(xpd = NA)
  
  # Main panel title: line=-2 pushes DOWN into the top of the plot area
  title(main = title_text, cex.main = 1.8, font.main = 2, line = -1)
  
  # Sub-labels: left and right tree names (placed just above the plot edge)
  mtext(label_left,  side = 3, line = 0.3, adj = 0.08, cex = 1.1, font = 2)
  mtext(label_right, side = 3, line = 0.3, adj = 0.92, cex = 1.1, font = 2)
  
  par(xpd = FALSE)
  
  nodelabels.cophylo(text = clean_bootstraps(cophy_obj$trees[[1]]$node.label),
                     which = "left",  frame = "none", cex = 0.65, adj = c(1.2, -0.4), col = "black")
  nodelabels.cophylo(text = clean_bootstraps(cophy_obj$trees[[2]]$node.label),
                     which = "right", frame = "none", cex = 0.65, adj = c(-0.2, -0.4), col = "black")
}

# --- OUTPUT: PDF ---
pdf("results/Figure_5_Combined.pdf", width = 14, height = 24)
par(mfrow = c(2, 1), oma = c(2, 2, 4, 2), mar = c(6, 2, 3, 2))
draw_combined_panel(cophy_nc_cp,   cols_nc_cp,   "A", "Nuclear Tree",       "Chloroplast Tree")
draw_combined_panel(cophy_mito_cp, cols_mito_cp, "B", "Mitochondrial Tree", "Chloroplast Tree")
dev.off()

# --- OUTPUT: PNG (600 dpi, publication quality) ---
png("results/Figure_5_Combined.png", width = 14, height = 24, units = "in", res = 900)
par(mfrow = c(2, 1), oma = c(2, 2, 4, 2), mar = c(6, 2, 3, 2))
draw_combined_panel(cophy_nc_cp,   cols_nc_cp,   "A", "Nuclear Tree",       "Chloroplast Tree")
draw_combined_panel(cophy_mito_cp, cols_mito_cp, "B", "Mitochondrial Tree", "Chloroplast Tree")
dev.off()

print("Done! results/Figure_5_Combined.pdf and results/Figure_5_Combined.png")