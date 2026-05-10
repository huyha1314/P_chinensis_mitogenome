# 1. Load required libraries
library(phytools)
library(ape)
library(grid)

# 2. Load the three tree files
tree_nc   <- read.tree("data/nc.treefile")
tree_cp   <- read.tree("data/tree_cp_3.treefile")
tree_mito <- read.tree("data/pro_mito.treefile")

# 3. BULLETPROOF LABEL CLEANER
clean_labels <- function(labels) {
  labels <- gsub("_", "", labels)
  labels <- gsub(" ", "", labels)
  labels <- gsub("MZ[0-9]+\\.?[0-9]*", "", labels)
  labels <- gsub("NC[0-9]+\\.?[0-9]*", "", labels)
  labels <- gsub("Shorea", "Shorea_", labels)
  labels <- gsub("Dipterocarpus", "Dipterocarpus_", labels)
  labels <- gsub("Vatica", "Vatica_", labels)
  labels <- gsub("Hopea", "Hopea_", labels)
  labels <- gsub("OutgroupTheobromacacao", "Outgroup_Theobroma_cacao", labels)
  labels <- gsub("CH1Parashoreachinensis", "CH1_Parashorea_chinensis", labels)
  labels <- gsub("CH2Parashoreachinensis", "CH2_Parashorea_chinensis", labels)
  labels <- gsub("VN1Parashoreachinensis", "VN1_Parashorea_chinensis", labels)
  return(labels)
}

tree_nc$tip.label   <- clean_labels(tree_nc$tip.label)
tree_mito$tip.label <- clean_labels(tree_mito$tip.label)
tree_cp$tip.label   <- clean_labels(tree_cp$tip.label)

# 4. ROOT THE TREES
out_name <- "Outgroup_Theobroma_cacao"
tree_nc   <- root(tree_nc, outgroup = out_name, resolve.root = TRUE)
tree_mito <- root(tree_mito, outgroup = out_name, resolve.root = TRUE)
tree_cp   <- root(tree_cp, outgroup = out_name, resolve.root = TRUE)

# 5. CREATE COPHYLOGENY OBJECTS
cophy_nc_cp   <- cophylo(tree_nc, tree_cp, rotate = FALSE)
cophy_mito_cp <- cophylo(tree_mito, tree_cp, rotate = FALSE)

# 6. BEAUTIFY LABELS
clean_spaces <- function(cophy_obj) {
  cophy_obj$trees[[1]]$tip.label <- gsub("_", " ", cophy_obj$trees[[1]]$tip.label)
  cophy_obj$trees[[2]]$tip.label <- gsub("_", " ", cophy_obj$trees[[2]]$tip.label)
  cophy_obj$assoc[,1] <- gsub("_", " ", cophy_obj$assoc[,1])
  cophy_obj$assoc[,2] <- gsub("_", " ", cophy_obj$assoc[,2])
  return(cophy_obj)
}
cophy_nc_cp   <- clean_spaces(cophy_nc_cp)
cophy_mito_cp <- clean_spaces(cophy_mito_cp)

cols_nc_cp   <- ifelse(grepl("Parashorea", cophy_nc_cp$assoc[,1]), "magenta", "mediumseagreen")
cols_mito_cp <- ifelse(grepl("Parashorea", cophy_mito_cp$assoc[,1]), "magenta", "mediumseagreen")

clean_bootstraps <- function(labels) {
  labels[is.na(labels)] <- ""
  labels[labels == "Root"] <- "Root    "
  return(labels)
}

# 7. PLOTTING FUNCTION
draw_combined_fig5 <- function() {
  # Giữ nguyên margin tổng thể
  par(mfrow = c(2, 1), mar = c(2, 0.5, 8, 0.5)) 
  
  # --- PANEL A ---
  plot(cophy_nc_cp, link.type="curved", link.col=cols_nc_cp, link.lwd=2.5, 
       fsize=1.8, font=3, pts=FALSE, mar=c(1, 0.2, 6, 0.2))
  par(xpd = NA) 
  nodelabels.cophylo(text = clean_bootstraps(cophy_nc_cp$trees[[1]]$node.label), 
                     which="left", frame="none", cex=0.9, adj=c(1.6, -0.5))
  nodelabels.cophylo(text = clean_bootstraps(cophy_nc_cp$trees[[2]]$node.label), 
                     which="right", frame="none", cex=0.9, adj=c(-0.6, -0.5))
  par(xpd = FALSE)
  
  # --- PANEL B ---
  plot(cophy_mito_cp, link.type="curved", link.col=cols_mito_cp, link.lwd=2.5, 
       fsize=1.8, font=3, pts=FALSE, mar=c(1, 0.2, 6, 0.2))
  par(xpd = NA)
  nodelabels.cophylo(text = clean_bootstraps(cophy_mito_cp$trees[[1]]$node.label), 
                     which="left", frame="none", cex=0.9, adj=c(1.6, -0.5))
  nodelabels.cophylo(text = clean_bootstraps(cophy_mito_cp$trees[[2]]$node.label), 
                     which="right", frame="none", cex=0.9, adj=c(-0.6, -0.5))
  par(xpd = FALSE)

  # 8. ADD TITLES USING GRID 
  # Chỉnh x=0.12 và x=0.88 để kéo các tiêu đề "Tree" lùi vào trong, tránh bị cắt chữ
  
  # Panel A Titles
  grid.text("A. Nuclear vs. Chloroplast Phylogeny", x=0.5, y=0.98, gp=gpar(cex=2.5, fontface="bold"))
  grid.text("Nuclear Tree", x=0.12, y=0.94, gp=gpar(cex=2.5, fontface="bold"))
  grid.text("Chloroplast Tree", x=0.88, y=0.94, gp=gpar(cex=2.5, fontface="bold"))
  
  # Panel B Titles
  grid.text("B. Mitochondrial vs. Chloroplast Phylogeny", x=0.5, y=0.48, gp=gpar(cex=2.5, fontface="bold"))
  grid.text("Mitochondrial Tree", x=0.12, y=0.44, gp=gpar(cex=2.5, fontface="bold"))
  grid.text("Chloroplast Tree", x=0.88, y=0.44, gp=gpar(cex=2.5, fontface="bold"))
}

# 9. OUTPUT
dir.create("results", showWarnings = FALSE)

pdf("results/Figure_5_Combined.pdf", width=20, height=13)
draw_combined_fig5()
dev.off()

png("results/Figure_5_Combined.png", width=20, height=13, units="in", res=600)
draw_combined_fig5()
dev.off()

cat("Success: Figure 5 Combined saved to results/\n")