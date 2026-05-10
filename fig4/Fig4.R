library(genoPlotR)
library(grid)

# 1. Load the backbone data
bbone_file <- "data/p_chinensis_rotated.backbone"
bbone <- read_mauve_backbone(bbone_file)

# 2. Assign unique colors to each LCB
n_blocks <- nrow(bbone$dna_segs[[1]])
colors <- rainbow(n_blocks)
for (i in 1:length(bbone$dna_segs)) {
  bbone$dna_segs[[i]]$col <- colors
  bbone$dna_segs[[i]]$fill <- colors
}

# 3. Define Labels and Plotting Function
dna_seg_labels <- expression(
  "VN1"~italic("Parashorea chinensis"),
  "CH1"~italic("Parashorea chinensis"),
  "CH2"~italic("Parashorea chinensis"),
  italic("Shorea roxburghii"),
  italic("Theobroma cacao")~"(Outgroup)"
)

draw_synteny_plot <- function() {
  # Use a viewport that is centered and leaves 10% space at top and bottom
  pushViewport(viewport(y=0.53, height=0.8, width=0.92))
  
  plot_gene_map(dna_segs=bbone$dna_segs, 
                comparisons=bbone$comparisons,
                main="",
                dna_seg_labels=dna_seg_labels,
                dna_seg_label_cex = 3.5,
                scale_cex = 15,
                annotation_height=1.2,
                scale=TRUE,
                plot_new=FALSE)
  popViewport()
  
  # Draw Title high up
  grid.text(expression(italic("Parashorea chinensis")~"Mitogenome Synteny (5-way comparison)"), 
            y=unit(0.97, "npc"), gp=gpar(cex=3, fontface="bold"))
  
  # Fix scale text: large and clearly below the line
  grid.edit("scale.text", gp = gpar(cex = 3), vjust = 2)
  
  # ---------------------------------------------------------
  # ADD LEGEND (Horizontal row, matching Mauve colors)
  # ---------------------------------------------------------
  
  # 1. Forward Alignment (Sử dụng "indianred" để khớp với màu đỏ nâu của plot)
  grid.polygon(x = c(0.25, 0.28, 0.28, 0.25), 
               y = c(0.055, 0.055, 0.065, 0.065), 
               gp = gpar(fill="indianred", alpha=0.6, col=NA))
  grid.text("Collinear alignment", 
            x = 0.29, y = 0.06, just="left", gp=gpar(cex=3, fontface="italic"))
  
  # 2. Inverted Alignment (Sử dụng "steelblue" để khớp với màu xanh xám của plot)
  grid.polygon(x = c(0.60, 0.63, 0.63, 0.60), 
               y = c(0.055, 0.055, 0.065, 0.065), 
               gp = gpar(fill="steelblue", alpha=0.6, col=NA))
  grid.text("Structural inversion", 
            x = 0.64, y = 0.06, just="left", gp=gpar(cex=3, fontface="italic"))
}

# 4. Generate Outputs
dir.create("results", showWarnings = FALSE)

# PDF
pdf("results/Figure_4_Synteny.pdf", width=17, height=12, onefile=FALSE)
draw_synteny_plot()
dev.off()

# PNG
png("results/Figure_4_Synteny.png", width=17, height=12, units="in", res=600)
draw_synteny_plot()
dev.off()

cat("Success: 5-way synteny plot saved to results/Figure_4_Synteny.pdf and results/Figure_4_Synteny.png\n")