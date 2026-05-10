library(genoPlotR)

# 1. Load the backbone data
bbone_file <- "data/p_chinensis_rotated.backbone"
bbone <- read_mauve_backbone(bbone_file)

# 2. Assign unique colors to each LCB (Local Collinear Block)
# This loops dynamically, so it handles all 5 genomes automatically
n_blocks <- nrow(bbone$dna_segs[[1]])
colors <- rainbow(n_blocks)

for (i in 1:length(bbone$dna_segs)) {
  bbone$dna_segs[[i]]$col <- colors
  bbone$dna_segs[[i]]$fill <- colors
}

# 3. Generate the plot
# Created directory if not exists (handled by run.sh but good to be safe)
dir.create("results", showWarnings = FALSE)

# Generate PDF
pdf("results/Figure_4_Synteny.pdf", width=15, height=12)
plot_gene_map(dna_segs=bbone$dna_segs, 
              comparisons=bbone$comparisons,
              main="Parashorea chinensis Mitogenome Synteny (5-way comparison)",
              # Labels MUST match the exact order of the .fasta files in your bash command
              dna_seg_labels=c("VN1_parashorea chinensis'", 
                               "CH1_parashorea chinensis", 
                               "CH2_parashorea chinensis", 
                               "Shorea_roxburghii", 
                               "Cacao_Outgroup"),
              annotation_height=1.2,
              scale=TRUE)
dev.off()

# Generate PNG (High-res)
png("results/Figure_4_Synteny.png", width=15, height=12, units="in", res=300)
plot_gene_map(dna_segs=bbone$dna_segs, 
              comparisons=bbone$comparisons,
              main="Parashorea chinensis Mitogenome Synteny (5-way comparison)",
              dna_seg_labels=c("VN1_parashorea chinensis'", 
                               "CH1_parashorea chinensis", 
                               "CH2_parashorea chinensis", 
                               "Shorea_roxburghii", 
                               "Cacao_Outgroup"),
              annotation_height=1.2,
              scale=TRUE)
dev.off()

cat("Success: 5-way synteny plot saved to results/Figure_4_Synteny.pdf and results/Figure_4_Synteny.png\n")