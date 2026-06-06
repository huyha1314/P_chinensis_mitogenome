#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Folder structure: 
# data/    -> input files
# results/ -> output files
mkdir -p results
echo "Step 1: Extracting homologous blocks to CSV..."
pixi run python Fig4.py

echo "Step 2: Mapping genes to blocks..."
pixi run python mapgenetoblock.py

echo "Step 3: Identifying genes in inverted regions..."
pixi run python get_ivert_gene.py > results/inverted_genes_report.txt

echo "Step 4: Finding repeats at breakpoints (Smoking Gun)..."
pixi run python find_smoking.py > results/smoking_gun_report.txt

echo "Step 5: Generating Figure 4 Visualization (R/genoPlotR)..."
pixi run Rscript Fig4.R

echo "Step 6: Converting to publication-ready CMYK TIFF..."
pixi run python ../convert_cmyk.py results/Figure_4_Synteny.png results/Figure_4_Synteny_CMYK.tiff

echo "----------------------------------------------------"
echo "All steps complete! Check the 'results/' folder."
