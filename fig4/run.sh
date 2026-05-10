#!/bin/bash

# Folder structure: 
# data/    -> input files
# results/ -> output files
mkdir -p results
echo "Step 1: Converting Mauve backbone to CSV..."
micromamba run -n visualiz python Fig4.py

echo "Step 2: Mapping genes to blocks..."
micromamba run -n visualiz python mapgenetoblock.py

echo "Step 3: Finding genes in inverted regions..."
micromamba run -n visualiz python get_ivert_gene.py > results/inverted_genes_report.txt

echo "Step 4: Searching for repeats at breakpoints (Smoking Gun)..."
micromamba run -n visualiz python find_smoking.py > results/smoking_gun_report.txt

echo "Step 5: Generating Figure 4 Visualization (PNG & PDF)..."
micromamba run -n visualiz python plot_fig4.py

echo "----------------------------------------------------"
echo "All steps complete! Check the 'results/' folder."
