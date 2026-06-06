#!/usr/bin/env bash
# ==============================================================
# Figure 5: Tanglegram (Cophylogeny) Visualization
# Nuclear vs Chloroplast & Mitochondrial vs Chloroplast
# ==============================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p results

echo "Running Figure 5 R script (tanglegram)..."
pixi run Rscript fig5.R
echo "Converting to publication-ready CMYK TIFF..."
pixi run python ../convert_cmyk.py results/Figure_5_Combined.png results/Figure_5_Combined_CMYK.tiff
echo "Done! Outputs saved to results/"
