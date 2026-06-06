#!/usr/bin/env bash
# ==============================================================
# Figure 6: Phylogenetic Tree (Deep-time & Late Pleistocene Zoom)
# ==============================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p results data

echo "Running Figure 6 R script..."
pixi run Rscript Fig6.R
echo "Converting to publication-ready CMYK TIFF..."
pixi run python ../convert_cmyk.py results/Figure_6.png results/Figure_6_CMYK.tiff
echo "Done! Outputs saved to results/"
