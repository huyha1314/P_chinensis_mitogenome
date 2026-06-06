#!/usr/bin/env bash
# ==============================================================
# Figure 8: RNA Editing & Positive Selection Landscape
# ==============================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p results

echo "Running Figure 8 R script..."
pixi run Rscript Fig8.R
echo "Converting to publication-ready CMYK TIFF..."
pixi run python ../convert_cmyk.py results/Figure_8.png results/Figure_8_CMYK.tiff
echo "Done! Outputs saved to results/"
