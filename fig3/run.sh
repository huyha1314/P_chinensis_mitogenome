#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "Running Figure 3 plotting script..."
pixi run python plot_fig3.py
echo "Converting to publication-ready CMYK TIFF..."
pixi run python ../convert_cmyk.py results/Fig3_Circular_Genomic_Map.png results/Fig3_Circular_Genomic_Map_CMYK.tiff
echo "Done!"
