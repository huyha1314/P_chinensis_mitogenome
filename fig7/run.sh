#!/usr/bin/env bash
# ==============================================================
# Figure 7: Intergenomic DNA Transfers (Circos)
# ==============================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p results

echo "Generating Fig7A..."
bash Fig7A.sh

echo "Generating Fig7B..."
bash Fig7B.sh

echo "Combining figures..."
pixi run python combine_fig7.py

echo "Converting to publication-ready CMYK TIFF..."
pixi run python ../convert_cmyk.py results/Fig7_Combined.png results/Fig7_Combined_CMYK.tiff

echo "Done!"