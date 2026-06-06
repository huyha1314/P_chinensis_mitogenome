#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p results
echo "Generating Geographic Map (Figure 1)..."
pixi run Rscript Fig1.R
echo "Converting to publication-ready CMYK TIFF..."
pixi run python ../convert_cmyk.py results/Fig1_Parashorea_Geographic_Map_Final.png results/Fig1_Parashorea_Geographic_Map_Final_CMYK.tiff
echo "Done! Results in results/"
