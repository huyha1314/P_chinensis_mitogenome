#!/usr/bin/env bash
# ==============================================================
# Figure S1: ONT Long Reads Spanning Recombination Breakpoint
# ==============================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p results

echo "Running Figure S1 Python script..."
pixi run python figS1.py
echo "Converting to publication-ready CMYK TIFF..."
pixi run python ../convert_cmyk.py results/FigureS1.png results/FigureS1_CMYK.tiff
echo "Done! Outputs saved to results/"
