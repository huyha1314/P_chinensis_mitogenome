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
micromamba run -n visualiz Rscript fig5.R
echo "Done! Outputs saved to results/"
