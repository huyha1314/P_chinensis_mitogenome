#!/usr/bin/env bash
# ==============================================================
# Figure 7: RNA Editing & Positive Selection Landscape
# ==============================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p results

echo "Running Figure 7 R script..."
micromamba run -n visualiz Rscript Fig7.R
echo "Done! Outputs saved to results/"
