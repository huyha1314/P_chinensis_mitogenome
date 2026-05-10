#!/usr/bin/env bash
# ==============================================================
# Figure S1: ONT Long Reads Spanning Recombination Breakpoint
# ==============================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p results

echo "Running Figure S1 Python script..."
micromamba run -n visualiz python figS1.py
echo "Done! Outputs saved to results/"
