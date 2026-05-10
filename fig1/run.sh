#!/usr/bin/env bash
set -e
mkdir -p results
echo "Generating Geographic Map (Figure 1)..."
micromamba run -n visualiz Rscript Fig1.R
echo "Done! Results in results/"
