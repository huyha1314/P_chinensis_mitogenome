#!/bin/bash
set -e

echo "Running Figure 3 plotting script..."
micromamba run -n visualiz python plot_fig3.py

echo "Done!"

