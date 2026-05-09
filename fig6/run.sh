#!/bin/bash
set -e

echo "Generating Fig6A..."
bash Fig6A.sh

echo "Generating Fig6B..."
bash Fig6B.sh

echo "Combining figures..."
micromamba run -n visualiz python combine_fig6.py

echo "Done!"