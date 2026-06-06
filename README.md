# Mitogenome Visualization Pipeline for *Parashorea chinensis*

[![Python](https://img.shields.io/badge/Python-3.14.4-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.5.3-blue.svg)](https://www.r-project.org/)
[![Circos](https://img.shields.io/badge/Circos-0.69.9-green.svg)](http://circos.ca/)
[![Bioconda](https://img.shields.io/badge/recipe-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the official visualization pipeline and source code for the genomic figures presented in the paper:  
**"Structural Hyper-Variation and Adaptive Evolution in the Bipartite Mitogenome of Parashorea chinensis"**

The pipeline is designed for high reproducibility, utilizing pinned conda environments and automated execution scripts to generate publication-quality figures (TIFF/PNG/PDF/SVG).

---

## 📂 Repository Structure

The repository is organized by figure panels. Each directory includes its own `data/` (inputs), `results/` (outputs), and execution scripts.

```text
.
├── env/                    # Environment configuration
│   ├── env.yml             # Pinned dependencies (Conda + Pip)
│   └── create_env.sh       # Environment creation helper
├── fig1/                   # Figure 1: Geographic Distribution Map (R)
├── fig2/                   # Figure 2: Assembly Graphs (GFA files)
├── fig3/                   # Figure 3: Circular Mitogenome Maps (pycirclize)
├── fig4/                   # Figure 4: Linear Synteny Visualization (Python)
├── fig5/                   # Figure 5: Tanglegrams / Cophylogeny (R)
├── fig6/                   # Figure 6: Phylogenetic Tree (R)
├── fig7/                   # Figure 7: Intergenomic DNA Transfers (Circos)
├── fig8/                   # Figure 8: RNA Editing & Positive Selection (R)
└── figS1/                  # Figure S1: Recombination Breakpoint Proof (pysam)
```

---

## ⚙️ Installation & Setup

We use **Pixi** (recommended) or **Micromamba** (or Conda) to manage the specialized bioinformatics environment. The environment pins all major dependencies to specific versions used in the study to ensure long-term reproducibility.

### Option A: Using Pixi (Recommended, Automated)
Simply run any of the execution scripts or commands. Pixi will automatically download, compile, and manage all R, Python, and bioinformatics dependencies (like Circos, ggtree, etc.) in a local environment:
```bash
# Run a specific script using the pixi environment:
pixi run Rscript fig6/Fig6.R
```

### Option B: Using Micromamba / Conda
1. **Clone the repository:**
   ```bash
   git clone https://github.com/username/P_chinensis_mitogenome_vis.git
   cd P_chinensis_mitogenome_vis
   ```

2. **Build the environment:**
   ```bash
   cd env
   ./create_env.sh
   micromamba activate visualiz
   ```
   *Note: This will install Python 3.10, R 4.3.3, and all necessary bioinformatics libraries including Circos and pysam.*

---

## 🚀 Execution Guide

Each figure can be regenerated using the provided `run.sh` scripts.

### Figure 1: Geographic Map
Generates the sampling location and distribution map.
```bash
cd fig1 && ./run.sh
```

### Figure 3: Circular Genomic Maps
Generates OGDRAW-style circular maps with GC content, GC skew, and gene orientations.
```bash
cd fig3 && ./run.sh
```

### Figure 4: Linear Synteny Map
Performs block-based synteny alignment between multiple mitogenome accessions.
```bash
cd fig4 && ./run.sh
```

### Figure 5: Tanglegrams (Phylogenetic Comparison)
Generates high-resolution cophylogeny plots (Nuclear vs. CP and Mito vs. CP) with curved links.
```bash
cd fig5 && ./run.sh
```

### Figure 6: Phylogenetic Tree
Plots the deep-time and late Pleistocene isolation phylogenetic trees.
```bash
cd fig6 && ./run.sh
```

### Figure 7: Intergenomic DNA Transfers
Runs the Circos engine for complex DNA transfer modeling across organelles.
```bash
cd fig7 && ./run.sh
```
*Or execute step-by-step:*
```bash
cd fig7
./Fig7A.sh
./Fig7B.sh
pixi run python combine_fig7.py
```

### Figure 8: RNA Editing & Positive Selection
Stitches RNA editing bar charts with phylogenetic trees showing episodic selection sites.
```bash
cd fig8 && ./run.sh
```

### Figure S1: Recombination Proof
Visualizes long-read (ONT) support for structural variants spanning recombination breakpoints.
```bash
cd figS1 && ./run.sh
```

---

## 🎨 Publication & Journal Submission Requirements

All color figures are fully compliant with standard journal print requirements:
*   **Color Space**: Set up as **CMYK** (converted from RGB to prevent color shifting during printing).
*   **Dimensions**: Sized to fit exactly within the full text width of **17.0 cm** (6.69 inches) or single column of **8.0 cm** (3.15 inches), with height not exceeding **22.5 cm**.
*   **Resolution**: Exported as high-resolution (**600 DPI**) TIFF files with lossless LZW compression.
*   **Scale Marks**: Coordinate scales, phylogenetic branch lengths, physical maps, and read depths are clearly labeled with appropriate units or scale bars.

### Generate CMYK TIFF Figures
To convert the generated figures to publication-ready CMYK TIFF format, run:
```bash
pixi run python convert_cmyk.py
```
This script will process all figure files and output `*_CMYK.tiff` files under the respective `results/` folders.
