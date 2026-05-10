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
├── fig6/                   # Figure 6: Intergenomic DNA Transfers (Circos)
├── fig7/                   # Figure 7: RNA Editing & Positive Selection (R)
└── figS1/                  # Figure S1: Recombination Breakpoint Proof (pysam)
```

---

## ⚙️ Installation & Setup

We use **Micromamba** (or Conda) to manage the specialized bioinformatics environment. The environment pins all major dependencies to specific versions used in the study to ensure long-term reproducibility.

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
   *Note: This will install Python 3.14.4, R 4.5.3, and all necessary bioinformatics libraries including Circos and pysam.*

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

### Figure 6: Intergenomic DNA Transfers
Runs the Circos engine for complex DNA transfer modeling across organelles.
```bash
cd fig6
./Fig6A.sh
./Fig6B.sh
python combine_fig6.py
```

### Figure 7: RNA Editing & Positive Selection
Stitches RNA editing bar charts with phylogenetic trees showing episodic selection sites.
```bash
cd fig7 && ./run.sh
```

### Figure S1: Recombination Proof
Visualizes long-read (ONT) support for structural variants spanning recombination breakpoints.
```bash
cd figS1 && ./run.sh
```


