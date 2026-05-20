# Kinesin Motor Protein Phylogenetic Analysis

<div align="center">

![Python](https://img.shields.io/badge/Python-3.8%2B-blue?style=for-the-badge&logo=python)
![Biopython](https://img.shields.io/badge/Biopython-1.81-green?style=for-the-badge)
![Jupyter](https://img.shields.io/badge/Jupyter-Notebook-orange?style=for-the-badge&logo=jupyter)
![License](https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge)
![Colab](https://img.shields.io/badge/Google-Colab-red?style=for-the-badge&logo=googlecolab)

**A comprehensive phylogenetic analysis of Kinesin motor proteins and the Kinesin-14 family using Biopython, MUSCLE, RAxML, and interactive Auspice visualizations.**


</div>

---

## Table of Contents

- [Background](#-background)
- [Repository Structure](#-repository-structure)
- [Installation](#-installation)
- [Notebooks](#-notebooks)
- [Data](#-data)
- [Results & Figures](#-results--figures)
- [Methods](#-methods)
- [References](#-references)

---

## Background

### Kinesin Motor Proteins

Kinesin motor proteins are specialized molecular motors in eukaryotic cells that "walk" along microtubules, acting as intracellular transport vehicles powered by ATP hydrolysis. They move cargo — such as organelles, vesicles, and mRNA — towards the **plus end** of microtubules (anterograde transport), essential for:

- Intracellular organization
- Cell division
- Axonal transport in neurons

The kinesin superfamily is classified into **14+ families (Kinesin-1 through Kinesin-14)** based on motor domain position and phylogenetic grouping.

### The Kinesin-14 Family

The **Kinesin-14** family (formerly C-terminal motor) proteins possess a unique C-terminal motor domain. At least four members (DmNcd, ScKAR3, CgCHO2, AtKCBP) are confirmed **minus-end directed motors**, making them functionally distinct from all other kinesin families.

### Phylogenetic Approach

| Tool | Purpose |
|------|---------|
| **Biopython** | Sequence fetching, parsing, alignment I/O |
| **MUSCLE** | Multiple sequence alignment |
| **RAxML / FastTree** | Maximum likelihood tree inference |
| **ETE3 / Matplotlib** | Static tree rendering and figure generation |
| **Auspice (Nextstrain)** | Interactive phylogeny visualization |

---

## Repository Structure

```
Kinesin_Phylogeny_RAxML_Auspice/
│
├── 📓 Notebooks/
│   ├── 01_Kinesin_Motor_Proteins_Phylogeny.ipynb
│   └── 02_Kinesin14_Family_Phylogeny.ipynb
│
├── 🧬 Data/
│   ├── Sequences/
│   │   ├── kinesin_all_families.fasta
│   │   ├── kinesin14_family.fasta
│   │   └── accession_numbers.txt
│   └── Alignments/
│       ├── kinesin_aligned.fasta
│       └── kinesin14_aligned.fasta
│
├── 📊 Figures/
│   ├── kinesin_phylotree.png
│   ├── kinesin14_phylotree.png
│   ├── kinesin_circular_tree.png
│   ├── sequence_length_distribution.png
│   └── bootstrap_support_heatmap.png
│
├── 🔧 Scripts/
│   ├── fetch_sequences.py
│   ├── run_alignment.py
│   ├── run_raxml.py
│   ├── plot_tree.py
│   └── auspice_json_builder.py
│
├── 🌐 Auspice/
│   ├── kinesin_auspice.json
│   └── kinesin14_auspice.json
│
├── 📖 Theory/
│   ├── methods.md
│   └── auspice_guide.md
│
└── README.md
```

---

## Installation

### Option 1: Google Colab (No Setup Required)

Click the **"Open in Colab"** badges above. All dependencies install automatically inside the notebooks.

### Option 2: Local Setup

```bash
git clone https://github.com/YOUR_USERNAME/kinesin-phylogenetics.git
cd kinesin-phylogenetics
```

#### Conda (recommended)
```bash
conda env create -f environment.yml
conda activate kinesin-phylo
jupyter notebook
```

#### pip
```bash
pip install -r requirements.txt
```

#### System tools (Ubuntu/Debian)
```bash
sudo apt-get install muscle raxml mafft fasttree
```

#### System tools (macOS)
```bash
brew install muscle raxml mafft fasttree
```

---

## Notebooks

### Notebook 1 — Kinesin Motor Proteins Phylogeny
`notebooks/01_Kinesin_Motor_Proteins_Phylogeny.ipynb`

| Section | Description |
|---------|-------------|
| 1. Setup & Dependencies | Install and import all libraries |
| 2. Sequence Retrieval | Fetch Kinesin-1 through Kinesin-14 from NCBI |
| 3. Sequence QC | Length distributions, composition stats |
| 4. Multiple Alignment | MUSCLE alignment with visualization |
| 5. Tree Inference | RAxML / FastTree ML tree |
| 6. Tree Visualization | Matplotlib + ETE3 figures |
| 7. Auspice Export | Build interactive JSON for Nextstrain |
| 8. Bootstrap Analysis | Support values mapped to branches |

### Notebook 2 — Kinesin-14 Family Phylogeny
`notebooks/02_Kinesin14_Family_Phylogeny.ipynb`

| Section | Description |
|---------|-------------|
| 1. Setup & Dependencies | Install and import all libraries |
| 2. Kinesin-14 Sequences | Fetch DmNcd, ScKAR3, CgCHO2, AtKCBP, HsKIFC2/3 |
| 3. Motor Domain Alignment | MUSCLE + conserved block extraction |
| 4. ML Tree (RAxML-NG) | 100 bootstrap replicates |
| 5. Directionality Annotation | Highlight minus-end vs plus-end clades |
| 6. Species Mapping | Organism-level color annotations |
| 7. Auspice Export | Interactive JSON visualization |

---

## Data

### Kinesin Families Analyzed

| Family | Key Members | Directionality | Function |
|--------|------------|----------------|----------|
| Kinesin-1 | KHC, KIF5A/B/C | Plus-end | Organelle/vesicle transport |
| Kinesin-2 | KRP85/95, KIF3A/B | Plus-end | Cilia, flagella (IFT) |
| Kinesin-3 | Unc104, KIF1A/B/C | Plus-end | Synaptic vesicle transport |
| Kinesin-4 | ChrKin, KIF4A/B | Plus-end | Chromosome positioning |
| Kinesin-5 | BimC, Eg5/KIF11 | Plus-end | Bipolar spindle assembly |
| Kinesin-6 | MKLP1, KIF23 | Plus-end | Cytokinesis |
| Kinesin-7 | CENP-E, KIF10 | Plus-end | Kinetochore–MT attachment |
| Kinesin-8 | Kip3, KIF18A/B | Plus-end | MT depolymerization |
| Kinesin-13 | MCAK, KIF2A/B/C | Depolymerizer | MT catastrophe |
| **Kinesin-14** | **NCD, KAR3, KIFC1/2/3** | **Minus-end** | **Spindle organization** |

### Sequence Sources
- **NCBI Protein** — primary sequence database
- **UniProt/SwissProt** — reviewed, curated entries
- **Duke Kinesin Site** — https://sites.duke.edu/kinesin/

### Species Abbreviations

| Prefix | Species |
|--------|---------|
| Hs | *Homo sapiens* |
| Dm | *Drosophila melanogaster* |
| Sc | *Saccharomyces cerevisiae* |
| At | *Arabidopsis thaliana* |
| Ce | *Caenorhabditis elegans* |
| Cg | *Cricetulus griseus* |

---

## Results & Figures

### Figure 1 — Full Kinesin Superfamily Tree
Maximum likelihood phylogenetic tree of the kinesin superfamily. Families are color-coded by clade. Bootstrap support values (100 replicates) shown at internal nodes.

### Figure 2 — Kinesin-14 Family Tree
Resolved phylogeny of the Kinesin-14 family. Red branches = confirmed minus-end directed motors. Bootstrap values shown at nodes.

### Figure 3 — Circular (Radial) Tree
Radial layout emphasizing relative evolutionary distances between all kinesin families.

### Figure 4 — Sequence Length Distribution
Distribution of sequence lengths across all kinesin families, grouped by family.

---

## Methods

### 1. Sequence Retrieval
Protein sequences retrieved from NCBI using Biopython `Entrez`:
```python
Entrez.esearch(db="protein", term="kinesin[Title] AND motor[Title]")
```

### 2. Multiple Sequence Alignment
Aligned with **MUSCLE v5**:
```bash
muscle -align kinesin_all_families.fasta -output kinesin_aligned.fasta
```

### 3. Tree Inference
Maximum likelihood tree with **RAxML**, LG+G4 model, 100 bootstrap replicates:
```bash
raxmlHPC -f a -m PROTGAMMALG -p 12345 -x 12345 -# 100 \
         -s kinesin_aligned.fasta -n kinesin_tree
```

### 4. Tree Visualization
Trees rendered with **ETE3** and **Matplotlib**. Interactive export to **Auspice JSON** format for Nextstrain viewer.

See [`Theory/methods.md`](Theory/methods.md) for complete details.

---

## References

1. Lawrence, C.J. et al. (2004). A standardized kinesin nomenclature. *J. Cell Biol.* 167(1):19–22.
2. Dagenbach, E.M. & Endow, S.A. (2004). A new kinesin tree. *J. Cell Sci.* 117(1):3–7.
3. Cock, P.J.A. et al. (2009). Biopython: freely available Python tools for computational molecular biology. *Bioinformatics* 25(11):1422–1423.
4. Stamatakis, A. (2014). RAxML version 8. *Bioinformatics* 30(9):1312–1313.
5. Edgar, R.C. (2004). MUSCLE: multiple sequence alignment with high accuracy. *Nucleic Acids Res.* 32(5):1792–1797.
6. Hadfield, J. et al. (2018). Nextstrain: real-time tracking of pathogen evolution. *Bioinformatics* 34(23):4121–4123.
7. Duke Kinesin Site: https://sites.duke.edu/kinesin/

---


