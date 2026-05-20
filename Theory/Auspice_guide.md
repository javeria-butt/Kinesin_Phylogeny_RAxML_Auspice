# Auspice Interactive Viewer — Guide

[Auspice](https://docs.nextstrain.org/projects/auspice/en/stable/) is the interactive phylogenetic visualisation component of Nextstrain. This guide explains how to view the kinesin trees interactively.

---

## Quick Start (No Installation)

1. Generate the Auspice JSON (already done if you ran the notebooks):
   - `auspice/kinesin_auspice.json` — full kinesin superfamily
   - `auspice/kinesin14_auspice.json` — Kinesin-14 family

2. Go to **[https://auspice.us](https://auspice.us)**

3. Drag and drop either JSON file onto the page.

The tree loads instantly in your browser — no account or installation required.

---

## What You Can Do in Auspice

| Feature | How |
|---------|-----|
| Colour by family / species / direction | Dropdown menu top-left |
| Zoom into a clade | Click any internal node |
| Search for a sequence | Search box top-right |
| Toggle rectangular / radial layout | Layout button |
| Show / hide branch labels | Options panel |
| Export as SVG | Download button |

---

## Local Auspice Installation (Optional)

If you want to run Auspice locally (e.g. for offline use):

```bash
# Install Node.js first: https://nodejs.org
npm install --global auspice

# Serve the JSON files
auspice view --datasetDir auspice/
```

Then open `http://localhost:4000` in your browser and select the dataset from the dropdown.

---

## Colour Keys

### Full Kinesin Superfamily (`kinesin_auspice.json`)

| Colour by | Options |
|-----------|---------|
| `family` | Kinesin-1 through Kinesin-14 |
| `species` | H. sapiens, D. melanogaster, S. cerevisiae, … |
| `directionality` | Plus-end / Minus-end |
| `gene` | KIF5A, KHC, NCD, KAR3, … |

### Kinesin-14 Family (`kinesin14_auspice.json`)

| Colour by | Options |
|-----------|---------|
| `directionality` | Minus-end (confirmed) / Unknown |
| `species` | H. sapiens, D. melanogaster, S. cerevisiae, A. thaliana, C. elegans, C. griseus |
| `gene` | KIFC1, KIFC2, KIFC3, NCD, KAR3, KCBP, … |

---

## Troubleshooting

**Tree does not load**
- Make sure the file ends in `.json`
- Check that you generated it with `auspice_json_builder.py` or the notebook
- Try a different browser (Chrome / Firefox recommended)

**No colours appear**
- Use the "Color By" dropdown to select an attribute

**Branch labels missing**
- Enable "Show branch labels" in the Options panel
