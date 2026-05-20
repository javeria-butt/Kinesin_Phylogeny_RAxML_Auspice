# Methods

## 1. Sequence Retrieval

Protein sequences were retrieved from the NCBI Protein database using Biopython's `Entrez` module. Representative members of each kinesin family were selected based on the phylogenetic classification published at the [Duke Kinesin Site](https://sites.duke.edu/kinesin/) (Dagenbach & Endow, 2004).

```python
from Bio import Entrez, SeqIO
Entrez.email = "user@example.com"
handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
```

Requests were rate-limited to ≤3 per second in accordance with NCBI usage guidelines.

---

## 2. Multiple Sequence Alignment

Sequences were aligned using **MUSCLE v5** with default parameters:

```bash
muscle -align input.fasta -output aligned.fasta
```

If MUSCLE was unavailable, **MAFFT** (`--auto` mode) was used as a fallback:

```bash
mafft --auto input.fasta > aligned.fasta
```

Alignment quality was assessed by inspecting conserved blocks around the kinesin motor domain (~350 aa), which is expected to align reliably across all family members.

---

## 3. Phylogenetic Tree Inference

### Primary: RAxML (Maximum Likelihood)

Maximum likelihood trees were inferred using **RAxML v8** with the LG substitution matrix and a discrete Gamma model of rate heterogeneity (PROTGAMMALG). Bootstrap support was assessed with 100 rapid bootstrap replicates (`-f a`):

```bash
raxmlHPC -f a \
  -m PROTGAMMALG \
  -p 42 -x 42 \
  -# 100 \
  -s aligned.fasta \
  -n kinesin \
  -w /path/to/results/
```

The best-scoring ML tree with bootstrap values mapped was used for all subsequent analyses (`RAxML_bipartitions.*`).

### Fallback: FastTree

When RAxML was unavailable, **FastTree** was used with the LG model and Gamma rate correction:

```bash
FastTree -lg -gamma -quiet aligned.fasta > tree.nwk
```

### Fallback: Neighbour-Joining (Biopython)

As a last resort, a Neighbour-Joining tree was computed using Biopython with BLOSUM62 distances:

```python
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
calc = DistanceCalculator("blosum62")
dm   = calc.get_distance(alignment)
ctor = DistanceTreeConstructor(calc, "nj")
tree = ctor.build_tree(alignment)
```

NJ trees do not carry bootstrap support and should be treated as preliminary.

---

## 4. Tree Visualization

Trees were rendered using Biopython's `Phylo.draw()` with custom colour assignments per kinesin family (full superfamily analysis) or per motor directionality and species (Kinesin-14 analysis).

Circular (radial) trees were drawn using Matplotlib's polar projection.

All figures were saved at 150 DPI in PNG format.

---

## 5. Interactive Visualization — Auspice

Trees were exported to **Auspice v2 JSON** format for interactive viewing at [https://auspice.us](https://auspice.us). Node attributes include:

- `family` — kinesin family assignment
- `species` — source organism
- `gene` — gene name
- `directionality` — plus-end / minus-end (where known)
- `bootstrap` — bootstrap support value (if available)

To view:
1. Run `python scripts/auspice_json_builder.py` to generate the JSON
2. Go to [https://auspice.us](https://auspice.us)
3. Drag and drop the JSON file onto the page

---

## 6. Bootstrap Interpretation

Branch support thresholds used:

| Bootstrap value | Interpretation |
|-----------------|---------------|
| ≥ 90%           | Strongly supported |
| 70–89%          | Well supported |
| 50–69%          | Moderately supported |
| < 50%           | Poorly supported |

Only nodes with ≥ 70% bootstrap support are considered reliable for biological interpretation.

---

## 7. Software Versions

| Software | Version | Reference |
|----------|---------|-----------|
| Python   | 3.10    | — |
| Biopython | 1.81  | Cock et al. 2009 |
| MUSCLE   | 5.1     | Edgar 2004 |
| RAxML    | 8.2.12  | Stamatakis 2014 |
| FastTree | 2.1.11  | Price et al. 2010 |
| Matplotlib | 3.7   | Hunter 2007 |
| NumPy    | 1.24    | Harris et al. 2020 |

---

## References

- Cock, P.J.A. et al. (2009). Biopython. *Bioinformatics* 25:1422–1423.
- Dagenbach, E.M. & Endow, S.A. (2004). A new kinesin tree. *J. Cell Sci.* 117:3–7.
- Edgar, R.C. (2004). MUSCLE. *Nucleic Acids Res.* 32:1792–1797.
- Lawrence, C.J. et al. (2004). A standardized kinesin nomenclature. *J. Cell Biol.* 167:19–22.
- Price, M.N. et al. (2010). FastTree 2. *PLoS ONE* 5:e9490.
- Stamatakis, A. (2014). RAxML version 8. *Bioinformatics* 30:1312–1313.
