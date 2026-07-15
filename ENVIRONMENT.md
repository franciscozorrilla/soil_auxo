# Software environment

This file lists the software needed to reproduce the analyses. External tool versions are
those stated in the manuscript Methods. Figures can be regenerated from the committed
intermediate data (`data/figure_data/`) with only R and Python; reproducing the pipeline from
raw reads additionally needs the tools in the *Metagenomic pipeline* table and an HPC cluster.

## Languages

- **Python 3.9** — the Jupyter notebooks in `code/` were authored under Python 3.9.13.
- **R ≥ 4.1** — the analysis notebooks and scripts under `analysis/`.

(R package and Python package versions were the current CRAN / PyPI releases at the time of
analysis in 2025; the pinned external-tool versions below are the authoritative record.)

## Python packages

Figure rendering and image analysis (`code/plots.ipynb`, `code/Image processing …ipynb`) — see
[`requirements.txt`](requirements.txt):

```
matplotlib  numpy  pandas  Pillow  pycirclize  scikit-image  scipy
```

GEM auxotrophy prediction (`analysis/02_gem_auxotrophy/get_auxo.py`):

```
cobra  reframed        # plus a linear-programming solver (e.g. CPLEX) at runtime
```

## R packages

```
tidyverse   ape        caret      ggalluvial  gggenes    gggenomes   ggnewscale
ggpubr      ggrepel     ggtree     ggtreeExtra gridExtra  MASS        phylolm
rstatix     factoextra
```

`ggtree` / `ggtreeExtra` are from Bioconductor; the rest are on CRAN. The phylogenetic
logistic regression uses `phyloglm()` from **phylolm**; the Box–Cox transform uses **MASS**.

## External bioinformatics tools (versions from the Methods)

| Step | Tool | Version |
|---|---|---|
| Quality filtering | fastp | 0.20.1 |
| Assembly | MEGAHIT | 1.2.9 |
| Read mapping | bwa | 0.7.17 |
| Isolate assembly | shovill | 1.1.0 |
| Binning | MaxBin2 / MetaBAT / CONCOCT | 2.2.7 / 2.15 / 1.1.0 |
| Bin refinement | metaWRAP | 1.3.2 |
| Dereplication | dRep | 3.0.0 |
| Taxonomy | GTDB-Tk | 2.0.0 (reference r207) |
| Functional annotation | eggNOG-mapper (eggNOG DB 5.0.2) | 2.1.9 |
| Mobile elements | geNomad | 1.5.0 |
| Insertion sequences | DIAMOND (vs the ISfinder database) | 2.0.15 |
| Metabolic models | CarveMe | 1.5.1 |
| Auxotrophy simulation | ReFramed | 1.4 |
| Cross-feeding | SMETANA (`--detailed`) | 1.1.0 |
| Pathway completeness | KEMET | — |
| Genome completeness | CheckM | — |
| Tree visualisation | iTOL | v5 |
| Pipeline orchestration | metaGEM (Snakemake) | — |

The reconstructed genome-scale metabolic models are deposited in `data/models/`, so the
downstream analyses do not require re-running CarveMe or the upstream pipeline.
