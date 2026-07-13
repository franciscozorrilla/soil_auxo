# 🧫 soil_auxo

Code and data for the publication *Obligate cross-feeding of metabolites is common in soil microbial communities*.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13847512.svg)](https://doi.org/10.5281/zenodo.13847512)

## 🧬 Abstract

Ecological theory predicts that obligate cross-feeding of metabolites is common in microbial communities. However, systematic studies testing this prediction in natural microbial communities like soil are missing, thus hindering the understanding of these ecologically important ecosystems. Here, we address this gap by analysing 6,931 bacterial isolates from 27 soil microbial communities sampled in Germany. We find that the growth of up to 50% of community members depended on supplementation with amino acids, vitamins, or nucleotides. For 73% of isolates, supplementation with multiple amino acids was necessary. Genomic analysis of 62 strains revealed that accumulation of insertion sequences and gene loss was associated with the observed auxotrophic phenotypes. Genome-scale metabolic models, computational analyses, and cocultivation experiments suggest other co-occurring genotypes complement the metabolic needs of auxotrophs, thus enabling their growth and survival. Together, our results indicate that soil bacteria likely exist within integrated metabolic networks rather than as metabolically autonomous units.

## 🗺️ Finding the code and data for a figure

Every figure, supplementary figure, table, and analysis is mapped to the exact code and data that produce it:

- **[REPO_MAP.md](REPO_MAP.md)** — human-readable map (per-figure tables, contributions, reproduction values).
- **[MANIFEST.tsv](MANIFEST.tsv)** — the same map, machine-readable.
- Each figure folder under **`figures/`** has its own `README.md`.

**Two code layers, both provided.** For the figures produced computationally, the **R** notebooks in `analysis/` compute the numbers and intermediate tables, and the main published panels (Fig 1b/c/d, 2b, 3, 4, 6) and Supplementary Fig. 8 were then re-rendered in **Python** (`code/plots.ipynb`); the other supplementary figures are the R renders. The experimental figures (coculture assays, LC-MS/MS, microscopy) are the wet-lab authors' work — see [Data availability](#-data-availability).

## 📁 Repository layout

```
soil_auxo/
├── analysis/            R/Python analysis code (computation behind every figure), in pipeline order
│   ├── 01_metagem_pipeline/      metaGEM Snakefile + config.yaml (raw reads → MAGs → GEMs)
│   ├── 02_gem_auxotrophy/        get_auxo.py → auxotrophy predictions from the GEMs
│   ├── 03_phylogenetics/         phylogenetic.Rmd → tree, phyloglm, MGE/IS (Supp Figs 4,5,11; Table 2)
│   ├── 04_genome_features/       metabolic_modeling_soil.Rmd → Fig 4, Supp Figs 5–7, 260-gene set
│   ├── 05_smetana_crossfeeding/  SMETANA cross-feeding (data under data/smetana/)
│   ├── 06_saturation/            sat_sim/sat_sim_notebook.Rmd → complementarity vs community size (Fig 6b)
│   └── 07_community_composition/ assembly/binning/qfilter/drep/compositionVis.R (Supp Fig 3)
├── code/                Python figure renders (plots.ipynb, image-analysis notebook)
├── figures/             one folder per figure (fig1–6, supp_fig_01–11), each with a README + outputs
├── data/
│   ├── models/          genome-scale metabolic models (SBML, 102 genomes)
│   ├── genomes/         assembled nucleotide + protein sequences, 16S
│   ├── annotations/     eggNOG, ISfinder, geNomad, MAG GFF (gzip-compressed)
│   ├── figure_data/     intermediate tables consumed directly by the plotting code
│   ├── gem_curation/    GEM curation inputs (GPRs, media, memote)
│   ├── qc/              assembly/binning statistics
│   └── smetana/         SMETANA cross-feeding output
├── supplementary_tables/  Supplementary_Tables.xlsx + 260-gene table
└── archive/             superseded drafts and exploratory analyses (not in the paper)
```

## ▶️ Reproducing the analysis

Raw sequences are in the European Nucleotide Archive under [**PRJEB80563**](https://www.ebi.ac.uk/ena/browser/view/PRJEB80563). To reproduce the metagenomic pipeline from raw reads you need an HPC cluster and [metaGEM](https://github.com/franciscozorrilla/metaGEM); to reproduce the figures from the committed intermediate data you only need R and Python. The full software list (R and Python packages, and the external tool versions from the Methods) is in [ENVIRONMENT.md](ENVIRONMENT.md); Python dependencies are in [requirements.txt](requirements.txt).

**Run all scripts from the repository root.** Reads resolve to `data/…`; regenerated figures are written to `outputs/` (gitignored), while the committed published versions live in `figures/`. Compressed inputs (`*.gz`) are read transparently. Stochastic steps are seeded (`set.seed(42)`) so the bootstrapped confidence intervals and complementarity curves are reproducible.

Pinned reproduction checks (from the committed data): Fig 3a Spearman ρ = **0.360**; Fig 4 medians **934 / 880** (genes·Mbp⁻¹) and **227 / 5.1** (IS·Mbp⁻¹); Fig 6b complementarity **29 / 65 / 87 / 95 %**.

## 📦 Data availability

The intermediate tables under `data/figure_data/` are sufficient to regenerate the computational figures (Fig. 3, 4, 6b; Supplementary Figs 3–7, 11) without the raw data. Large annotation and cross-feeding outputs are gzip-compressed in place and read transparently.

Several figures report experimental measurements (coculture assays for Fig. 5a/b and Supplementary Fig. 8; LC-MS/MS for Supplementary Fig. 9; microscopy) or were contributed by co-authors (Fig. 1a, 2a, 6a; Supplementary Figs 1, 2, 10). Their underlying measurements are provided as source data with the paper or are available from the corresponding authors; see [REPO_MAP.md](REPO_MAP.md) for the per-figure breakdown.

## 📄 Licence

Code is released under the MIT Licence ([LICENSE](LICENSE)); data under CC-BY-4.0 ([DATA_LICENSE](DATA_LICENSE)).

## 🖋️ Citation

If you use this resource, please cite:

> Obligate cross-feeding of metabolites is common in soil microbial communities. Ghada Yousif*, Francisco Zorrilla*, Swagatika Dash, Leonardo Oña, Aditi Shekhar, Samir Giri, Rui Guan, Sharvari Harshe, Michael Itermann, Daphne Welter, Vladimir Benes, Kiran R. Patil, Christian Kost. bioRxiv 2025.01.29.635426; doi: https://doi.org/10.1101/2025.01.29.635426
