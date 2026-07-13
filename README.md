# 🧫 soil_auxo

Code and data for the publication *Obligate cross-feeding of metabolites is common in soil microbial communities*.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13847512.svg)](https://doi.org/10.5281/zenodo.13847512)

This README maps **every figure, supplementary figure, table, and analysis** to the code that produces it, the data it consumes, and who produced it, so a reader can locate the code and data behind any panel without guesswork. [`MANIFEST.tsv`](MANIFEST.tsv) is the same map in machine-readable form, and each folder under `figures/` has its own README.

## Contents

- [Abstract](#-abstract)
- [Contributions](#-contributions)
- [Repository layout](#-repository-layout)
- [Reproducing the analysis](#️-reproducing-the-analysis)
- [Main figures](#-main-figures)
- [Supplementary figures](#-supplementary-figures)
- [Supplementary tables](#-supplementary-tables)
- [Analysis pipeline](#-analysis-pipeline-upstream-of-the-figures)
- [Data availability](#-data-availability)
- [Not in the paper → archive/](#-not-in-the-paper--archive)
- [Licence & citation](#-licence)

## 🧬 Abstract

Ecological theory predicts that obligate cross-feeding of metabolites is common in microbial communities. However, systematic studies testing this prediction in natural microbial communities like soil are missing, thus hindering the understanding of these ecologically important ecosystems. Here, we address this gap by analysing 6,931 bacterial isolates from 27 soil microbial communities sampled in Germany. We find that the growth of up to 50% of community members depended on supplementation with amino acids, vitamins, or nucleotides. For 73% of isolates, supplementation with multiple amino acids was necessary. Genomic analysis of 62 strains revealed that accumulation of insertion sequences and gene loss was associated with the observed auxotrophic phenotypes. Genome-scale metabolic models, computational analyses, and cocultivation experiments suggest other co-occurring genotypes complement the metabolic needs of auxotrophs, thus enabling their growth and survival. Together, our results indicate that soil bacteria likely exist within integrated metabolic networks rather than as metabolically autonomous units.

**Numbering.** This repository uses the **Supplementary Figure 1–11 / Supplementary Table 1–4** numbering of the manuscript.

**Two code layers.** For the figures produced computationally, the **R** notebooks in `analysis/` compute the numbers and intermediate tables, and the main published panels (Fig. 1b/c/d, 2b, 3, 4, 6) and Supplementary Fig. 8 were re-rendered in **Python** (`code/plots.ipynb`, by L. Oña) from those computed values; the remaining supplementary figures are the R renders. Both layers are provided. The experimental figures (coculture assays, LC-MS/MS, microscopy) are the wet-lab authors' work; their measurements are provided as source data with the paper (see [Data availability](#-data-availability)).

## 👥 Contributions

- **Figures 1 & 2** (isolation screen, auxotrophy prevalence, 16S phylogenetic distribution, community richness): experimental work led by **G. Yousif**; panels 1b/c/d and 2b were rendered in Python by **L. Oña**.
- **Figures 3 & 4** (GEM predictions, genome coding density and insertion sequences): analysis by **F. Zorrilla** (R), rendered in Python by **L. Oña**.
- **Figure 5** (pairwise cocultures, cross-feeding networks): **L. Oña / G. Yousif**. Panel **5c** was assembled by **G. Yousif** from SMETANA output produced by **F. Zorrilla**.
- **Figure 6**: panel **6a** (27-community complementarity) by **L. Oña**; panel **6b** (58 focal isolates) by **F. Zorrilla**, rendered in Python by **L. Oña**.
- **Supplementary Figures 3, 4, 5, 6, 7, 11**: **F. Zorrilla** (R). Supp Fig **1** (map) and **9** (LC-MS) by **G. Yousif**; Supp Fig **8** by **L. Oña**; Supp Fig **2** (metagenome composition) and **10** (SMETANA network) by **G. Yousif** (Supp Fig 10 from F. Zorrilla's SMETANA output).
- Everything under `analysis/` (metaGEM pipeline, GEM auxotrophy prediction, phylogenetics/phyloglm, genome-feature analysis, SMETANA runs, saturation, assembly/binning QC) is **F. Zorrilla's**.

## 📁 Repository layout

```
soil_auxo/
├── README.md                 this file — figure/analysis → code+data map + quickstart
├── MANIFEST.tsv              machine-readable version of the map
├── ENVIRONMENT.md            software environment (R/Python packages, tool versions)
├── requirements.txt          Python dependencies
├── LICENSE                   MIT (code)  ·  DATA_LICENSE  CC-BY-4.0 (data)  ·  CITATION.bib
├── analysis/                 analysis code, in pipeline order
│   ├── 00_load_data.R            shared data loader sourced by the notebooks below
│   ├── 01_metagem_pipeline/      metaGEM Snakefile + config.yaml (raw reads → MAGs → GEMs)
│   ├── 02_gem_auxotrophy/        get_auxo.py → auxotrophy predictions from GEMs
│   ├── 03_phylogenetics/         phylogenetic.Rmd → tree, phyloglm, MGE/IS (Supp Figs 4,5,11; Table 2)
│   ├── 04_genome_features/       metabolic_modeling_soil.Rmd → Fig 4, eggNOG, geNomad, gene counts, 260-gene set
│   ├── 05_smetana_crossfeeding/  SMETANA cross-feeding data + README (Fig 5c, Supp Fig 10)
│   ├── 06_saturation/            sat_sim/sat_sim_notebook.Rmd → complementarity vs community size (Fig 6b)
│   └── 07_community_composition/ assembly / binning / qfilter / drep Vis scripts (Supp Fig 3)
├── code/                     Python figure renders (plots.ipynb, image-analysis notebook)
├── figures/                  one folder per display item (fig1–6, supp_fig_01–11), each with a README
├── data/
│   ├── models/                   genome-scale metabolic models (SBML .xml, 102 genomes)
│   ├── genomes/                  assembled nucleotide + protein sequences, 16S
│   ├── annotations/              eggNOG, ISfinder, geNomad, MAG GFF (gzip-compressed)
│   ├── figure_data/              intermediate tables consumed directly by plotting code
│   ├── gem_curation/             GEM curation inputs (GPRs, media, memote)
│   ├── qc/                       assembly / binning statistics
│   └── smetana/                  SMETANA cross-feeding output
├── supplementary_tables/     Supplementary_Tables.xlsx (S1–S4) + 260-gene table
└── archive/                  superseded drafts and exploratory analyses (not in the paper)
```

## ▶️ Reproducing the analysis

Raw sequences are in the European Nucleotide Archive under [**PRJEB80563**](https://www.ebi.ac.uk/ena/browser/view/PRJEB80563). To reproduce the metagenomic pipeline from raw reads you need an HPC cluster and [metaGEM](https://github.com/franciscozorrilla/metaGEM); to reproduce the figures from the committed intermediate data you only need R and Python. The full software list (R and Python packages, external tool versions) is in [ENVIRONMENT.md](ENVIRONMENT.md); Python dependencies are in [requirements.txt](requirements.txt).

**Run everything from the repository root.** The three R notebooks (`analysis/03`, `04`, `06`) each begin by sourcing **`analysis/00_load_data.R`** — a shared loader that reads the committed data under `data/` and builds the objects they share (`checkm`, `auxo2`, `kemet_auxo_all`, `eggnog`, `isfinder_count`, `genes_count`, …) — so each notebook runs top-to-bottom from a defined set of inputs. Data reads resolve to `data/…`; regenerated tables and figures are written to `outputs/` (gitignored), while the committed published versions live in `figures/` and `data/`. Compressed inputs (`*.gz`) are read transparently. Stochastic steps are seeded (`set.seed(42)`) so the bootstrapped confidence intervals and complementarity curves are reproducible.

The committed `data/figure_data/figure3panelA–B.tsv` (Fig 3) and `figure4panelA–B.tsv` (Fig 4) are the intermediate tables the Python renders consume. A slow per-gene phyloglm loop in `phylogenetic.Rmd` is gated behind `recompute_genes_plr` (it defaults to loading the committed `genes_plr.tsv`). Pinned reproduction checks (from the committed data): Fig 3a Spearman ρ = **0.360**; Fig 4 medians **934 / 880** (genes·Mbp⁻¹) and **227 / 5.1** (IS·Mbp⁻¹); Fig 6b complementarity **29 / 65 / 87 / 95 %**; the differential-gene set = **260** genes.

## 🖼️ Main figures

| Panel | Shows | Code | Input data | Contributors | Reproduce / notes |
|---|---|---|---|---|---|
| **Fig 1a** | Nested sampling design (9 sites → soil column → 3 particles → 27 communities) | Schematic (BioRender) | — | **G. Yousif** | schematic |
| **Fig 1b** | Proportion of isolates auxotrophic per community | `code/plots.ipynb` | full auxotrophy screen (6,931 isolates; **G. Yousif**, not committed) | **G. Yousif** (data); **L. Oña** (render) | experimental; per-community proportion 10–50% (values in the render) |
| **Fig 1c** | Relative frequency of the number of amino-acid auxotrophies per strain | `code/plots.ipynb` | full auxotrophy screen (**G. Yousif**, not committed) | **G. Yousif** (data); **L. Oña** (render) | experimental; 73% require multiple amino acids |
| **Fig 1d** | Relative frequency of each individual amino-acid auxotrophy | `code/plots.ipynb` | full auxotrophy screen (**G. Yousif**, not committed) | **G. Yousif** (data); **L. Oña** (render) | experimental; Val 59%, Leu 47%, Ile 27% |
| **Fig 2a** | Phylogenetic distribution of auxotrophs (n=390) vs prototrophs (n=240), 16S | 16S rRNA circular phylogeny | 16S sequences | **G. Yousif** | Firmicutes / Actinobacteria / Proteobacteria |
| **Fig 2b** | Community richness vs auxotroph proportion | `code/plots.ipynb` | richness + auxotroph proportion (**G. Yousif**, not committed) | **G. Yousif** (data); **L. Oña** (render) | Spearman ρ = 0.60, P = 9.7×10⁻⁴, n = 27 |
| **Fig 3a** | GEM-predicted vs observed auxotrophies | `code/plots.ipynb` (cells 15/16) | `figure3panelA.tsv` ← `auxopred_default.tsv` + `auxolab_summary.tsv` | **F. Zorrilla** (analysis); **L. Oña** (render) | Spearman ρ = 0.360, N = 62 |
| **Fig 3b** | Annotation-based (KEMET) vs observed auxotrophies | `code/plots.ipynb` (cells 20/21) | `figure3panelB.tsv` ← `kemet_summary.tsv` + `auxolab_summary.tsv` | **F. Zorrilla** (analysis); **L. Oña** (render) | Spearman ρ = −0.071, N = 62 |
| **Fig 4a** | Annotated genes per Mbp (metagenomic / prototroph / auxotroph) | `metabolic_modeling_soil.Rmd`; `code/plots.ipynb` | `figure4panelA.tsv` | **F. Zorrilla** (analysis); **L. Oña** (render) | medians **934** (auxo, n=53) / **880** (proto, n=9); Mann-Whitney BH P = 0.013, W = 106 |
| **Fig 4b** | Insertion sequences per Mbp | `metabolic_modeling_soil.Rmd`; `code/plots.ipynb` | `figure4panelB.tsv` ← `ISfinder_k1000.tsv` | **F. Zorrilla** (analysis); **L. Oña** (render) | medians **227** (auxo, n=53) / **5.1** (proto, n=9); Mann-Whitney BH P = 0.016, W = 99 |
| **Fig 5a** | Pairwise cross-feeding networks (19 auxo + 22 proto) | experimental coculture analysis | coculture growth data | **L. Oña / G. Yousif** | experimental |
| **Fig 5b** | Proportion of auxotroph recipients that grew, by donor type | `code/Image processing and colony isolation_all species.ipynb` | coculture growth data | **L. Oña / G. Yousif** | experimental |
| **Fig 5c** | SMETANA-predicted metabolic interactions (focal community) | `analysis/05_smetana_crossfeeding/` | `data/smetana/` (SMETANA v1.1.0 `--detailed`, AA + vitamins) | **G. Yousif** (figure); **F. Zorrilla** (SMETANA data) | source data provided |
| **Fig 6a** | Complementarity vs community size across the 27 communities (per-community grey, mean in colour; 746 auxotrophic strains, sizes 2–5) | `code/plots.ipynb` (cell 1) | per-community complementarity summary (**L. Oña**, not committed) | **L. Oña** | contributed by L. Oña; see Data availability |
| **Fig 6b** | Complementarity of the 58 focal isolates: computational (red) vs experimental (blue) | `sat_sim/sat_sim_notebook.Rmd`; `code/plots.ipynb` (cell 14) | `analysis/06_saturation/sat_sim/*_auxo.tsv`, `data/figure_data/auxo_df.tsv` | **F. Zorrilla** (analysis); **L. Oña** (render) | computational **29/65/87/95 %**; experimental 44/57/75 % |

## 📑 Supplementary figures

| Supp Fig | Shows | Code | Input data | Contributors |
|---|---|---|---|---|
| **S1** | Locations of the 9 soil sampling sites (map) | sampling-site map | site coordinates | **G. Yousif** |
| **S2** | Community-level phylum profile of the 9 metagenomic samples | metagenome taxonomic composition | metagenome taxonomic profiles | **G. Yousif** |
| **S3** | De-novo assembly & binning statistics (102 genomes): genomes reconstructed, completeness vs contamination, basepairs/contigs | `analysis/07_community_composition/{assembly,binning,qfilter,drep}Vis.R` | assembly/binning/qfilter QC (CheckM, GTDB-Tk) | **F. Zorrilla** |
| **S4** | Phylogenetic tree of 102 genomes + auxotrophy heatmaps | `phylogenetic.Rmd` | protein tree; auxotrophy calls | **F. Zorrilla** |
| **S5** | Mobile genetic elements (plasmid / virus / insertion sequences per Mbp) | `phylogenetic.Rmd`; `metabolic_modeling_soil.Rmd` (geNomad) | `ISfinder_k1000.tsv`; geNomad summaries | **F. Zorrilla** |
| **S6** | Genome size (Mbp) and annotated gene number (102 genomes) | `metabolic_modeling_soil.Rmd` | `genome_contigs_bp.tsv`; `gtdbtk_summary.tsv` | **F. Zorrilla** |
| **S7** | PCA of gene copy number (eggNOG annotations) | `metabolic_modeling_soil.Rmd` | `eggnog_summary.tsv.gz` | **F. Zorrilla** |
| **S8** | Auxotrophy number vs tendency to grow in coculture | `code/plots.ipynb` | coculture growth data | **L. Oña** |
| **S9** | Auxotroph growth vs amount of amino acid produced (LC-MS/MS) | LC-MS/MS analysis | LC-MS/MS measurements | **G. Yousif** |
| **S10** | Predicted metabolite-exchange network (focal community) | `analysis/05_smetana_crossfeeding/` | `data/smetana/` | **G. Yousif** (figure); **F. Zorrilla** (SMETANA data) |
| **S11** | Distribution of unadjusted phyloglm P-values | `phylogenetic.Rmd` | `genes_plr.tsv` | **F. Zorrilla** |

## 📋 Supplementary tables

| Table | Contents | Source | Contributors |
|---|---|---|---|
| **S1** | Genomic evidence for AA biosynthesis in false-negative model predictions | `Supplementary_Tables.xlsx` (tab 1) | **F. Zorrilla** |
| **S2** | Summary of phyloglm models (λ, P, half-life + 95% CI, 1000 bootstraps, `set.seed(42)`) | `Supplementary_Tables.xlsx` (tab 2); `phylogenetic.Rmd` → `genes_plr.tsv` | **F. Zorrilla** |
| **S3** | Isolates used in the pairwise coculture assay (19 auxo / 22 proto) | `Supplementary_Tables.xlsx` (tab 3) | **G. Yousif** |
| **S4** | Isolates used in the complementarity experiments, Fig 6b | `Supplementary_Tables.xlsx` (tab 4) | **G. Yousif** |
| **260-gene table** | Differentially abundant gene set (230 depleted + 30 enriched) | `Supp_Table_260_differential_genes.xlsx`; `metabolic_modeling_soil.Rmd` → `diff_genes.tsv` | **F. Zorrilla** |

The 260-gene set is built on n = 51 genomes; the two additional genomes used in Fig. 4 are polymorphism clones and are excluded here.

## 🔬 Analysis pipeline (upstream of the figures)

| Stage | Code | Produces | Notes |
|---|---|---|---|
| Metagenomic pipeline | `analysis/01_metagem_pipeline/Snakefile` + `config.yaml` (metaGEM) | assembly, binning, GEMs | reconstructed GEMs are deposited in `data/models/`; rerunning the full pipeline needs an HPC cluster and the raw reads (ENA PRJEB80563) |
| GEM auxotrophy prediction | `analysis/02_gem_auxotrophy/get_auxo.py` | `auxopred_default.tsv` (Fig 3a) | default/canonical calls; sensitivity-analysis variants are in `archive/` |
| Phylogenetics + phyloglm | `analysis/03_phylogenetics/phylogenetic.Rmd` | tree (S4), MGE/IS (S5), gene associations, Supp Table 2, S11 | bootstraps use `set.seed(42)` |
| Genome features | `analysis/04_genome_features/metabolic_modeling_soil.Rmd` | Fig 4, S6, S7, geNomad (S5), 260-gene set, SMETANA plots | master analysis notebook |
| SMETANA cross-feeding | `analysis/05_smetana_crossfeeding/` | Fig 5c, S10 | site-based community models: `data/smetana/` (SMETANA v1.1.0 `--detailed`, AA + vitamins; see analysis/05 README) |
| Saturation | `analysis/06_saturation/sat_sim/sat_sim_notebook.Rmd` | Fig 6b | sampling uses `set.seed(42)` |
| Assembly / binning QC | `analysis/07_community_composition/{assembly,binning,qfilter,drep}Vis.R` | Supp Fig 3 | genomes reconstructed, completeness/contamination, basepairs/contigs |

## 📦 Data availability

- Raw sequencing reads are in the ENA under **[PRJEB80563](https://www.ebi.ac.uk/ena/browser/view/PRJEB80563)** (not stored here).
- The intermediate tables under `data/figure_data/` are sufficient to regenerate the computational figures (Fig. 3, 4, 6b; Supplementary Figs 3–7, 11) without the raw data. Large annotation and cross-feeding outputs are gzip-compressed in place (e.g. `data/annotations/eggnog_summary.tsv.gz`) and read transparently. The SMETANA detailed output is committed under `data/smetana/` (per-site runs in `site_based_smetana_d/` and the combined `global_detailed.tsv.gz`); the full raw output is available on request.
- Several figures report experimental measurements (coculture assays for Fig. 5a/b and Supplementary Fig. 8; LC-MS/MS for Supplementary Fig. 9; microscopy) or were contributed by co-authors (Fig. 1a–d, 2a/2b, 6a; Supplementary Figs 1, 2, 10). Their underlying measurements are provided as source data with the paper or are available from the corresponding authors.
- Data are released under CC-BY-4.0; code under MIT.

## 🗄️ Not in the paper → `archive/`

Superseded drafts and exploratory material are kept in `archive/` rather than presented as figures:

| Item | Reason |
|---|---|
| `archive/plots.Rmd` (+ `plots.nb.html`) | earlier plotting notebook, superseded by `metabolic_modeling_soil.Rmd` + `code/plots.ipynb` |
| `archive/newfig3.*`, `archive/updated_fig3.*` | superseded composite figure drafts |
| `archive/scatter_annotations_GEMs.*` | earlier GEM-vs-annotation scatter, superseded by the published Fig. 3 |
| `archive/enriched_COGs.pdf` | gene-enrichment plot supporting the Supplementary Notes; not a display item |
| `archive/get_auxo_*lenient.py`, `archive/auxopred_*.tsv` | lenient auxotrophy-call variants used for sensitivity analysis |
| `archive/compositionVis.R` (+ `assembled_vs_binned*.pdf`, `percent_mapping.pdf`) | exploratory community-composition / mapping QC script; its outputs are not published-figure panels |

## 📄 Licence

Code is released under the MIT Licence ([LICENSE](LICENSE)); data under CC-BY-4.0 ([DATA_LICENSE](DATA_LICENSE)).

## 🖋️ Citation

If you use this resource, please cite:

> Obligate cross-feeding of metabolites is common in soil microbial communities. Ghada Yousif*, Francisco Zorrilla*, Swagatika Dash, Leonardo Oña, Aditi Shekhar, Samir Giri, Rui Guan, Sharvari Harshe, Michael Itermann, Daphne Welter, Vladimir Benes, Kiran R. Patil, Christian Kost. bioRxiv 2025.01.29.635426; doi: https://doi.org/10.1101/2025.01.29.635426
