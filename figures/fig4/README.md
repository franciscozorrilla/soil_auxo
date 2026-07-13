# Figure 4 — Coding density and insertion sequences in auxotrophic genomes

| Panel | Shows | Contributors | Code | Inputs | Reproduce / notes |
|---|---|---|---|---|---|
| Fig. 4a | Annotated genes per Mbp for metagenomic, prototrophic and auxotrophic genomes | F. Zorrilla (analysis); L. Oña (Python render) | analysis/04_genome_features/metabolic_modeling_soil.Rmd; code/plots.ipynb | data/figure_data/figure3panelC.tsv | Medians 934 (auxotroph, n=53) / 880 (prototroph, n=9); Mann-Whitney BH P = 0.013, W = 106 |
| Fig. 4b | Insertion sequences per Mbp | F. Zorrilla (analysis); L. Oña (Python render) | analysis/04_genome_features/metabolic_modeling_soil.Rmd; code/plots.ipynb | data/figure_data/figure3panelD.tsv (from ISfinder_k1000.tsv) | Medians 227 (auxotroph, n=53) / 5.1 (prototroph, n=9); Mann-Whitney BH P = 0.016, W = 99 |

*Run the analysis scripts from the repository root; the figure regenerates into `outputs/`. This folder documents the figure-to-code mapping; the final published image is not committed here — see the Contributors and Code columns for how each panel was produced (R render, Python render by L. Oña, or an experimental / co-author figure).*
*See `../../REPO_MAP.md` and `../../MANIFEST.tsv` for the full map.*
