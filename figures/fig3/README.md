# Figure 3 — Metabolic modelling and annotation-based predictions of auxotrophies

| Panel | Shows | Contributors | Code | Inputs | Reproduce / notes |
|---|---|---|---|---|---|
| Fig. 3a | GEM-predicted vs experimentally observed amino-acid auxotrophies | F. Zorrilla (analysis); L. Oña (Python render) | code/plots.ipynb | data/figure_data/figure3panelA.tsv (from auxopred_default.tsv + auxolab_summary.tsv) | Spearman rho = 0.360, N = 62 |
| Fig. 3b | Annotation-based (KEMET) predicted vs observed amino-acid auxotrophies | F. Zorrilla (analysis); L. Oña (Python render) | code/plots.ipynb | data/figure_data/figure3panelB.tsv (from kemet_summary.tsv + auxolab_summary.tsv) | Spearman rho = -0.071, N = 62 |

*Run the analysis scripts from the repository root; the figure regenerates into `outputs/`. This folder documents the figure-to-code mapping; the final published image is not committed here — see the Contributors and Code columns for how each panel was produced (R render, Python render by L. Oña, or an experimental / co-author figure).*
*See `../../README.md` and `../../MANIFEST.tsv` for the full map.*
