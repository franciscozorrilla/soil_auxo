# Figure 1 — Auxotrophs are prevalent in soil microbial communities

| Panel | Shows | Contributors | Code | Inputs | Reproduce / notes |
|---|---|---|---|---|---|
| Fig. 1a | Nested sampling design (nine sites; soil column; three particles per column; 27 communities) | G. Yousif | Schematic (BioRender) | — | Schematic illustration |
| Fig. 1b | Proportion of isolates that are auxotrophic in each community | G. Yousif (data); L. Oña (Python render) | code/plots.ipynb | data/figure_data/auxolab_summary.tsv | Per-community auxotroph proportion (range 10-50%) |
| Fig. 1c | Relative frequency of the number of amino-acid auxotrophies per strain | G. Yousif (data); L. Oña (Python render) | code/plots.ipynb | data/figure_data/auxolab_summary.tsv | Distribution of auxotrophies per strain (73% require multiple amino acids) |
| Fig. 1d | Relative frequency of each individual amino-acid auxotrophy | G. Yousif (data); L. Oña (Python render) | code/plots.ipynb | data/figure_data/auxolab_summary.tsv | Per-amino-acid auxotrophy frequency (Val 59%, Leu 47%, Ile 27%) |

*Run the analysis scripts from the repository root; the figure regenerates into `outputs/`. This folder documents the figure-to-code mapping; the final published image is not committed here — see the Contributors and Code columns for how each panel was produced (R render, Python render by L. Oña, or an experimental / co-author figure).*
*See `../../REPO_MAP.md` and `../../MANIFEST.tsv` for the full map.*
