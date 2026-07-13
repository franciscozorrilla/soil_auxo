# Figure 2 — Species richness correlates with the number of auxotrophic strains

| Panel | Shows | Contributors | Code | Inputs | Reproduce / notes |
|---|---|---|---|---|---|
| Fig. 2a | Phylogenetic distribution of auxotrophic (n=390) and prototrophic (n=240) isolates from 16S rRNA | G. Yousif | 16S rRNA circular phylogeny | data/genomes/16S sequences | Phylum-level distribution (Firmicutes, Actinobacteria, Proteobacteria) |
| Fig. 2b | Community species richness vs the proportion of auxotrophic strains | G. Yousif (data); L. Oña (Python render) | code/plots.ipynb | richness and auxotroph proportion per community | Spearman rho = 0.60, P = 9.7e-4, n = 27 |

*Run the analysis scripts from the repository root; the figure regenerates into `outputs/`. This folder documents the figure-to-code mapping; the final published image is not committed here — see the Contributors and Code columns for how each panel was produced (R render, Python render by L. Oña, or an experimental / co-author figure).*
*See `../../REPO_MAP.md` and `../../MANIFEST.tsv` for the full map.*
