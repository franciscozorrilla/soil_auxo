# SMETANA cross-feeding analysis (Fig 5c, Supplementary Fig 10)

Community metabolic interactions were predicted with **SMETANA v1.1.0** (`--detailed`),
combining the genome-scale models of isolates and metagenome-assembled genomes extracted
from the **same physical site** (the focal community, site 7 / Meadow_2), focusing on the
exchange of amino acids and vitamins.

- **Output data (committed):** `../../data/smetana/site_based_smetana_d/` (per-community
  `*_detailed.tsv`) and `../../data/smetana/global_detailed.tsv.gz` (full run, gzip-compressed).
- **Plotting:** the alluvial visualisation is a section of
  `../04_genome_features/metabolic_modeling_soil.Rmd`. The published Fig 5c / Supp Fig 10
  network figures were assembled by G. Yousif from this SMETANA output.

SMETANA is part of the [metaGEM](https://github.com/franciscozorrilla/metaGEM) pipeline
(see `../01_metagem_pipeline/`).
