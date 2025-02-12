# 🧫 soil_auxo

Repository with code & data for the publication *Obligate cross-feeding of metabolites is common in soil microbial communities*.

### System requirements

To reproduce the analysis starting from raw sequencing data, you will need access to a high performance computer cluster. To reproduce plots starting from intermediate files, a standard laptop or desktop is required.

### Installation guide

A detailed [installation guide for metaGEM](https://github.com/franciscozorrilla/metaGEM/tree/master/config) can be found in the GitHub repo. To set up Rstudio for plotting, follow [this guide](https://rstudio-education.github.io/hopr/starting.html).

### Demo

Examples of how to run metaGEM are available as a [tutorial](https://github.com/franciscozorrilla/unseenbio_metaGEM) and [interactive demo with google colab](https://colab.research.google.com/drive/1I1S8AoGuJ9Oc2292vqAGTDmZcbnolbuj#scrollTo=awiAaVwSF5Fz).

### Instructions for use

Please see the [file](https://github.com/franciscozorrilla/soil_auxo/blob/main/scripts/plots.Rmd) under `scripts/plots.Rmd` for detailed and commentented code on how to reproduce plots.

## 🧬 Abstract

Many microorganisms are refractory to laboratory cultivation. One possible explanation, known as the great plate count anomaly, is metabolic dependencies among community members. However, systematic studies of these interactions in communities like soil are missing, hindering advances in understanding these ecologically important ecosystems. Here, we address this issue by systematically analysing 6,931 bacterial isolates of 27 soil microbial communities. We find that the growth of up to 50% of all community members depended essentially on supplementation with amino acids, vitamins, or nucleotides. In 73% of cases, supplementation with multiple amino acids was necessary. Genomic analysis of 62 strains revealed that accumulation of insertion sequences and specific gene loss was associated with the observed auxotrophies. Finally, genome-scale metabolic models, computational analyses, and cocultivation experiments demonstrated that other co-occurring genotypes complemented the metabolic needs of auxotrophs, thus facilitating their growth. Our results demonstrate that soil bacteria exist within integrated metabolic networks, which hampers their cultivation.

## 🖋️ Citation

If you use this resource, please cite:

 > Obligate cross-feeding of metabolites is common in soil microbial communities. Ghada Yousif*, Francisco Zorrilla*, Swagatika Dash, Leonardo Oña, Aditi Shekhar, Samir Giri, Rui Guan, Sharvari Harshe, Michael Itermann, Daphne Welter, Vladimir Benes, Kiran R. Patil, Christian Kost
bioRxiv 2025.01.29.635426; doi: https://doi.org/10.1101/2025.01.29.635426 

## 🧪 Usage & description

This repo contains supplementary files, sequences, and metabolic models associated with isolates and metagenomic samples from our study. Raw sequences have been desposited in the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/home) (ENA) under accession code [PRJEB80563](https://www.ebi.ac.uk/ena/browser/view/PRJEB80563). 

Clone repo:

```
$ git clone https://github.com/franciscozorrilla/soil_auxo.git
```

See repo structure:

```
$ tree soil_auxo/ -L 1
.
├── 16S
├── LICENSE
├── README.md
├── binning_stats
├── dna
├── final_figs_data
├── humann3
├── media
├── models
├── protein
└── supplementary_files

10 directories, 2 files
```

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13847512.svg)](https://doi.org/10.5281/zenodo.13847512)
