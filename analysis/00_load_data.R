# =============================================================================
# Shared data loader for the soil_auxo analysis notebooks.
#
# This script reads the committed intermediate data under data/ and builds the
# objects that are shared across analysis/03_phylogenetics/phylogenetic.Rmd,
# analysis/04_genome_features/metabolic_modeling_soil.Rmd, and
# analysis/06_saturation/sat_sim/sat_sim_notebook.Rmd. Each of those notebooks
# calls source("analysis/00_load_data.R") in its setup chunk, so the objects are
# available before any analysis chunk runs. This removes the cross-notebook
# object dependencies and forward references, so each notebook runs top-to-bottom.
#
# Run from the repository root (paths are relative to it). The code here is
# copied verbatim from the notebooks (cited as file:line) so it produces exactly
# the same objects; the notebooks may re-create some of these objects later in
# their own chunks, which is harmless.
# =============================================================================

library(dplyr)
library(tidyr)

# ---- Tier 0: base committed reads ----
# metabolic_modeling_soil.Rmd:35,37,39
bigg_mets = read.delim("data/figure_data/bigg_classes.tsv")
metadata  = read.delim("data/figure_data/soil_auxo_metadata_full.tsv")
ions = c("zn2","so4","k","pi","o2","mn2","mg2","ca2","cl","h","cobalt2","cu2","fe2")
# metabolic_modeling_soil.Rmd:277
contigs_bp = read.delim("data/figure_data/genome_contigs_bp.tsv")

# lab-observed auxotrophies (metabolic_modeling_soil.Rmd:42-69)
auxolab = read.delim("data/figure_data/auxolab_summary.tsv") %>%
  replace(is.na(.), 0) %>%
  select(c(Sample_ID,AAs_auxo:Nucleotide_auxo,Ala:Val)) %>%
  pivot_longer(!Sample_ID,values_to = "auxo_lab", names_to = "metabolite") %>%
  mutate(name=gsub("Ala","L-Alanine",metabolite)) %>%
  mutate(name=gsub("Arg","L-Arginine",name)) %>%
  mutate(name=gsub("Asp","L-Aspartic acid",name)) %>%
  mutate(name=gsub("Asn","L-Asparagine",name)) %>%
  mutate(name=gsub("Cys","L-Cysteine",name)) %>%
  mutate(name=gsub("Glu","L-Glutamate",name)) %>%
  mutate(name=gsub("Gln","L-Glutamine",name)) %>%
  mutate(name=gsub("Gly","Glycine",name)) %>%
  mutate(name=gsub("His","L-Histidine",name)) %>%
  mutate(name=gsub("Ile","L-Isoleucine",name)) %>%
  mutate(name=gsub("Leu","L-Leucine",name)) %>%
  mutate(name=gsub("Lys","L-Lysine",name)) %>%
  mutate(name=gsub("Met","L-Methionine",name)) %>%
  mutate(name=gsub("Phe","L-Phenylalanine",name)) %>%
  mutate(name=gsub("Pro","L-Proline",name)) %>%
  mutate(name=gsub("Ser","L-Serine",name)) %>%
  mutate(name=gsub("Thr","L-Threonine",name)) %>%
  mutate(name=gsub("Trp","L-Tryptophan",name)) %>%
  mutate(name=gsub("Tyr","L-Tyrosine",name)) %>%
  mutate(name=gsub("Val","L-Valine",name)) %>%
  mutate(name=gsub("AAs_auxo","*aminoacids",name)) %>%
  mutate(name=gsub("Vitamin_auxo","*vitamins",name)) %>%
  mutate(name=gsub("Nucleotide_auxo","*nucleotides",name)) %>%
  mutate(auxo_lab=as.factor(auxo_lab))

# KEGG modules for the 14 amino acids with dedicated biosynthesis pathways
# (metabolic_modeling_soil.Rmd:123-125)
AA_KEGG_modules_ID = c("M00020","M00018","M00021","M00338","M00609","M00017","M00019","M00535","M00570","M00432","M00016","M00525","M00526","M00527","M00030","M00433","M00031","M00844","M00845","M00015","M00026","M00022","M00023","M00024","M00910","M00025","M00040")
AA_KEGG_modules_meta = c("ser/thr","ser/thr","cys/met","cys/met","cys/met","cys/met","BCAA","BCAA","BCAA","BCAA","lysine","lysine","lysine","lysine","lysine","lysine","lysine","arg/pro","arg/pro","arg/pro","his","aromatic","aromatic","aromatic","aromatic","aromatic","aromatic")
AA_KEGG_modules = data.frame(AA_KEGG_modules_ID,AA_KEGG_modules_meta) %>% rename(KEGG_Module=AA_KEGG_modules_ID)

# ---- Tier 1: model predictions + genome metadata ----
# GEM-predicted auxotrophies (default/canonical calls) (metabolic_modeling_soil.Rmd:74-82)
auxopred2 = read.delim("data/figure_data/auxopred_default.tsv") %>%
    mutate(present=1) %>% filter(metabolite!="") %>%
    pivot_wider(names_from = "metabolite",values_from = "present",values_fill = 0) %>%
    mutate(ala__L=0,arg__L=0,asp__L=0,cys__L=0,glu__L=0,gly=0,leu__L=0,met__L=0,pro__L=0,ser__L=0,tyr__L=0) %>%
    pivot_longer(!model,values_to = "auxo_model",names_to = "compound") %>%
    rename(genome=model) %>%
    left_join(.,bigg_mets %>% select(name,compound),by="compound") %>%
    mutate(Sample_ID = gsub("_.*$","",genome), auxo_model = as.factor(auxo_model)) %>%
    mutate(name=gsub("L-Glutamic acid","L-Glutamate",name)) %>% filter(!compound %in% ions)

# enriched genome table (checkm + metadata + contig/size) (metabolic_modeling_soil.Rmd:280-284)
checkm = read.delim("data/figure_data/checkm_summary.tsv") %>%
  dplyr::rename(genome=Bin.Id) %>%
  mutate(Quality=ifelse(Completeness>=90&Contamination<=5,"HQ","MQ")) %>%
  mutate(Sample_ID=gsub("_.*$","",genome)) %>%
  left_join(.,metadata,by="Sample_ID") %>% left_join(.,contigs_bp,by="genome")

# ---- Tier 2: joins onto the above ----
# GEM predictions vs lab observations (metabolic_modeling_soil.Rmd:85-87)
auxo2 = left_join(auxolab,auxopred2 %>% filter(!compound%in%ions) ,by=c("name","Sample_ID")) %>%
    drop_na() %>%
    mutate(agree=as.numeric(ifelse(auxo_lab==auxo_model,"1","0")))

# taxonomy (metabolic_modeling_soil.Rmd:332-335)
gtdbtk = read.delim("data/figure_data/gtdbtk_summary.tsv") %>%
  rename(genome=user_genome) %>%
  separate(.,classification,sep = ";",c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  left_join(.,checkm,by="genome")

# KEMET pathway-completeness table (metabolic_modeling_soil.Rmd:127-131)
kemet = read.delim("data/figure_data/kemet_summary.tsv",header = FALSE) %>% select(-V8) %>%
    rename(genome=V1,KEGG_Module=V2,module_name=V3,completeness=V4,nPresent__nTotal=V5,missing_KEGG_ko=V6,present_KEGG_ko=V7) %>%
    separate(.,nPresent__nTotal,into=c("n_present","n_total"),sep="__") %>%
    mutate(n_present=as.numeric(n_present),n_total=as.numeric(n_total),n_missing=n_total-n_present) %>%
    left_join(.,checkm,by="genome")

# insertion-sequence hits (phylogenetic.Rmd:438-440)
isfinder_hits = read.delim("data/annotations/ISfinder_k1000.tsv") %>%
  mutate(genome=gsub("_shovill.*$","_shovill",reference),genome=gsub("_NODE.*$","",genome),genome=gsub("_k.*$","",genome)) %>%
  left_join(.,checkm,by="genome")

# ---- Tier 3: annotation-derived (eggNOG is large; this step is the slow one) ----
# eggNOG annotations with COG-category names (metabolic_modeling_soil.Rmd:342-368)
eggnog = read.delim("data/annotations/eggnog_summary.tsv.gz") %>%
  mutate(cog_name=gsub("^A$","rna processing and modification",COG_category)) %>%
  mutate(cog_name=gsub("^B$","chromatin structure and dynamics",cog_name)) %>%
  mutate(cog_name=gsub("^C$","energy production and conversion",cog_name)) %>%
  mutate(cog_name=gsub("^D$","cell cycle control and mitosis",cog_name)) %>%
  mutate(cog_name=gsub("^E$","amino acid metabolism and transport",cog_name)) %>%
  mutate(cog_name=gsub("^F$","nucleotide metabolism and transport",cog_name)) %>%
  mutate(cog_name=gsub("^G$","carbohydrate metabolism and transport",cog_name)) %>%
  mutate(cog_name=gsub("^H$","coenzyme metabolism",cog_name)) %>%
  mutate(cog_name=gsub("^I$","lipid metabolism",cog_name)) %>%
  mutate(cog_name=gsub("^J$","tranlsation",cog_name)) %>%
  mutate(cog_name=gsub("^K$","transcription",cog_name)) %>%
  mutate(cog_name=gsub("^L$","replication and repair",cog_name)) %>%
  mutate(cog_name=gsub("^M$","cell wall/membrane/envelop biogenesis",cog_name)) %>%
  mutate(cog_name=gsub("^N$","cell motility",cog_name)) %>%
  mutate(cog_name=gsub("^O$","post-translational modification, protein turnover, chaperone functions",cog_name)) %>%
  mutate(cog_name=gsub("^P$","inorganic ion transport and metabolism",cog_name)) %>%
  mutate(cog_name=gsub("^Q$","secondary structure",cog_name)) %>%
  mutate(cog_name=gsub("^T$","signal transduction",cog_name)) %>%
  mutate(cog_name=gsub("^U$","intracellular trafficing and secretion",cog_name)) %>%
  mutate(cog_name=gsub("^V$","defense mechanism",cog_name)) %>%
  mutate(cog_name=gsub("^Y$","nuclear structure",cog_name)) %>%
  mutate(cog_name=gsub("^Z$","cytoskeleton",cog_name)) %>%
  mutate(cog_name=gsub("^R$","general functional prediction only",cog_name)) %>%
  mutate(cog_name=gsub("^S$","function unknown ",cog_name)) %>%
  mutate(cog_name=gsub("-","function unknown ",cog_name)) %>%
  left_join(.,gtdbtk,by="genome")

# ---- Tier 4+: objects derived from the annotation reads ----
# annotation/KEMET-based auxotrophy calls (metabolic_modeling_soil.Rmd:155-200; built then the `name` column is added)
kemet %>% filter(KEGG_Module %in% AA_KEGG_modules$KEGG_Module) %>%
    mutate(metabolite=ifelse(KEGG_Module=="M00017","Met",
                      ifelse(KEGG_Module=="M00015","Pro",
                      ifelse(grepl("M00016|M00030|M00031|M00525|M00526|M00527|M00433",KEGG_Module),"Lys",
                      ifelse(KEGG_Module=="M00018","Thr",
                      ifelse(KEGG_Module=="M00019","Ile,Val",
                      ifelse(KEGG_Module=="M00020","Ser",
                      ifelse(grepl("M00021|M00338|M00609",KEGG_Module),"Cys",
                      ifelse(KEGG_Module=="M00023","Trp",
                      ifelse(grepl("M00024|M00910",KEGG_Module),"Phe",
                      ifelse(grepl("M00025|M00040",KEGG_Module),"Tyr",
                      ifelse(KEGG_Module=="M00026","His",
                      ifelse(KEGG_Module=="M00570|M00535","Ile",
                      ifelse(grepl("M00844|M00845",KEGG_Module),"Arg",
                      ifelse(grepl("M00432",KEGG_Module),"Leu","other"))))))))))))))) %>%
    relocate(metabolite) %>% filter(metabolite!="other") %>%
    separate_rows(metabolite,sep=",") %>% mutate(comp_num=ifelse(completeness=="COMPLETE",1,0)) %>%
    group_by(genome,metabolite) %>% summarize(proto=sum(comp_num)) %>% mutate(auxo_kemet=ifelse(proto>0,0,1)) %>%
    mutate(Sample_ID = gsub("_.*$","",genome), auxo_kemet=as.factor(auxo_kemet)) %>%
  left_join(.,auxo2,by=c("genome","metabolite")) -> kemet_auxo_all
kemet_auxo_all %>%
  mutate(name=gsub("Ala","L-Alanine",metabolite)) %>%
  mutate(name=gsub("Arg","L-Arginine",name)) %>%
  mutate(name=gsub("Asp","L-Aspartic acid",name)) %>%
  mutate(name=gsub("Asn","L-Asparagine",name)) %>%
  mutate(name=gsub("Cys","L-Cysteine",name)) %>%
  mutate(name=gsub("Glu","L-Glutamate",name)) %>%
  mutate(name=gsub("Gln","L-Glutamine",name)) %>%
  mutate(name=gsub("Gly","Glycine",name)) %>%
  mutate(name=gsub("His","L-Histidine",name)) %>%
  mutate(name=gsub("Ile","L-Isoleucine",name)) %>%
  mutate(name=gsub("Leu","L-Leucine",name)) %>%
  mutate(name=gsub("Lys","L-Lysine",name)) %>%
  mutate(name=gsub("Met","L-Methionine",name)) %>%
  mutate(name=gsub("Phe","L-Phenylalanine",name)) %>%
  mutate(name=gsub("Pro","L-Proline",name)) %>%
  mutate(name=gsub("Ser","L-Serine",name)) %>%
  mutate(name=gsub("Thr","L-Threonine",name)) %>%
  mutate(name=gsub("Trp","L-Tryptophan",name)) %>%
  mutate(name=gsub("Val","L-Valine",name)) %>%
  mutate(name=gsub("AAs_auxo","*aminoacids",name)) %>%
  mutate(name=gsub("Vitamin_auxo","*vitamins",name)) %>%
  mutate(name=gsub("Nucleotide_auxo","*nucleotides",name)) -> kemet_auxo_all

# normalized insertion-sequence counts per genome (phylogenetic.Rmd:443)
isfinder_hits %>% group_by(genome)%>% filter(seqid>=40) %>% dplyr::summarize(count=n()) %>% full_join(checkm,by="genome") %>% mutate_at("count", ~replace_na(.,0)) %>% mutate(count_norm=count/size) -> isfinder_count

# per-IS-family counts across genomes (metabolic_modeling_soil.Rmd:580)
isfinder_hits %>% group_by(genome,query) %>% summarise(count=n()) %>% pivot_wider(names_from = "query",values_from = "count",values_fill = 0) %>% pivot_longer(!genome,names_to = "query",values_to = "count") %>% left_join(.,checkm) %>% drop_na() %>% ungroup() -> is_sum

# total annotated genes per genome, normalized per Mbp (metabolic_modeling_soil.Rmd:378)
eggnog %>% group_by(genome) %>% summarize(count=n()) %>% left_join(checkm) %>% mutate(count_norm=1000000*count/size) -> genes_count

# geNomad plasmid/virus predictions (used by the Supp Fig 5 MGE analysis in both notebooks)
# (metabolic_modeling_soil.Rmd:1106-1108, 1167-1169)
genomad_plasmid_genes = read.delim("data/annotations/genomad/cat_mags_plasmid_genes.tsv")
genomad_plasmid = read.delim("data/annotations/genomad/cat_mags_plasmid_summary.tsv") %>% mutate(Sample_ID=gsub("_.*$","",seq_name),genome=gsub("^([^_]+_[^_]+).*", "\\1",seq_name)) %>% left_join(.,metadata)
genomad_virus_genes = read.delim("data/annotations/genomad/cat_mags_virus_genes.tsv")
genomad_virus = read.delim("data/annotations/genomad/cat_mags_virus_summary.tsv") %>% mutate(Sample_ID=gsub("_.*$","",seq_name),genome=gsub("^([^_]+_[^_]+).*", "\\1",seq_name)) %>% left_join(.,metadata)

# differentially abundant gene set (phylogenetic.Rmd:288) -- committed result
diff_genes = read.delim("data/figure_data/diff_genes.tsv")
