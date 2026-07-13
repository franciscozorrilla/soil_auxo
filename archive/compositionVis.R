library(tidyverse)
library(tidytext)
library(ggpubr)

taxonomy=read.delim("data/qc/GTDBTk.stats",header=TRUE) %>% 
  select(user_genome,classification) %>% 
  separate(.,classification,into = c("kingdom","phylum","class","order","family","genus","species"),sep = ";")

abundance=read.delim("data/qc/abundance.stats",header=FALSE) %>% unique()
colnames(abundance)=c("user_genome","absolute_ab","rel_ab")
abundance %>% mutate(user_genome=gsub("\\.bin","_bin",user_genome)) -> abundance

metadata = read.delim("data/figure_data/soil_auxo_metadata_full.tsv")
taxonomy %>% mutate(Sample_ID = gsub("_.*$","",user_genome)) %>% left_join(.,metadata) -> taxonomy

taxab = left_join(taxonomy,abundance,by="user_genome")
taxab$sample = gsub("\\..*$","",taxab$user_genome)
taxab$species = gsub("s__$","Undefined sp.",taxab$species)
taxab$species = gsub("s__","",taxab$species)
taxab$genus = gsub("g__$","Undefined gen.",taxab$genus)
taxab$genus = gsub("g__","",taxab$genus)

ggplot(taxab%>% filter(species!="Undefined sp.")) +
  geom_bar(aes(x=reorder_within(Sample_ID,-rel_ab,Sample_ID),y=rel_ab*100,fill=species),stat="identity") + 
  scale_x_reordered() +
  facet_wrap(~sample,scales = "free") + 
  ylab("Relative abundance (%)") + 
  xlab("Species") +
  coord_flip()

ggsave("outputs/compositionVis.pdf",width = 30,height=30)

# metagenomic samples kingdom view
ggplot(taxab%>% filter(Type=="Metagenomic")) +
  geom_bar(aes(x=reorder_within(Sample_ID,-rel_ab,Sample_ID),y=rel_ab*100,fill=kingdom),stat="identity") + 
  scale_x_reordered() +
  facet_wrap(~Habitats,scales = "free") + 
  ylab("Relative abundance (%)") + 
  xlab("Sample") + theme(axis.text.x = element_text(angle = 45, hjust=1))

# metagenomic samples phyun view
ggplot(taxab%>% filter(Type=="Metagenomic")) +
  geom_bar(aes(x=reorder_within(Sample_ID,-rel_ab,Sample_ID),y=rel_ab*100,fill=phylum),stat="identity") + 
  scale_x_reordered() +
  facet_wrap(~Habitats,scales = "free") + 
  ylab("Relative abundance (%)") + 
  xlab("Sample")

# metagenomic samples class view
ggplot(taxab%>% filter(Type=="Metagenomic")) +
  geom_bar(aes(x=reorder_within(Sample_ID,-rel_ab,Sample_ID),y=rel_ab*100,fill=class),stat="identity") + 
  scale_x_reordered() +
  facet_wrap(~Habitats,scales = "free") + 
  ylab("Relative abundance (%)") + 
  xlab("Sample") 

# metagenomic samples order view
ggplot(taxab%>% filter(Type=="Metagenomic")) +
  geom_bar(aes(x=reorder_within(Sample_ID,-rel_ab,Sample_ID),y=rel_ab*100,fill=order),stat="identity") + 
  scale_x_reordered() +
  facet_wrap(~Habitats,scales = "free") + 
  ylab("Relative abundance (%)") + 
  xlab("Sample") 

# metagenomic samples family view
ggplot(taxab%>% filter(Type=="Metagenomic")) +
  geom_bar(aes(x=reorder_within(Sample_ID,-rel_ab,Sample_ID),y=rel_ab*100,fill=family),stat="identity") + 
  scale_x_reordered() +
  facet_wrap(~Habitats,scales = "free") + 
  ylab("Relative abundance (%)") + 
  xlab("Sample") 

# metagenomic samples genus view
ggplot(taxab%>% filter(Type=="Metagenomic")) +
  geom_bar(aes(x=reorder_within(Sample_ID,-rel_ab,Sample_ID),y=rel_ab*100,fill=genus),stat="identity") + 
  scale_x_reordered() +
  facet_wrap(~Habitats,scales = "free") + 
  ylab("Relative abundance (%)") + 
  xlab("Sample") 

# metagenomic samples species view
ggplot(taxab%>% filter(Type=="Metagenomic")) +
  geom_bar(aes(x=reorder_within(Sample_ID,-rel_ab,Sample_ID),y=rel_ab*100,fill=species),stat="identity") + 
  scale_x_reordered() +
  facet_wrap(~Habitats,scales = "free") + 
  ylab("Relative abundance (%)") + 
  xlab("Sample") 



taxonomy %>% select(user_genome:species) %>% 
  pivot_longer(!user_genome,values_to = "taxa",names_to = "level") %>% 
  select(-user_genome) %>% group_by(taxa,level) %>% summarise(count=n()) %>% 
  ggplot() + geom_bar(aes(x=reorder(taxa,count),y=count),stat="identity") + facet_wrap(~level,scales = "free",ncol = 2) + coord_flip()

ggsave("outputs/taxonomy.pdf",height = 12,width = 8)

smag_map=read.delim("data/qc/percent_mapped.stats") %>% left_join(.,metadata)
assembly_map = read.delim("data/qc/percent_assembly.tsv") %>% left_join(smag_map) %>% pivot_longer(.,cols = c(percent_mapped,percent_assembly),names_to = "mapping",values_to = "value")

ggplot(assembly_map) + geom_jitter(aes(x=mapping,y=value)) +
  geom_boxplot(aes(x=mapping,y=value,fill=Type),alpha=0.4,outlier.shape = NA) +
  ylab("% of reads mapped to assemblies/MAGs/SAGs of samples") +
  theme(legend.title = element_blank())+ 
  theme(legend.position = "none")+
  facet_wrap(~Type) +
  coord_flip()

ggsave("outputs/percent_mapping.pdf")

#plot this, not such a good viz, try with % diffs
reassembled %>% 
  group_by(Sample_ID) %>% 
  summarize(mags_contigs = sum(contigs),mags_size = sum(size)) %>% 
  unique() %>% left_join(.,assembly) %>% rename(assembly_contigs = contigs, assembly_size = length_total) %>%
  pivot_longer(.,cols=c(mags_contigs,mags_size,assembly_contigs,assembly_size),names_to = "item",values_to = "value") %>% 
  unique() %>% ggplot() + geom_jitter(aes(x=item,y=value)) +
  geom_boxplot(aes(x=item,y=value,fill=Type),alpha=0.4,outlier.shape = NA) +
  ylab("% of reads mapped to assemblies/MAGs/SAGs of samples") +
  theme(legend.title = element_blank())+ 
  theme(legend.position = "none")+
  facet_grid(item~Type,scales="free") +
  coord_flip()

# very interesting, try to see if isolates that show chaos are also the lowest depth sequenced (nope, try GC content instead)
reassembled %>% 
  group_by(Sample_ID) %>% 
  summarize(mags_contigs = sum(contigs),mags_size = sum(size)) %>% unique() %>% 
  left_join(.,assembly) %>% rename(assembly_contigs = contigs, assembly_size = length_total) %>% 
  mutate(percent_bp=100*mags_size/assembly_size,percent_contigs = 100*mags_contigs/assembly_contigs) %>%
  unique() %>% ggplot() + geom_point(aes(percent_bp,percent_contigs,color=Type)) + xlab("% assembled bp retained in MAGs") + 
  ylab("% assembled contigs retained in MAGs") + facet_wrap(~Type,scales="free")

ggsave("outputs/assembled_vs_binned_percent.pdf")

reassembled %>% 
  group_by(Sample_ID) %>% 
  summarize(mags_contigs = sum(contigs),mags_size = sum(size)) %>% unique() %>% 
  left_join(.,assembly) %>% rename(assembly_contigs = contigs, assembly_size = length_total) %>% 
  mutate(percent_bp=100*mags_size/assembly_size,percent_contigs = 100*mags_contigs/assembly_contigs) %>%
  left_join(.,smag_map) %>%
  unique() %>% ggplot() + geom_point(aes(mags_contigs,assembly_contigs,color=Type)) + facet_wrap(~Type,scales="free",ncol=1) + theme(legend.position = "none") -> panel1
reassembled %>% 
  group_by(Sample_ID) %>% 
  summarize(mags_contigs = sum(contigs),mags_size = sum(size)) %>% unique() %>% 
  left_join(.,assembly) %>% rename(assembly_contigs = contigs, assembly_size = length_total) %>% 
  mutate(percent_bp=100*mags_size/assembly_size,percent_contigs = 100*mags_contigs/assembly_contigs) %>%
  left_join(.,smag_map) %>%
  unique() %>% ggplot() + geom_point(aes(mags_size,assembly_size,color=Type)) + facet_wrap(~Type,scales="free",ncol=1)-> panel2

ggarrange(panel1,panel2,widths = c(1,1.25))

ggsave("outputs/assembled_vs_binned.pdf",height = 6,width = 10)



# smetana
soilsmet=read.delim("detailed_summary.tsv",header = TRUE,sep="\t",row.names = NULL)

# manipulate data to summarize across simulations by calculating mean, sd, and median
soilsmet %>% 
  group_by(community,receiver,donor,compound) %>% 
  summarize(smet_ave=mean(smetana),smet_sd=sd(smetana),smet_med=median(smetana)) -> soilsmet_df

taxab$tax = paste(taxab$genus, taxab$species, sep=" ")
taxab%>% mutate(tax=ifelse(species=="Undefined sp.",tax,species)) ->taxab # bad coding

ggplot(soilsmet_df %>% 
         mutate(compound = gsub("M_","",compound),compound = gsub("_e","",compound)) %>% filter(., !grepl('^cu2$|^pi$|^fe2$|^fe3$|^h$|^s$|^h2s$', compound))%>% 
         left_join(.,taxab%>% rename(receiver=user_genome,receiver_tax=tax) %>% select(receiver,receiver_tax)) %>% 
         left_join(.,taxab%>% rename(donor=user_genome,donor_tax=tax) %>% select(donor,donor_tax,Type,Habitats)) %>% 
         left_join(.,bigg_mets) %>%
         filter(smet_ave>= 0.75) %>% unique(),
       aes(axis1 = donor_tax, axis2 = name, axis3 = receiver_tax,
           y = smet_ave)) +
  scale_x_discrete(limits = c("Donor", "Metabolite", "Reciever")) +
  xlab("Interaction") +
  geom_alluvium(aes(fill = name)) +
  geom_stratum(width=0.3) +
  theme_minimal() + geom_text(stat = "stratum", aes(label = after_stat(stratum)),min.y=0.2)+theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none") + 
  facet_wrap(community~Type~Habitats,scales="free")

ggsave("outputs/smetana_detailed_molweight.pdf",height = 10,width = 18)

