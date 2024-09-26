# libs
library(tidyverse)

# load metadata
metadata = read.delim("../soil_auxo_metadata_full.tsv")

# drep checkm files
drep=read.delim("drep.stats") %>% 
  mutate(Sample_ID=gsub("\\.fa","",Bin.Id)) %>% 
  mutate(Sample_ID=gsub("_.*$","",Sample_ID)) %>% 
  mutate(method=gsub("^.*_","",Bin.Id)) %>% 
  mutate(method=gsub("\\.fa","",method)) %>% 
  mutate(method=ifelse(grepl("lane",method),"MEGAHIT",ifelse(grepl("shovill",method),"shovill","metaGEM"))) %>%
  left_join(.,metadata)

ggplot(drep) + geom_point(aes(y=Completeness,x=Contamination+0.0001,color=method,shape=Type),stat="identity") + scale_x_log10()


ggplot(drep) + geom_boxplot(aes(y=Contamination)) +facet_wrap(~method) + scale_y_log10() -> boxplot_cont
ggplot(drep) + geom_boxplot(aes(y=Completeness)) +facet_wrap(~method)-> boxplot_comp
ggarrange(boxplot_cont,boxplot_comp,nrow = 2)

#drep winners file
drep_winners=read.delim("drep_winners.stats")%>% 
  mutate(Sample_ID=gsub("\\.fa","",genome)) %>% 
  mutate(Sample_ID=gsub("_.*$","",Sample_ID)) %>% 
  mutate(method=gsub("^.*_","",genome)) %>% 
  mutate(method=gsub("\\.fa","",method)) %>% 
  mutate(method=ifelse(grepl("lane",method),"MEGAHIT",ifelse(grepl("shovill",method),"shovill","metaGEM"))) %>%
  filter(genome!="low") %>% left_join(.,metadata)


ggplot(drep_winners) + geom_bar(aes(y=Sample_ID,x=score,fill=method),stat="identity") + facet_wrap(~Type,ncol=1,scales="free")

drep_winners %>% select(genome,method,Type) %>% ggplot() + geom_bar(aes(method,fill=method),stat="count")

# 5 samples with more than 1 MAG generated!
drep_winners %>% group_by(sample) %>% summarize(count=n()) %>% unique()

#drep scores
drep_scores=read.delim("drep_scores.stats")

#check by how much drepping improves, i.e. winners vs not per sample
drep %>% 
  mutate(winner=ifelse(Bin.Id %in% drep_winners$genome,"TRUE","FALSE")) %>% 
  mutate(winner=as.logical(winner)) %>% 
  mutate(quality=ifelse(Contamination<=5&Completeness>=90,"HQ","MQ")) %>%
  left_join(drep_scores) %>%
  mutate(drep_score=ifelse(is.na(drep_score),0,drep_score)) -> drep_df

# nice colors and dodge showing best method and samples where multple bins were reconstructed
drep_df %>%
  group_by(Sample_ID,winner) %>% 
  summarize(ave_drep=ave(drep_score)) %>% 
  unique() %>% pivot_wider(.,names_from = winner,values_from = ave_drep) %>% 
  mutate(diff=`TRUE`-`FALSE`) %>% left_join(.,drep_winners %>% select(Sample_ID,method)) %>%
  ggplot() + 
  geom_bar(aes(x=reorder(Sample_ID,-diff),y=diff,fill=method),stat="identity",position="dodge2",color="black")+ coord_flip() + ylab("Gain in drep score relative to losers") + xlab("Sample") + scale_y_log10()
