library(gridExtra)
library(dplyr)
library(ggplot2)
library(ggpubr)

concoctCheckm = read.delim("data/qc/concoct.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(concoctCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size","set")
concoctBins= read.delim("data/qc/concoct_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(concoctBins) = c("bin","contigs","length")
concoct = left_join(concoctCheckm,concoctBins%>%select(-length),by="bin") %>% filter(contamination<=10,completeness>=50) %>% distinct() %>% select(-set)
concoct$sample = gsub("\\..*$","",concoct$bin)
concoct$binner = "CONCOCT"

metabatCheckm = read.delim("data/qc/metabat.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(metabatCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size","set")
metabatBins= read.delim("data/qc/metabat_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(metabatBins) = c("bin","contigs","length")
metabat = left_join(metabatCheckm,metabatBins%>%select(-length),by="bin") %>% filter(contamination<=10,completeness>=50)%>% distinct() %>% select(-set)
metabat$sample = gsub("\\..*$","",metabat$bin)
metabat$binner = "MetaBAT2"

maxbinCheckm = read.delim("data/qc/maxbin.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(maxbinCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size","set")
maxbinBins= read.delim("data/qc/maxbin_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(maxbinBins) = c("bin","contigs","length")
maxbin = left_join(maxbinCheckm,maxbinBins%>%select(-length),by="bin") %>% filter(contamination<=10,completeness>=50)%>% distinct() %>% select(-set)
maxbin$contigs = as.numeric(maxbin$contigs)
maxbin$sample = gsub("\\..*$","",maxbin$bin)
maxbin$binner = "MaxBin2"

refinedCheckm = read.delim("data/qc/refined.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(refinedCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size","set")
refinedBins= read.delim("data/qc/refined_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(refinedBins) = c("bin","contigs","length")
refined = left_join(refinedCheckm,refinedBins%>%select(-length),by="bin") %>% filter(contamination<=10,completeness>=50)%>% distinct() %>% select(-set)
refined$sample = gsub("\\..*$","",refined$bin)
refined$binner = "metaWRAP_refined"

reassembledCheckm = read.delim("data/qc/reassembled.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(reassembledCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size")
reassembledBins= read.delim("data/qc/reassembled_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(reassembledBins) = c("bin","contigs","length")
reassembled = left_join(reassembledCheckm,reassembledBins%>%select(-length),by="bin") %>% filter(contamination<=10,completeness>=50)%>% distinct()
reassembled$sample = gsub("\\..*$","",reassembled$bin)
reassembled$binner = "metaWRAP_reassembled"
colnames(reassembled)[9] = "Sample_ID"

metadata = read.delim("data/figure_data/soil_auxo_metadata_full.tsv")

#bins <- as.data.frame(matrix(0,nrow = 5,ncol=2))
#colnames(bins) = c("variable","value")
#bins$variable = c("maxbin2","refined","CONCOCT","metabat2","reassembled")
#bins$value = c(as.numeric(dim(maxbin)[1]),as.numeric(dim(refined)[1]),as.numeric(dim(concoct)[1]),as.numeric(dim(metabat)[1]),as.numeric(dim(reassembled)[1]))
#rbind(concoct,metabat,maxbin,refined,reassembled) %>% group_by(binner,sample) %>% summarize(count=n()) -> bins
#colnames(bins)[1]

binplot = ggplot(reassembled %>% 
                   mutate(quality=ifelse(contamination<=5&completeness>=90,"HQ","MQ")) %>% 
                   group_by(Sample_ID,quality) %>% 
                   summarize(count=n()) %>% 
                   left_join(.,metadata),aes(x=reorder(Codes,-count),y=count,fill= quality)) +
  geom_bar(stat = "identity",color="black") +
  ylab("Generated bins") + 
  xlab("Sample") +
  theme(legend.title = element_blank()) +
  ggtitle("Number of bins") + 
  coord_flip() + 
  theme(legend.position = "bottom") + facet_grid(vars(Type),scales="free_y",space="free")

compcont = ggplot(reassembled %>% 
                    mutate(quality=ifelse(contamination<=5&completeness>=90,"HQ","MQ")) %>% 
                    left_join(.,metadata)) + 
  geom_point(aes(completeness,contamination,color=quality,shape=Type))

sizeN50 = ggplot(reassembled %>% 
                   mutate(quality=ifelse(contamination<=5&completeness>=90,"HQ","MQ")) %>% 
                   left_join(.,metadata)) + 
  geom_point(aes(size,N50,color=quality,shape=Type)) +
  scale_y_log10()

contigsN50 = ggplot(reassembled %>% 
                      mutate(quality=ifelse(contamination<=5&completeness>=90,"HQ","MQ")) %>% 
                      left_join(.,metadata)) + 
  geom_point(aes(contigs,N50,color=quality,shape=Type)) +
  scale_y_log10() +
  scale_x_log10()

contigsSize = ggplot(reassembled %>% 
                       mutate(quality=ifelse(contamination<=5&completeness>=90,"HQ","MQ")) %>% 
                       left_join(.,metadata)) + 
  geom_point(aes(contigs,size,color=quality,shape=Type)) +
  scale_y_log10() +
  scale_x_log10()


plot=ggarrange(binplot,ggarrange(compcont,sizeN50,contigsN50,contigsSize,ncol=1,common.legend = TRUE, legend="bottom"),nrow=1)

ggsave("outputs/binningVis.pdf",plot=plot, height = 14, width = 12)


### 

