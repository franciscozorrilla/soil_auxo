library(ggplot2)
library(gridExtra)
library(dplyr)

assembly = read.delim("data/qc/assembly.stats",stringsAsFactors = FALSE,header = FALSE,sep = " ")
colnames(assembly) = c("Sample_ID","contigs","length_total")
metadata = read.delim("data/figure_data/soil_auxo_metadata_full.tsv")
assembly %>% left_join(.,metadata) -> assembly

assembly$ave = assembly$length_total/assembly$contigs

aveplot = ggplot(data=assembly) +
  geom_jitter(aes(x=ave,y=Type)) + 
  geom_boxplot(aes(x=ave,y=Type),alpha=0.4,outlier.shape = NA) + 
  xlab("Average contig length") +
  ggtitle("Average contig length across samples") +
  theme(legend.title = element_blank())+ 
  theme(legend.position = "none") + 
  scale_x_log10()

contplot = ggplot(data=assembly) +
  geom_jitter(aes(x=contigs,y=Type)) + 
  geom_boxplot(aes(x=contigs,y=Type),alpha=0.4,outlier.shape = NA) + 
  ggtitle("Contigs across samples") +
  xlab("Contigs across samples") +
  theme(legend.title = element_blank()) + 
  scale_x_log10() + 
  theme(legend.position = "none")

scatplot = ggplot(data=assembly) +
  geom_point(aes(x=length_total,y=contigs,color=Type)) +
  ggtitle("Total length vs number of contigs")+
  xlab("Total length") +
  ylab("Number of contigs") +
  scale_y_log10() + 
  scale_x_log10() +
  theme(legend.position = "bottom")

barplot = ggplot(data=assembly) +
  geom_bar(aes(x=reorder(Codes,-length_total),y=length_total),stat = "identity",color="black",size=0.2,alpha=0.5) +
  ggtitle("Total length across assemblies") + 
  ylab("Length (bp)") +
  xlab("Sample ID") +
  scale_y_log10() +
  coord_flip() + 
  theme(legend.position = "none") + facet_grid(vars(Type),scales="free_y",space="free")

assemblyplot=grid.arrange(barplot,arrangeGrob(scatplot,aveplot,contplot,nrow=3),ncol=2,nrow=1)

ggsave("outputs/assemblyVis.pdf",plot= assemblyplot,device = "pdf",height = 14,width = 12)
