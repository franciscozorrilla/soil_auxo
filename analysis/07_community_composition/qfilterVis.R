library(gridExtra)
library(dplyr)
library(ggplot2)

qfilter = read.delim("data/qc/qfilter.stats",stringsAsFactors = FALSE, header = FALSE, sep = " ")
colnames(qfilter) = c("Sample_ID","readsBF","readsAF","basesBF","basesAF","percentReads","q20BF","q20AF","q30BF","q30AF")
metadata = read.delim("data/figure_data/soil_auxo_metadata_full.tsv")
qfilter %>% left_join(.,metadata) -> qfilter

reads = ggplot(data = qfilter) + 
  geom_jitter(aes(x=readsAF,y=Type)) + 
  geom_boxplot(aes(x=readsAF,y=Type),alpha=0.4,outlier.shape = NA) + 
  ggtitle("Post-flitering number of reads across samples") + 
  xlab("Total reads") + 
  theme(legend.title = element_blank()) + scale_x_log10()

bases = ggplot(data = qfilter) + 
  geom_jitter(aes(x=basesAF,y=Type)) + 
  geom_boxplot(aes(x=basesAF,y=Type),alpha=0.4,outlier.shape = NA) + 
  ggtitle("Post-flitering number of bases across samples") + 
  xlab("Total reads") + 
  theme(legend.title = element_blank()) + scale_x_log10()

q20 = ggplot(data = qfilter) + 
  geom_jitter(aes(x=q20AF*100,y=Type)) + 
  geom_boxplot(aes(x=q20AF*100,y=Type),alpha=0.4,outlier.shape = NA) +
  ggtitle("Post-filtering percent Q20 bases across samples") + 
  xlab("Percent of Q20 bases") + 
  theme(legend.title = element_blank())

q30 = ggplot(data = qfilter) + 
  geom_jitter(aes(x=q30AF*100,y=Type)) + 
  geom_boxplot(aes(x=q30AF*100,y=Type),alpha=0.4,outlier.shape = NA) +
  ggtitle("Post-filtering percent Q30 bases across samples") + 
  xlab("Percent of Q30 bases") + 
  theme(legend.title = element_blank())

bar = ggplot(data = qfilter) +
  geom_bar(aes(x=reorder(Codes,-basesBF),y=basesBF,fill="Pre-filtering"),color="black",stat = "identity") + 
  geom_bar(aes(x=reorder(Codes,-basesBF),y=basesAF,fill="Post-filtering"),color="black",stat = "identity") + 
  geom_bar(aes(x=reorder(Codes,-basesBF),y=q20AF*basesAF,fill="Q20 bases"),color="black",stat = "identity") + 
  geom_bar(aes(x=reorder(Codes,-basesBF),y=q30AF*basesAF,fill="Q30 bases"),color="black",stat = "identity") + 
  coord_flip() +
  ggtitle("QC summary stacked bar plot") + 
  xlab("Sample ID") + 
  ylab("Base pairs") +
  theme(legend.title = element_blank(),legend.position = "bottom")+ facet_grid(vars(Type),scales="free_y",space="free")

qfilt=grid.arrange(bar,arrangeGrob(reads,bases,q20,q30,nrow=4,ncol=1),ncol =2,nrow=1)
ggsave("outputs/qfilterVis.pdf",plot= qfilt,device = "pdf",height = 14, width=12)
