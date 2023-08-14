library(ggsci)
library(ggalluvial)
library(ggplot2)

dat01<-read.csv("P7.tsv",header = T,sep = "\t")
dat01<-dat01[,c('mutation_id','sample_id','variant_freq')]
dat01$mutation_id<-as.factor(dat01$mutation_id)
dat01$sample_id<-as.factor(dat01$sample_id)
library(alluvial)
pdf("P7_dynamic.pdf",width=8,height=4)
ggplot(dat01, aes(x=sample_id,y=variant_freq,alluvium= mutation_id))+geom_alluvium(aes(fill = mutation_id, colour = mutation_id),alpha = .75, decreasing = FALSE)+scale_fill_brewer(type = "qual", palette = "Set3") +scale_color_brewer(type = "qual", palette = "Set3")+theme(panel.grid =element_blank(),axis.title.x=element_blank(),panel.background = element_rect(fill = "transparent",colour = NA),axis.line=element_line(colour="black"))
dev.off()