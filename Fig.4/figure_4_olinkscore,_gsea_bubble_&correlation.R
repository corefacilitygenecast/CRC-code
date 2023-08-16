#correlation for significant protein and PFS/OS

remove(list=ls())
datacli=read.csv("26out.csv")
datacli=datacli %>% mutate(P_value_sig2=paste0(round(P_value,2),' ',P_value_sig))
theme_heatmap <- 
  theme_bw(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(5, "pt"),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(face='bold'),
        axis.text.x = element_text(angle=45,vjust=1, hjust=1))
theme_heatmap <- 
  theme_bw(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(5, "pt"),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(face='bold'),
        axis.text.x = element_text(angle=45,vjust=1, hjust=1))

ggplot(datacli,aes(x=gene, y=index))+
  geom_point(aes(color = r, size = abs(r))) +
  geom_text(aes(label = P_value_sig2), color = 'black', size = 1) +
  scale_size_area(max_size = 6) +
  ggsci::scale_color_material('teal') +
  theme_heatmap +labs(y='outcome', x ='Protein', title = 'Spearman Correlation between protein and clinical outcome', size = '|rho|', color = 'rho')

#gsea bubble diagram
rm(list = ls())
options(stingsAsFactors = F)
library(OlinkAnalyze)
library(dplyr)
library(ggplot2)
library(stringr)
library(clusterProfiler)
data=read.csv("data1.csv")
ttest_results=olink_ttest(df=data1,variable = "Group3",alternative="two.sided")
gsea_resutls_ttest=olink_pathway_enrichment(data=data1,test_results = ttest_results)
write.csv(gsea_resutls_ttest,"gsea_resutls_ttest.csv")
rt=read.csv("gsea_resutls_ttest.csv")
head(rt)
p=ggplot(rt,aes(NES,ID))+
  geom_point(aes(size=log10.pvalue,color=NES))+coord_cartesian(xlim = c(-2.2,1.8))
p1=p+scale_colour_gradient2(low ="#1C86EE",mid="white",high="#CD0000")+labs(color="NES",size="-log10pvalue",x="NES",y="Pathways")+theme(axis.text.x=element_text(color="black", size=10),axis.text.y=element_text(color="black", size=10)) + 
  theme_bw()+theme(panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank())
p1

# olink score 
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(DESeq2)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(apeglm)))
input <- fread("84_olink_Group_tr_sig.txt",header = T,sep= "\t",stringsAsFactors = F,data.table = F,encoding = "UTF-8",check.names = F)
group <- fread("Group1_condition.xls",header = T,sep ="\t",stringsAsFactors = F,data.table = F,encoding = "UTF-8",check.names=F)

colData <-group
data<-input
## cancel proteins with negative values
input[input<0] = NA
nrow(input)
input <- na.omit(input)
nrow(input)
gene_name<-input$Gene_id
input$Gene_id <-NULL

input <-input*10
input <-round(input, digits = 0)
nrow(input)
factors <-colnames(colData)
design <- eval(parse(text=paste("formula(~",factors,")",sep="")))
design
dds <- DESeqDataSetFromMatrix(countData=input,colData=colData,design=design)
dds <-DESeq(dds)

res.ape <- lfcShrink(dds=dds, coef="groups_PD_vs_non_PD", type="apeglm")
## calculate betaCoeff
beta_conf<-coef(dds, SE = FALSE)

beta_conf <-as.data.frame(beta_conf)

result<-cbind(gene_name,beta_conf$groups_PD_vs_non_PD)
colnames(result) <-c("Gene_id","groups_PD_vs_non_PD")

write.table(result,"betaCoeff.txt",sep="\t", quote=F, row.name=F)

result <-as.data.frame(result)

betaCoeffSubset <-result
##calculate score for each sample
reducedSignature = apply(input*as.numeric(betaCoeffSubset$groups_PD_vs_non_PD), 2, sum)
stand.fun = function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}
reducedGeneSignature = stand.fun(reducedSignature)
score <-as.data.frame(reducedGeneSignature)

write.table(score,"Score_all_sample.txt",sep="\t", quote=F, row.name=T)
##0.33 quantile
reducedGeneSignatureGroup = ifelse(reducedGeneSignature >= as.numeric(quantile(reducedGeneSignature,0.33)), "High", "Low")
Group<-as.data.frame(reducedGeneSignatureGroup)
write.table(Group,"0.33_Score_Group.txt",sep="\t", quote=F, row.name=T)

