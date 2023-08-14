
clinical<-read.table('clinical.txt',sep='\t',header=TRUE)
clinical$Group2<-factor(clinical$Group2,levels=c('PR','SD','PD'))
clinical_T1<-subset(clinical,Time=='1')
clinical_T2<-subset(clinical,Time=='2')
clinical_T3<-subset(clinical,Time=='3')

#AtT1(Pretreatment)
clinical_T1<-clinical_T1[order(clinical_T1$Group2),]
patient_order<-clinical_T1$Sample_ID
library(ComplexHeatmap)
orr.assign.1 <- setNames(c("#2196F3","#FFC107","#F44336"),c('PR','SD','PD'))
clini<-HeatmapAnnotation(ORR=clinical_T1$Group2,col=list(ORR=orr.assign.1))
mat<-read.csv('snv.landscape.txt',sep='\t',header=TRUE)
library(tidyr)
library(dplyr)
mat<-subset(mat,gene  %in% patient_order)

mat$gene<-factor(mat$gene,levels=patient_order)
mat<-mat[order(mat$gene),]
head(mat)
write.csv(mat,'snv.landscape.order.csv',row.names=F)
mat<-read.table('snv.landscape.order.csv',sep=",")
mat<-t(mat)
rownames(mat)=mat[,1]
mat<-mat[,-1]
colnames(mat)=mat[1,]
mat<-mat[-1,]
dim(mat)
library(ComplexHeatmap)
library(ggplot2)

myfun<-function(x){
  32-sum(is.na(x))
}
num<-apply(mat,1,myfun)
mat<-cbind(mat,num)
mat<-subset(mat,num>1)

mat<-subset(mat,select=-c(num))
mat
gene_name<-rownames(mat)
library(ComplexHeatmap)
muttype = c("nonsynonymous_SNV","frameshift_deletion","nonframeshift_deletion",
            "frameshift_insertion","nonframeshift_insertion","frameshift_substitution","nonframeshift_substitution","stopgain","stoploss",'UTR','splicing',"UNKNOWN")
col = RColorBrewer::brewer.pal(n = length(muttype), name = 'Paired')
col <- col[c(2, 1, 7, 4, 3, 6, 5, 8, 9, 10,11,12)]
names(col) = muttype

alter_fun= function(x, y, w, h, v) {
  n = sum(v)
  w = w * 0.95
  h = h * 0.9
  grid.rect(x, y, w*0.95, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA))
  if(n){
    grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h,gp = gpar(fill = col[names(which(v))], col = NA), just = "top")
  }
}

pdf('landscape_snv_T1.pdf',width=6,height=8)
oncoPrint(mat,alter_fun = alter_fun,col = col,column_order=patient_order,remove_empty_columns = FALSE, remove_empty_rows = TRUE,row_names_side="left",pct_gp = gpar(fontsize = 7),row_names_gp = gpar(fontsize = 8),column_names_gp=gpar(fontsize = 7),pct_side="right",show_column_names=TRUE,bottom_annotation = clini)
dev.off()

#AtT2(C2D15)
clinical_T2<-clinical_T2[order(clinical_T2$Group2),]
patient_order<-clinical_T2$Sample_ID
library(ComplexHeatmap)
orr.assign.1 <- setNames(c("#2196F3","#FFC107","#F44336"),c('PR','SD','PD'))
clini<-HeatmapAnnotation(ORR=clinical_T2$Group2,col=list(ORR=orr.assign.1))
mat<-read.csv('snv.landscape.txt',sep='\t',header=TRUE)
library(tidyr)
library(dplyr)
mat<-subset(mat,gene  %in% patient_order)

mat$gene<-factor(mat$gene,levels=patient_order)
mat<-mat[order(mat$gene),]
gene_name<-c('APC','TP53','CARD11','ALK','ARID1A','FBXW7','PTPRT','BRAF','EGFR','EPHA5','FOXO1','HMCN1','KMT2C','KRAS','NTRK3','PIK3CA','PREX2','AR','ATRX','CHEK1','DICER1','EPHA7','ESR1','FLCN','GLI1','IRS2','KMT2D','MTOR','NCOR1','NF1','PIK3CG','POLE','PTPRS','SMAD3','SPEN','STAG2','TET1','TOP1')
mat<-mat[,which(names(mat)%in%c('gene',gene_name))]
write.csv(mat,'snv.landscape.order.csv',row.names=F)
mat<-read.table('snv.landscape.order.csv',sep=",")
mat<-t(mat)
rownames(mat)=mat[,1]
mat<-mat[,-1]
colnames(mat)=mat[1,]
mat<-mat[-1,]
dim(mat)
library(ComplexHeatmap)
library(ggplot2)
library(ComplexHeatmap)
muttype = c("nonsynonymous_SNV","frameshift_deletion","nonframeshift_deletion",
            "frameshift_insertion","nonframeshift_insertion","frameshift_substitution","nonframeshift_substitution","stopgain","stoploss",'UTR','splicing',"UNKNOWN")
col = RColorBrewer::brewer.pal(n = length(muttype), name = 'Paired')
col <- col[c(2, 1, 7, 4, 3, 6, 5, 8, 9, 10,11,12)]
names(col) = muttype

alter_fun= function(x, y, w, h, v) {
  n = sum(v)
  w = w * 0.95
  h = h * 0.9
  grid.rect(x, y, w*0.95, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA))
  if(n){
    grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h,gp = gpar(fill = col[names(which(v))], col = NA), just = "top")
  }
}

pdf('landscape_snv_T2.pdf',width=6,height=8)
oncoPrint(mat,alter_fun = alter_fun,col = col,column_order=patient_order,row_order=gene_name,remove_empty_columns = FALSE, remove_empty_rows = FALSE,row_names_side="left",pct_gp = gpar(fontsize = 7),row_names_gp = gpar(fontsize = 8),column_names_gp=gpar(fontsize = 7),pct_side="right",show_column_names=TRUE,bottom_annotation = clini)
dev.off()

#AtT3(C4D15)
clinical_T3<-clinical_T3[order(clinical_T3$Group2),]
patient_order<-clinical_T3$Sample_ID
library(ComplexHeatmap)
orr.assign.1 <- setNames(c("#2196F3","#FFC107","#F44336"),c('PR','SD','PD'))
clini<-HeatmapAnnotation(ORR=clinical_T3$Group2,col=list(ORR=orr.assign.1))
mat<-read.csv('snv.landscape.txt',sep='\t',header=TRUE)
library(tidyr)
library(dplyr)
mat<-subset(mat,gene  %in% patient_order)

mat$gene<-factor(mat$gene,levels=patient_order)
mat<-mat[order(mat$gene),]
gene_name<-c('APC','TP53','CARD11','ALK','ARID1A','FBXW7','PTPRT','BRAF','EGFR','EPHA5','FOXO1','HMCN1','KMT2C','KRAS','NTRK3','PIK3CA','PREX2','AR','ATRX','CHEK1','DICER1','EPHA7','ESR1','FLCN','GLI1','IRS2','KMT2D','MTOR','NCOR1','NF1','PIK3CG','POLE','PTPRS','SMAD3','SPEN','STAG2','TET1','TOP1')
mat<-mat[,which(names(mat)%in%c('gene',gene_name))]
write.csv(mat,'snv.landscape.order.csv',row.names=F)
mat<-read.table('snv.landscape.order.csv',sep=",")
mat<-t(mat)
rownames(mat)=mat[,1]
mat<-mat[,-1]
colnames(mat)=mat[1,]
mat<-mat[-1,]
dim(mat)
library(ComplexHeatmap)
library(ggplot2)
library(ComplexHeatmap)
muttype = c("nonsynonymous_SNV","frameshift_deletion","nonframeshift_deletion",
            "frameshift_insertion","nonframeshift_insertion","frameshift_substitution","nonframeshift_substitution","stopgain","stoploss",'UTR','splicing',"UNKNOWN")
col = RColorBrewer::brewer.pal(n = length(muttype), name = 'Paired')
col <- col[c(2, 1, 7, 4, 3, 6, 5, 8, 9, 10,11,12)]
names(col) = muttype

alter_fun= function(x, y, w, h, v) {
  n = sum(v)
  w = w * 0.95
  h = h * 0.9
  grid.rect(x, y, w*0.95, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA))
  if(n){
    grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h,gp = gpar(fill = col[names(which(v))], col = NA), just = "top")
  }
}

pdf('landscape_snv_T3.pdf',width=6,height=8)
oncoPrint(mat,alter_fun = alter_fun,col = col,column_order=patient_order,row_order=gene_name,remove_empty_columns = FALSE, remove_empty_rows = FALSE,row_names_side="left",pct_gp = gpar(fontsize = 7),row_names_gp = gpar(fontsize = 8),column_names_gp=gpar(fontsize = 7),pct_side="right",show_column_names=TRUE,bottom_annotation = clini)
dev.off()
