NPX_CRC<-read.csv('NPX_BL_sig.txt',sep='\t',header=TRUE,check.names = F)
cytokine<-colnames(NPX_CRC[2:27])
clinical<-read.csv('clinical.txt',header=TRUE,sep="\t")
NPX_clinical<-merge(x=NPX_CRC,y=clinical,by="Patient_ID",all=FALSE)

NPX_clinical$Group1<-factor(NPX_clinical$Group3,levels=c('non-PD','PD'))
NPX_clinical$Tumor_site<-as.factor(NPX_clinical$Tumor_site)
NPX_clinical$treated_lines<-as.factor(NPX_clinical$treated_lines)
NPX_clinical$number_organs<-as.factor(NPX_clinical$number_organs)
NPX_clinical$Liver_metastasis<-as.factor(NPX_clinical$Liver_metastasis)
NPX_clinical$Lung_metastasis<-as.factor(NPX_clinical$Lung_metastasis)
NPX_clinical$Colon_Rectum<-as.factor(NPX_clinical$Colon_Rectum)

NPX_clinical<-NPX_clinical[order(NPX_clinical$Group1,NPX_clinical$PFS),]
Sample_order<-NPX_clinical$Sample_ID
NPX_clinical_mat<-NPX_clinical[,1:27]
NPX_clinical_mat$Sample_ID<-as.character(NPX_clinical_mat$Sample_ID)
parameters<-colnames(NPX_clinical_mat)
colnames(NPX_clinical_mat)<-NULL
NPX_clinical_mat<-rbind(parameters,NPX_clinical_mat)
NPX_clinical_mat<-t(NPX_clinical_mat)
colnames(NPX_clinical_mat)=NPX_clinical_mat[1,]
NPX_clinical_mat<-NPX_clinical_mat[-1,]
rownames_NPX<-NPX_clinical_mat[,1]
NPX_clinical_mat<-NPX_clinical_mat[,-1]
NPX_clinical_mat<-apply(NPX_clinical_mat,2,as.numeric)
rownames(NPX_clinical_mat)<-rownames_NPX
exp<-apply(NPX_clinical_mat,1,scale)
rownames(exp)<-colnames(NPX_clinical_mat)
exp <- t(exp)


library(ComplexHeatmap)
library(circlize)
col_anno<-HeatmapAnnotation(Group3=NPX_clinical$Group3,Tumor_site=NPX_clinical$Tumor_site,treated_lines=NPX_clinical$treated_lines,number_organs=NPX_clinical$number_organs,Liver_metastasis=NPX_clinical$Liver_metastasis,Lung_metastasis=NPX_clinical$Lung_metastasis,PFS=NPX_clinical$PFS)
pdf('NPX_BL_heatmap_sig.pdf',width=8,height=8)
Heatmap(exp,col=colorRamp2(c(-2,0,2), c("blue", "#EEEEEE", "red")),cluster_rows = FALSE,cluster_columns=FALSE,column_split=NPX_clinical$Group3,rect_gp = gpar(col = "white", lwd = 2),bottom_annotation = col_anno)
dev.off()
