NPX_CRC<-read.csv('NPX_CRC.txt',sep='\t',header=TRUE,check.names = F)
cytokine<-colnames(NPX_CRC[2:85])
clinical<-read.csv('clinical.txt',header=TRUE,sep="\t")
NPX_clinical<-merge(x=NPX_CRC,y=clinical,by="Patient_ID",all=FALSE)
NPX_clinical$Group2<-factor(NPX_clinical$Group2,levels=c('PR','SD','PD'))
Sample_order<-c('P28','P14','P7','P4','P13','P23','P32','P33','P22','P12','P18','P19','P5','P16','P30','P8','P2','P3','P29','P35','P17','P36','P25','P11','P10','P1','P9','P34','P20','P6','P27','P24')
NPX_clinical$Patient_ID<-factor(NPX_clinical$Patient_ID,levels=Sample_order)#自定义排序
NPX_clinical<-NPX_clinical[order(NPX_clinical$Patient_ID),]
NPX_clinical_mat<-NPX_clinical[,1:85]
NPX_clinical_mat$Patient_ID<-as.character(NPX_clinical_mat$Patient_ID)
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
col_anno<-HeatmapAnnotation(Group2=NPX_clinical$Group2,PFS=NPX_clinical$PFS)

pdf('NPX_BL_heatmap.pdf',width=12,height=20)
Heatmap(exp,col=colorRamp2(c(-2,0,2), c("blue", "#EEEEEE", "red")),cluster_rows = FALSE,cluster_columns=FALSE,column_order=Sample_order,rect_gp = gpar(col = "white", lwd = 2),bottom_annotation = col_anno)
dev.off()

