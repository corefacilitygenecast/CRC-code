library(haven)
library(tidyverse)
library(reshape2)

x<-read.table("merge_plot.xls",header=T,sep="\t")
data <- x %>%
  mutate(text = case_when(  
    pvalue == 0 ~ paste(round(Correlation, 2), "\n "), # round() 只保留两位小数
    pvalue <= 0.05 & pvalue >0.01 ~ paste(round(Correlation,2), "\n*"),
    pvalue <= 0.01 & pvalue >0.001 ~ paste(round(Correlation,2), "\n**"),
    pvalue <= 0.001 ~ paste(round(Correlation,2), "\n***"),
pvalue > 0.05 ~ paste(round(Correlation,2), "\n ")))


data$new_cell1 <- factor(data$cell1,levels = c("IL8","TNFRSF9","TIE2","MCP_3","CD40_L","CD244","EGF","ANGPT1","IL7","PGF","IL6","ADGRG1","MCP_1","CRTAM","CXCL11","MCP_4","TRAIL","FGF2","CXCL9","CD8A","CAIX","MUC_16","ADA","CD4","Gal_9","VEGFR_2","CD40","IL18","GZMH","KIR3DL1","LAP_TGF_beta_1","CXCL1","TNFSF14","TWEAK","PDGF_subunit_B","PDCD1","FASLG","CD28","CCL19","MCP_2","CCL4","IL15","Gal_1","PD_L1","CD27","CXCL5","IL5","HGF","GZMA","HO_1","CX3CL1","CXCL10","CD70","IL10","TNFRSF12A","CCL23","CD5","CCL3","MMP7","ARG1","NCR1","DCN","TNFRSF21","TNFRSF4","MIC_A_B","CCL17","ANGPT2","PTN","IFN_gamma","LAMP3","CASP_8","ICOSLG","MMP12","CXCL13","PD_L2","VEGFA","LAG3","CCL20","TNF","KLRD1","GZMB","CD83","IL12","CSF_1","CF_1_meanvaf","CF_2_meanvaf","CF_3_meanvaf","CF_4_meanvaf","CF_1_hGE","CF_2_hGE","CF_3_hGE","CF_4_hGE","delta_hGE_CF_2","delta_hGE_CF_3","delta_hGE_CF_4","delta_meanvaf_CF_2","delta_meanvaf_CF_3","delta_meanvaf_CF_4"))
data$new_cell2 <- factor(data$cell2,levels = rev(c("IL8","TNFRSF9","TIE2","MCP_3","CD40_L","CD244","EGF","ANGPT1","IL7","PGF","IL6","ADGRG1","MCP_1","CRTAM","CXCL11","MCP_4","TRAIL","FGF2","CXCL9","CD8A","CAIX","MUC_16","ADA","CD4","Gal_9","VEGFR_2","CD40","IL18","GZMH","KIR3DL1","LAP_TGF_beta_1","CXCL1","TNFSF14","TWEAK","PDGF_subunit_B","PDCD1","FASLG","CD28","CCL19","MCP_2","CCL4","IL15","Gal_1","PD_L1","CD27","CXCL5","IL5","HGF","GZMA","HO_1","CX3CL1","CXCL10","CD70","IL10","TNFRSF12A","CCL23","CD5","CCL3","MMP7","ARG1","NCR1","DCN","TNFRSF21","TNFRSF4","MIC_A_B","CCL17","ANGPT2","PTN","IFN_gamma","LAMP3","CASP_8","ICOSLG","MMP12","CXCL13","PD_L2","VEGFA","LAG3","CCL20","TNF","KLRD1","GZMB","CD83","IL12","CSF_1","CF_1_meanvaf","CF_2_meanvaf","CF_3_meanvaf","CF_4_meanvaf","CF_1_hGE","CF_2_hGE","CF_3_hGE","CF_4_hGE","delta_hGE_CF_2","delta_hGE_CF_3","delta_hGE_CF_4","delta_meanvaf_CF_2","delta_meanvaf_CF_3","delta_meanvaf_CF_4")))
pdf("Extended_Fig.6.pdf",width=33,height=27)
ggplot(data, aes(new_cell1,new_cell2)) + 
  geom_tile(aes(fill = Correlation), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "blue",mid = "white",high = "red")+
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() + # 不要背景
  theme(axis.title.x=element_blank(), # 去掉 title
        axis.ticks.x=element_blank(), # 去掉x 轴
        axis.title.y=element_blank(), # 去掉 y 轴
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"), # 调整x轴文字，字体加粗
        axis.text.y = element_text(size = 14, face = "bold")) + #调整y轴文字
  labs(fill =paste0(" *  p <= 0.05","\n\n","**  p <= 0.01","\n\n","***  p <= 0.001","\n\n","Correlation"))    # 修改 legend 内容
dev.off()

