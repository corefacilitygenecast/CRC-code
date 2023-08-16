#correlation matrix
data=read.csv("26cf1-1.csv")
library(corrplot)
library(ggcorrplot)
library(eyedroppeR)
str(data)
set.seed(1000)
path="c.jpg"
extract_pal(1, path, label = "heatmap", sort = "auto")
cor_matrix=round(cor(data,use="na.or.complete",method = "pearson"),2)
p.mat_matrix=round(cor_pmat(data,method="pearson"),2)
col=COL2('RdBu', 100)
col=rev(col)
corrplot(cor_matrix,method = "circle",type="lower",order ="AOE",tl.cex = 0.5,cl.cex = 0.5,
         tl.col = "black",tl.srt = 45,col=col,p.mat = p.mat_matrix, insig = "label_sig", 
         sig.level = c(.01, .05), pch.cex = 0.2, pch.col = "black")

# univariate cox analysis
remove(list=ls())
library(ezcox)
library(survival)
data=read.csv("key.csv")
results=ezcox(data,time="OS",status = "Osyes",covariates = c("Age","Gender","Tumor_site", "treated_lines","number_organs", "Liver_metastasis","Lung_metastasis","Liver_lung_metastasis","Colon_Rectum","CF_1_TMB","CF_2_TMB","CF_3_TMB","CF_1_meanvaf","CF_2_meanvaf","CF_3_meanvaf","delta_meanvaf_CF_2","delta_meanvaf_CF_3","olink_T1","olinkT1_score","olink_T2","olink_T3","olink_T2b","olink_T3b","olinkT2_score","olinkT3_score","olinkT2b_score","olinkT3b_score"))
write.csv(results,"singlecox.csv")

#ggplot for univariate cox analysis
remove(list=ls())
data=read.csv("singlecox.csv")
library(ggplot2)
colnames(data)
data$ID=factor(data$ID,levels = rev(data$ID))
single=ggplot(data,aes(color=p.value,x=HR,xmin=lower_95,xmax=upper_95,y=ID))+
  geom_pointrange(size=0,fatten = 0)+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.y = element_text(color = rev(data$text.color)))+geom_pointrange(
          size=1, 
          fatten = 2)+ 
  geom_vline(xintercept = 1, color = "black",size = 0.4) 

#multivariate cox analysis
remove(list=ls())
data4=read.csv("finalmodel.csv")
colnames(data4)
data4$olinkT1_score=factor(data4$olinkT1_score,levels = c("low","high"))
multimodel=coxph(Surv(OS,Osyes)~.,data=data4)
cox=summary(multimodel)
cox$coefficients
cox$conf.int
mul_HR=round(cox$coefficients[,2],3)
mul_Pvalue=round(cox$coefficients[,5],4)
mul_CI1=round(cox$conf.int[,3],2)
mul_CI2=round(cox$conf.int[,4],2)
mul_CI95=paste(mul_CI1,'-',mul_CI2)
mul_cox=data.frame("HR"=mul_HR,"lower95"=mul_CI1,"up95"=mul_CI2,"Pvalue"=mul_Pvalue)
write.csv(mul_cox,"multiplecox.csv")

#ggplot for multivariate cox analysis
remove(list=ls())
data=read.csv("multiplecox.csv")
library(ggplot2)
colnames(data)
data$ID=factor(data$ID,levels = rev(data$ID))
multiple=ggplot(data,aes(color=Pvalue,x=HR,xmin=lower95,xmax=up95,y=ID))+
  geom_pointrange(size=0,fatten = 0)+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.y = element_text(color = rev(data$text.color)))+geom_pointrange(
          size=1, 
          fatten = 2)+
  geom_vline(xintercept = 1, color = "black",size = 0.4) 

#multivariate cox model prediction
remove(list=ls())
data=read.csv("finalmodelmodify.csv")
multicox=coxph(Surv(OS,Osyes)~.,data = data)

#risk score
riskscore=predict(multicox,type = "risk", newdata = data)
riskscore=as.data.frame(riskscore)
x
# riskscore KM
library(dplyr)
data$riskscore=riskscore$riskscore
data$riskscore2 <- ifelse(data$riskscore > median(data$riskscore),
                          "High","Low")
write.csv(data,"multimodelwithrisk.csv")

fit <- survfit(Surv(OS, Osyes) ~ riskscore2, data=data)
ggsurvplot(fit,data = data,palette =c("red", "#27408B"),surv.median.line ="hv",
           risk.table = T,
           xlab="Time(months)",
           ylab="OS",
           pval = TRUE, title="finalmodelos")
sdiff=survdiff(Surv(OS, Osyes) ~ riskscore2, data=data)
p.val = 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
data.survdiff=sdiff
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))

#pfs
remove(list=ls())
data=read.csv("multimodelpfs.csv")
colnames(data)
multicox=coxph(Surv(PFS,PFSyes)~.,data = data)
#risk score
riskscore=predict(multicox,type = "risk", newdata = data)
riskscore=as.data.frame(riskscore)

# riskscore KM
library(dplyr)
data$riskscore=riskscore$riskscore
data$riskscore2 <- ifelse(data$riskscore > median(data$riskscore),
                          "High","Low")
write.csv(data,"multimodelPFSwithrisk.csv")

fit <- survfit(Surv(PFS, PFSyes) ~ riskscore2, data=data)
ggsurvplot(fit,data = data,palette =c("red", "#27408B"),surv.median.line ="hv",
           risk.table = T,
           xlab="Time(months)",
           ylab="PFS",
           pval = TRUE, title="finalmodelPFs")
sdiff=survdiff(Surv(PFS, PFSyes) ~ riskscore2, data=data)
p.val = 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
data.survdiff=sdiff
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))


# AUC value
library(timeROC)
data$age_num <- 0
data$age_num[which(data$Age == "more_58")] <- 1

timeres_age <- timeROC(data$OS, data$Osyes, -data[,13], cause = 1, times = 12,iid = T)
timeres_treatlines <- timeROC(data$OS, data$Osyes, -data[,4], cause = 1, times = 12,iid = T)
timeres_CF_1 <- timeROC(data$OS, data$Osyes, data[,5], cause = 1, times = 12,iid = T)
timeres_delta_CF_2 <- timeROC(data$OS, data$Osyes, data[,6], cause = 1, times = 12,iid = T)
timeres_delta_CF_3 <- timeROC(data$OS, data$Osyes, data[,7], cause = 1, times = 12,iid = T)

data$olinkT1_score_num <- 0
data$olinkT1_score_num[which(data$olinkT1_score == "high")] <- 1
timeres_olinkt1_score <- timeROC(data$OS, data$Osyes, data[,11], cause = 1, times = 12, iid = T)

timeres_riskscore <- timeROC(data$OS, data$Osyes, data[,9], cause = 1, times = 12, iid = T)

data$riskscore2_num <- 0
data$riskscore2_num[which(data$riskscore2 == "High")] <- 1
timeres_riskscore2 <- timeROC(data$OS, data$Osyes, data[,12], cause = 1, times = 12, iid = T)




