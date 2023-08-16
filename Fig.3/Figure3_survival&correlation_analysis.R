# for survival analysis
data=read.csv("vaf.csv")
library(survival)
library(survminer)

dvaf2os=survfit(Surv(data$OS_months,data$OS_status)~delta.vaf2_status,data=data)
ggsurvplot(dvaf2os,data = data,palette =c("#EE7600", "#27408B"),surv.median.line ="hv",
           risk.table = T,
           xlab="Time(months)",
           ylab="OS",
           pval = TRUE)
sdiff=survdiff(Surv(data$OS_months,data$OS_status)~delta.vaf2_status,data  = data)
p.val = 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
data.survdiff=sdiff
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 

dvaf2pfs=survfit(Surv(data$PFS_months_x_x,data$PFS_status_x_x)~delta.vaf2_status,data=data)
ggsurvplot(dvaf2pfs,data = data,palette =c("#EE7600", "#27408B"),surv.median.line ="hv",
           risk.table = T,
           xlab="Time(months)",
           ylab="PFS",
           pval = TRUE)
sdiff=survdiff(Surv(data$PFS_months_x_x,data$PFS_status_x_x)~delta.vaf2_status,data  = data)
p.val = 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
data.survdiff=sdiff
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 

#for correlation analysis
getwd()
cor=read.csv("correlation.csv")
df1=as.data.frame(cor)
corT=cor.test(cor$deltasize,cor$delta_meanvaf_CF_2, method = "pearson")
cor1=corT$estimate
pValue=corT$p.value
library(ggplot2)
library(ggpubr)
library(ggExtra)
x=cor$deltasize
y=cor$delta_meanvaf_CF_2
p1=ggplot(df1, aes(x, y)) + 
  xlab(x)+ylab(y)+
  geom_point(shape = 21, colour = "#1874CD", fill = "#1874CD", size = 5, stroke = .5,alpha=0.8)+ geom_smooth(method="lm",formula = y ~ x,linetype=2,color="#6495ED",fill="#D3D3D3") + theme_bw()+
  stat_cor(method = 'pearson', aes(x =x, y =y))
p1