## 
### ---------------
###
### Create: Noe Lee
### Email: 2021585592@qq.com
###
### ---------------

rm(list = ls())
options(stringsAsFactors = F)
setwd("./HCC")
pdata_raw<-read.table(file="LIHC_survival.txt.gz",sep="\t")
head(pdata_raw)
rownames(pdata_raw)<-pdata_raw[,1]
pdata_raw<-t(pdata_raw)
rownames(pdata_raw)<-pdata_raw[,1]
pdata_raw<-pdata_raw[-1,-1]
pdata_raw<-t(pdata_raw)
save(pdata_raw,file="pdata_raw.Rdata")
load("step1-output.Rdata")
pat<-colnames(exprSet)
pat<-substr(pat,1,12)
pat<-gsub(".","-",pat,fixed = T)
table(pdata_raw[,1]%in%pat)
pdata<-pdata_raw[pdata_raw[,1]%in%pat,]
table(duplicated(pdata[,1]))
pdata<-pdata[!duplicated(pdata[,1]),]
colnames(exprSet)<-pat
exprSet<-exprSet[,colnames(exprSet)%in%pdata[,1]]
exprSet<-exprSet[,!duplicated(colnames(exprSet))]

rownames(group_list)<-pat
group_list<-group_list[rownames(group_list)%in%pdata[,1],]
group_list<-group_list[!duplicated(rownames(group_list)),]
pdata<-as.data.frame(pdata)
pdata$OS.time<-round(as.numeric(pdata$OS.time)/30,3)
pdata$OS<-as.numeric(pdata$OS)
pdata$group<-group_list[,2]

save(exprSet,pdata,file="survival.Rdata")

## survival analysis
library(survival)
library(survminer)
load("survival.Rdata")
# 利用ggsurvplot快速绘制漂亮的生存曲线图
sfit <- survfit(Surv(OS.time, OS)~group, data=pdata)
sfit
summary(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)
## 挑选感兴趣的基因做差异分析
pdata$TSPAN6=ifelse(exprSet['TSPAN6',]>median(exprSet['TSPAN6',]),'high','low')
table(pdata$TSPAN6)
ggsurvplot(survfit(Surv(OS.time, OS)~TSPAN6, data=pdata), conf.int=F, pval=TRUE)

pdata$SCYL3=ifelse(exprSet['SCYL3',]>median(exprSet['SCYL3',]),'high','low')
table(pdata$SCYL3)
ggsurvplot(survfit(Surv(OS.time, OS)~SCYL3, data=pdata), conf.int=F, pval=TRUE)
## 批量生存分析 使用  logrank test 方法
mySurv=with(pdata,Surv(OS.time, OS))
log_rank_p <- apply(exprSet , 1 , function(gene){
  #gene=exprSet[1,]
  pdata$group=ifelse(gene>median(gene),'high','low')  
  data.survdiff=survdiff(mySurv~group,data=pdata)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})
log_rank_p=sort(log_rank_p)
head(log_rank_p)
boxplot(log_rank_p) 
pdata$ZNF709=ifelse(exprSet['ZNF709',]>median(exprSet['ZNF709',]),'high','low')
table(pdata$ZNF709)
ggsurvplot(survfit(Surv(OS.time, OS)~ZNF709, data=pdata), conf.int=F, pval=TRUE)


## 批量生存分析 使用 coxph 回归方法
colnames(pdata)
mySurv=with(pdata,Surv(OS.time, OS))
cox_results <-apply(exprSet , 1 , function(gene){
  group=ifelse(gene>median(gene),'high','low')
  survival_dat <- data.frame(group=group,stringsAsFactors = F)
  m=coxph(mySurv ~ group, data =  survival_dat)
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['grouplow',])
  
})
cox_results=t(cox_results)
table(cox_results[,4]<0.05)
cox_results[cox_results[,4]<0.05,]

length(setdiff(rownames(cox_results[cox_results[,4]<0.05,]),
               names(log_rank_p[log_rank_p<0.05])
))
length(setdiff( names(log_rank_p[log_rank_p<0.05]),
                rownames(cox_results[cox_results[,4]<0.05,])
))
length(unique( names(log_rank_p[log_rank_p<0.05]),
               rownames(cox_results[cox_results[,4]<0.05,])
))
save(log_rank_p,cox_results,exprSet,pdata,file = 'survival_analysis.Rdata')
