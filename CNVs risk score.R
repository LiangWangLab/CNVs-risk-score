library(ggplot2)
library(tidyverse)
library(dplyr)
library(knitr)
library(survminer)
library(survival)
library(ggpubr)
####
normalized_counts <- read.table("Normalized_count_All.txt",header = T)
GOI<-read.table("Genes of interest.txt")

normalized_counts <- data.frame(t(normalized_counts))
normalized_counts <- normalized_counts[,GOI$V1]
#normalized_counts_log <- log2(normalized_counts)
control <- normalized_counts[grepl("CPD",rownames(normalized_counts)),]
control_median <- apply(control,2,median)
patient <- normalized_counts[!grepl("CPD",rownames(normalized_counts)),]
patient_ratio <- data.frame(t(apply(patient, 1, function(x) x/control_median)))
patient_ratio <- log2(patient_ratio)
#Geneinfo
Geneinfo <- read.table("Geneinfo.txt",header = T)
Geneinfo <- Geneinfo %>%slice(match(colnames(normalized_counts), Gene))
nGene <- nrow(Geneinfo)
amp <- which(Geneinfo$CNV=="Amplification")
del <- which(Geneinfo$CNV=="Deletion")
cutoff=0.3

promote_visit1 <- read.table("promote_visit1_meta.txt",header = T)
promote_visit1$OSMonth <- promote_visit1$OSDay/30
promote_visit1$PFSMonth <- promote_visit1$PFSDay/30
promote_visit2 <- read.table("promote_visit2_meta.txt",header = T)
promote_visit2$OSMonth <- promote_visit2$OSDay/30
promote_visit2$PFSMonth <- promote_visit2$PFSDay/30
SDMS_visit1 <- read.table("SDMS_visit1_meta.txt",header = T)
SDMS_visit1$OSMonth <- SDMS_visit1$OSDay/30
SDMS_visit2 <- read.table("SDMS_visit2_meta.txt",header = T)
SDMS_visit2$OSMonth <- SDMS_visit2$OSDay/30

cohort <- list(pv1OS=promote_visit1, pv2OS=promote_visit2, sv1OS=SDMS_visit1,sv2OS=SDMS_visit2,pv1PFS=NULL,pv2PFS=NULL)
cohortname<-c("PROMOTE pre-treatment","PROMOTE post-treatment","HR","HR visit 2")
####Waterfall####
for (l in 1:4){
  sample <- nrow(cohort[[l]])
  meta <- cohort[[l]]
  ratio <- patient_ratio[cohort[[l]]$PatientID,]
  #AMP
  for (i in amp) {
    png(filename = paste("CNV waterfall/",names(cohort)[l],colnames(ratio)[i],"copy number variation.png"),width=800,height = 800)
    barplot(sort(ratio[,i], decreasing = T),cex.axis=1.6, cex.lab=1.6,cex.main=2,
            col="blue", space=0.5, ylim=c(-1,2), axisnames = F, ylab="log2 ratio",
            main=paste(colnames(ratio)[i],"copy number variation"))
    abline(h = cutoff,col="red", lwd=2, lty=2)
    percentage <- sum(ratio[,i]>cutoff)/nrow(ratio)*100
    percentage <- round(percentage,2)
    text(sample-5,1,labels =paste("Amplification =",percentage,"%"),cex = 1.6) 
    dev.off()
  }
  #DEL
  for (i in del) {
    png(filename = paste("CNV waterfall/",names(cohort)[l],colnames(ratio)[i],"copy number variation.png"),width=800,height = 800)
    barplot(sort(ratio[,i], decreasing = T),cex.axis=1.6, cex.lab=1.6,cex.main=2,
            col="blue", space=0.5, ylim=c(-1,1), axisnames = F, ylab="log2 ratio",
            main=paste(colnames(ratio)[i],"copy number variation"))
    abline(h = -cutoff,col="red", lwd=2, lty=2)
    percentage <- sum(ratio[,i]<(-cutoff))/nrow(ratio)*100
    percentage <- round(percentage,2)
    text(sample-5,0.5,labels =paste("Deletion =",percentage,"%"),cex = 1.6) 
    dev.off()
  }
}
####indel####
cohort_indel <-list()
for (l in 1:4){
  sample <- nrow(cohort[[l]])
  meta <- cohort[[l]]
  ratio <- patient_ratio[cohort[[l]]$PatientID,]
  
  indel <- matrix("NC",nrow(ratio),nGene)
  indel <- data.frame(indel)
  rownames(indel)<-rownames(ratio)
  colnames(indel)<-colnames(ratio)
  for (i in 1:nrow(ratio)) {
    for (j in 1:ncol(ratio)) {
      if (j %in% amp & ratio[i,j]>cutoff) indel[i,j]="Amp"
      if (j %in% del & ratio[i,j]<(-cutoff)) indel[i,j]="Del"
    }
  }
  indel <- cbind(indel,meta)
  cohort_indel[[l]]<-indel
}
names(cohort_indel)=names(cohort)[1:4]
# write.table(cohort_indel[[1]],"pv1_indel.txt",col.names = NA, sep="\t",quote = F)
# write.table(cohort_indel[[2]],"pv2_indel.txt",col.names = NA, sep="\t",quote = F)
# write.table(cohort_indel[[3]],"sv1_indel.txt",col.names = NA, sep="\t",quote = F)
# write.table(cohort_indel[[4]],"sv2_indel.txt",col.names = NA, sep="\t",quote = F)
####Single gene survival####
for (l in 1:4){
  indel <- cohort_indel[[l]]
  median_survival_OS <- list()#run l=1 for PROMOTE median OS
  HR_OS <- list()
  for (g in 1:nGene) {
    df <- data.frame(OSMonth=indel$OSMonth,Death=indel$Death, CNV=indel[,g])
    if (g %in% amp) {df$CNV<- factor(df$CNV, levels = c("NC","Amp"))}
    if (g %in% del) {df$CNV<- factor(df$CNV, levels = c("NC","Del"))}
    cox <- coxph(Surv(OSMonth, Death) ~ CNV, data=df)
    HR <- round(exp(cox$coefficients),2)
    HR_OS[[g]] <- summary(cox)$conf.int
    OSfit <- survfit(Surv(OSMonth, Death) ~ CNV, data=df)
    table<-data.frame(summary(OSfit)$table)
    median_survival_OS[[g]]<-table
    png(file=paste("Single gene Survival/",cohortname[l]," ",colnames(ratio)[g],"-OS.png",sep = ""))
    OS_plot <- ggsurvplot(OSfit, pval=TRUE, 
                          xlab = "Time in months",title=paste(cohortname[l],"OS -",colnames(ratio)[g]),
                          legend="bottom", legend.title = colnames(ratio)[g],
                          risk.table = T, tables.theme = theme_cleantable(), risk.table.fontsize=5,
                          font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14)
    OS_plot$plot <- OS_plot$plot+annotate("text", x = 0, y = 0.1, label = paste("Hazard ratio = ",HR,sep = ""),size=5,hjust = 0)
    print(OS_plot)
    dev.off()
  }
  names(median_survival_OS)=colnames(ratio)
  names(HR_OS)=colnames(ratio)
}
for (l in 1:2){
  indel <- cohort_indel[[l]]
  median_survival_PFS <- list()#run l=1 for PROMOTE median PFS
  HR_PFS <- list()
  for (g in 1:nGene) {
    df <- data.frame(PFSMonth=indel$PFSMonth,Progression=indel$Progression, CNV=indel[,g])
    if (g %in% amp) {df$CNV<- factor(df$CNV, levels = c("NC","Amp"))}
    if (g %in% del) {df$CNV<- factor(df$CNV, levels = c("NC","Del"))}
    cox <- coxph(Surv(PFSMonth, Progression) ~ CNV, data=df)
    HR <- round(exp(cox$coefficients),2)
    HR_PFS[[g]] <- summary(cox)$conf.int
    PFSfit <- survfit(Surv(PFSMonth, Progression) ~ CNV, data=df)
    table<-data.frame(summary(PFSfit)$table)
    median_survival_PFS[[g]]<-table
    png(file=paste("Single gene Survival/",cohortname[l]," ",colnames(ratio)[g],"-PFS.png",sep = ""))
    PFS_plot <- ggsurvplot(PFSfit, pval=TRUE, 
                           xlab = "Time in months",title=paste(cohortname[l],"PFS -",colnames(ratio)[g]),
                           legend="bottom", legend.title = colnames(ratio)[g],
                           risk.table = T, tables.theme = theme_cleantable(), risk.table.fontsize=5,
                           font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14)
    PFS_plot$plot <- PFS_plot$plot+annotate("text", x = 0, y = 0.1, label = paste("Hazard ratio = ",HR,sep = ""),size=5,hjust = 0)
    print(PFS_plot)
    dev.off()
  }
  names(median_survival_PFS)=colnames(ratio)
  names(HR_PFS)=colnames(ratio)
}

#####Multi-genes survival#####
N=8
cohortname<-c("PROMOTE","PROMOTE","RWHR","RWHR","PROMOTE","RWHR")
####Overall Survival####
cox_result_list <- list()
risk.score_list <- list()
for (l in 1:4){
  indel <- cohort_indel[[l]]
  sample <- nrow(indel)
  cox_result <-  data.frame(coefficients=rep(NA,nGene), HazardRatio=rep(NA,nGene), pvalue=rep(NA,nGene), row.names = colnames(indel[1:nGene]))
  
  for (i in 1:nGene) {
    df <-data.frame(OSMonth=indel$OSMonth,Death=indel$Death,CNV=indel[,i])
    if (i %in% amp) {df$CNV<- factor(df$CNV, levels = c("NC","Amp"))}
    if (i %in% del) {df$CNV<- factor(df$CNV, levels = c("NC","Del"))}
    result <- coxph(Surv(OSMonth, Death) ~ CNV, data=df)
    summary(result)
    cox_result[i,1] <- round(result$coefficients,2)
    cox_result[i,2] <- round(exp(result$coefficients),2)
    cox_result[i,3] <- round(1-pchisq(result$score, df=1),5)
  }
  cox_result_list[[l]]<-cox_result
  # cox_result_top <- subset(cox_result,coefficients>0)
  # cox_result_top<-top_n(cox_result_top,-10,pvalue)
  cox_result_top<-top_n(cox_result,-10,pvalue)
  indel_top <- indel[,rownames(cox_result_top)]
  
  risk.score <- setNames(rep(0,sample), rownames(indel))
  weight <- 0
  for (i in 1:sample) {
    for (j in 1:N)  {
      if (indel_top[i,j]=="Del" | indel_top[i,j]=="Amp") {weight <- weight + cox_result_top$coefficients[j]}
      if (indel_top[i,j]=="NC")                          {weight <- weight}
    }
    risk.score[i] <- weight
    weight <- 0
  }
  risk.score_list[[l]] <- risk.score
  ####FIGURE
  table(risk.score)
  survival.data <- data.frame(OSMonth=indel$OSMonth,Death=indel$Death,row.names = rownames(indel))
  survival.data$risk <- ifelse (risk.score > median(risk.score) , "high" , "low")
  fit <- survfit(Surv(OSMonth, Death) ~ risk, data=survival.data)
  png(file=paste("Top/",names(cohort)[l],"median risk score.png"))
  plot <- ggsurvplot(fit, pval=TRUE, xlab = "Time in months", linetype = "strata",
                     title=paste(cohortname[l],"OS - Top 8 genes risk score"),legend="bottom", legend.title = "Risk score", 
                     risk.table = T, tables.theme = theme_cleantable(), risk.table.fontsize=5,
                     font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14) 
  print(plot)
  dev.off()
}

cox_result_list[[1]]
cox_result_list[[3]]
####OS LOOCV####
coefficient_CV_list <-list()
pvalue_CV_list <-list()
for (l in 1:4){
  indel <- cohort_indel[[l]]
  sample <- nrow(indel)
  #coxph
  coefficient_CV <- matrix(0, sample,nGene)
  pvalue_CV <- matrix(0, sample,nGene)
  colnames(coefficient_CV)<-colnames(indel)[1:nGene]
  rownames(coefficient_CV)<-rownames(indel)
  colnames(pvalue_CV)<-colnames(indel)[1:nGene]
  rownames(pvalue_CV)<-rownames(indel)
  for (j in 1:sample) { 
    subdata <- indel[-j,] # LOO sample
    for (i in 1:nGene) {
      df <-data.frame(OSMonth=indel$OSMonth[-j],Death=indel$Death[-j],CNV=subdata[,i],row.names = rownames(subdata))
      if (i %in% amp & length(unique(subdata[,i]))==2)
      {
        df$CNV<- factor(df$CNV, levels = c("NC","Amp"))
        result_CV <- coxph(Surv(OSMonth, Death) ~ CNV, data=df)
      }
      if (i %in% del & length(unique(subdata[,i]))==2) 
      {
        df$CNV<- factor(df$CNV, levels = c("NC","Del"))
        result_CV <- coxph(Surv(OSMonth, Death) ~ CNV, data=df)
      }
      coefficient_CV[j,i] <- result_CV$coefficient
      pvalue_CV[j,i] <- 1-pchisq(result_CV$score, df=1)
    }
  } 
  coefficient_CV_list[[l]]<-coefficient_CV
  pvalue<-apply(pvalue_CV,2,mean) 
  pvalue_CV_list[[l]]<-pvalue
  pvalue<-pvalue[pvalue!=0] #remove those p value=0
  #risk score
  # coefficient<-apply(coefficient_CV,2,mean) 
  # coefficient<-coefficient[coefficient>0]
  # pvalue<-pvalue[names(coefficient)]
  top.gene <- sort(pvalue)[1:N]
  indel_top <- indel[,names(top.gene)]
  coefficient_CV_top <- coefficient_CV[,names(top.gene)]
  risk.score.CV <- setNames(rep(0,sample), rownames(indel_top))
  weight <- 0
  for (i in 1:sample) {
    for (j in 1:N)  {
      if (indel_top[i,j]=="Del" | indel_top[i,j]=="Amp")      {weight <- weight + coefficient_CV_top[i,j]}
      if (indel_top[i,j]=="NC")  {weight <- weight}
    }
    risk.score.CV[i] <- weight
    weight <- 0
  }
  ####FIGURE
  table(risk.score.CV)
  survival.data <- data.frame(OSMonth=indel$OSMonth,Death=indel$Death,row.names = rownames(indel))
  survival.data$risk <- ifelse (risk.score.CV > median(risk.score.CV) , "high" , "low")
  fit <- survfit(Surv(OSMonth, Death) ~ risk, data=survival.data)
  png(file=paste("Top/",names(cohort)[l],"LOOCV median risk score.png"))
  plot <- ggsurvplot(fit, pval=TRUE, xlab = "Time in months", linetype = "strata",
                     title=paste(cohortname[l],"OS - Top 8 genes LOOCV risk score"),legend="bottom", legend.title = "Risk score", 
                     risk.table = T, tables.theme = theme_cleantable(), risk.table.fontsize=5,
                     font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14) 
  print(plot)
  dev.off()
}
####Progression free survival####
for (l in 5:6){
  indel <- cohort_indel[[l-4]]
  sample <- nrow(indel)
  cox_result <-  data.frame(coefficients=rep(NA,nGene), HazardRatio=rep(NA,nGene), pvalue=rep(NA,nGene), row.names = colnames(indel[1:nGene]))
  for (i in 1:nGene) {
    df <-data.frame(PFSMonth=indel$PFSMonth,Progression=indel$Progression,CNV=indel[,i])
    if (i %in% amp) {df$CNV<- factor(df$CNV, levels = c("NC","Amp"))}
    if (i %in% del) {df$CNV<- factor(df$CNV, levels = c("NC","Del"))}
    result <- coxph(Surv(PFSMonth, Progression) ~ CNV, data=df)
    summary(result)
    cox_result[i,1] <- round(result$coefficients,2)
    cox_result[i,2] <- round(exp(result$coefficients),2)
    cox_result[i,3] <- round(1-pchisq(result$score, df=1),5)
  }
  cox_result_list[[l]]<-cox_result
  # cox_result_top <- subset(cox_result,coefficients>0)
  # cox_result_top<-top_n(cox_result_top,-10,pvalue)
  cox_result_top<-top_n(cox_result,-10,pvalue)
  indel_top <- indel[,rownames(cox_result_top)]
  risk.score <- setNames(rep(0,sample), rownames(indel))
  weight <- 0
  for (i in 1:sample) {
    for (j in 1:N)  {
      if (indel_top[i,j]=="Del" | indel_top[i,j]=="Amp") {weight <- weight + cox_result_top$coefficients[j]}
      if (indel_top[i,j]=="NC")                          {weight <- weight}
    }
    risk.score[i] <- weight
    weight <- 0
  }
  risk.score_list[[l]] <- risk.score
  ####FIGURE
  table(risk.score)
  survival.data <- data.frame(PFSMonth=indel$PFSMonth,Progression=indel$Progression,row.names = rownames(indel))
  survival.data$risk <- ifelse (risk.score > median(risk.score) , "high" , "low")
  fit <- survfit(Surv(PFSMonth, Progression) ~ risk, data=survival.data)
  png(file=paste("Top/",names(cohort)[l],"median risk score.png"))
  plot <- ggsurvplot(fit, pval=TRUE, xlab = "Time in months", linetype = "strata",
                     title=paste(cohortname[l],"PFS - Top 8 genes risk score"),legend="bottom", legend.title = "Risk score", 
                     risk.table = T, tables.theme = theme_cleantable(), risk.table.fontsize=5,
                     font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14) 
  print(plot)
  dev.off()
}
####PFS LOOCV####
for (l in 5:6){
  indel <- cohort_indel[[l-4]]
  sample <- nrow(indel)
  #coxph
  coefficient_CV <- matrix(0, sample,nGene)
  pvalue_CV <- matrix(0, sample,nGene)
  colnames(coefficient_CV)<-colnames(indel)[1:nGene]
  rownames(coefficient_CV)<-rownames(indel)
  colnames(pvalue_CV)<-colnames(indel)[1:nGene]
  rownames(pvalue_CV)<-rownames(indel)
  for (j in 1:sample) { 
    subdata <- indel[-j,] # LOO sample
    for (i in 1:nGene) {
      df <-data.frame(PFSMonth=indel$PFSMonth[-j],Progression=indel$Progression[-j],CNV=subdata[,i],row.names = rownames(subdata))
      if (i %in% amp & length(unique(subdata[,i]))==2)
      {
        df$CNV<- factor(df$CNV, levels = c("NC","Amp"))
        result_CV <- coxph(Surv(PFSMonth, Progression) ~ CNV, data=df)
      }
      if (i %in% del & length(unique(subdata[,i]))==2) 
      {
        df$CNV<- factor(df$CNV, levels = c("NC","Del"))
        result_CV <- coxph(Surv(PFSMonth, Progression) ~ CNV, data=df)
      }
      coefficient_CV[j,i] <- result_CV$coefficient
      pvalue_CV[j,i] <- 1-pchisq(result_CV$score, df=1)
    }
  } 
  coefficient_CV_list[[l]]<-coefficient_CV
  pvalue<-apply(pvalue_CV,2,mean) 
  pvalue_CV_list[[l]]<-pvalue
  pvalue<-pvalue[pvalue!=0] #remove those p value=0
  #risk score
  # coefficient<-apply(coefficient_CV,2,mean) 
  # coefficient<-coefficient[coefficient>0]
  # pvalue<-pvalue[names(coefficient)]
  top.gene <- sort(pvalue)[1:N]
  indel_top <- indel[,names(top.gene)]
  coefficient_CV_top <- coefficient_CV[,names(top.gene)]
  risk.score.CV <- setNames(rep(0,sample), rownames(indel_top))
  weight <- 0
  for (i in 1:sample) {
    for (j in 1:N)  {
      if (indel_top[i,j]=="Del" | indel_top[i,j]=="Amp")      {weight <- weight + coefficient_CV_top[i,j]}
      if (indel_top[i,j]=="NC")  {weight <- weight}
    }
    risk.score.CV[i] <- weight
    weight <- 0
  }
  ####FIGURE
  table(risk.score.CV)
  survival.data <- data.frame(PFSMonth=indel$PFSMonth,Progression=indel$Progression,row.names = rownames(indel))
  survival.data$risk <- ifelse (risk.score.CV > median(risk.score.CV) , "high" , "low")
  fit <- survfit(Surv(PFSMonth, Progression) ~ risk, data=survival.data)
  png(file=paste("Top/",names(cohort)[l],"LOOCV median risk score.png"))
  plot <- ggsurvplot(fit, pval=TRUE, xlab = "Time in months", linetype = "strata",
                     title=paste(cohortname[l],"PFS - Top 8 genes LOOCV risk score"),legend="bottom", legend.title = "Risk score", 
                     risk.table = T, tables.theme = theme_cleantable(), risk.table.fontsize=5,
                     font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14) 
  print(plot)
  dev.off()
}
####
names(cox_result_list)=names(cohort)
names(risk.score_list)=names(cohort)
names(coefficient_CV_list)<-names(cohort)
names(pvalue_CV_list)<-names(cohort)
pvalue_CV_list<-lapply(pvalue_CV_list,function(x) sort(x))
####Gene of interest####
GOI<-c(names(pvalue_CV_list[[1]])[1:5],names(pvalue_CV_list[[3]])[c(1:4)],names(pvalue_CV_list[[5]])[1:4])
GOI<-GOI[!duplicated(GOI)]
####GOI Overall Survival####
#####Multi-genes survival
cox_result_GOI_list <- list()
risk.score_GOI_list <- list()
for (l in 1:4){
  indel <- cohort_indel[[l]]
  sample <- nrow(indel)
  indel_GOI<-indel[,GOI]
  cox_result <-  data.frame(coefficients=rep(0,length(GOI)), pvalue=rep(NA,length(GOI)), row.names = colnames(indel_GOI))
  
  for (i in 1:length(GOI)) {
    df <-data.frame(OSMonth=indel$OSMonth,Death=indel$Death,CNV=indel_GOI[,i])
    if ("Amp" %in% df$CNV & length(unique(indel_GOI[,i]))==2) {
      df$CNV<- factor(df$CNV, levels = c("NC","Amp"))
      result <- coxph(Surv(OSMonth, Death) ~ CNV, data=df)
    }
    if ("Del" %in% df$CNV & length(unique(indel_GOI[,i]))==2) {
      df$CNV<- factor(df$CNV, levels = c("NC","Del"))
      result <- coxph(Surv(OSMonth, Death) ~ CNV, data=df)
    }
    cox_result[i,1] <- result$coefficients
    cox_result[i,2] <- 1-pchisq(result$score, df=1)
  }
  cox_result[is.na(cox_result)]<-0
  cox_result_GOI_list[[l]]<-cox_result
  
  risk.score <- setNames(rep(0,sample), rownames(indel_GOI))
  weight <- 0
  for (i in 1:sample) {
    for (j in 1:length(GOI))  {
      if (indel_GOI[i,j]=="Del" | indel_GOI[i,j]=="Amp") {weight <- weight + cox_result$coefficients[j]}
      if (indel_GOI[i,j]=="NC")                          {weight <- weight}
    }
    risk.score[i] <- weight
    weight <- 0
  }
  risk.score_GOI_list[[l]] <- risk.score
  ####FIGURE
  table(risk.score)
  survival.data <- data.frame(OSMonth=indel$OSMonth,Death=indel$Death,row.names = rownames(indel))
  survival.data$risk <- ifelse (risk.score > median(risk.score) , "high" , "low")
  fit <- survfit(Surv(OSMonth, Death) ~ risk, data=survival.data)
  table<-data.frame(summary(fit)$table)
  png(file=paste("GOI/",names(cohort)[l],"median risk score.png"))
  plot <- ggsurvplot(fit, pval=TRUE, xlab = "Time in months", linetype = "strata",
                     title=paste(cohortname[l],"OS - 11 genes panel risk score"),legend="bottom", legend.title = "Risk score", 
                     risk.table = T, tables.theme = theme_cleantable(), risk.table.fontsize=5,
                     font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14) 
  print(plot)
  dev.off()
}
####GOI OS LOOCV ####
coefficient_GOI_CV_list <-list()
risk.score_GOI_CV_list <- list()
for (l in 1:4){
  indel <- cohort_indel[[l]]
  sample <- nrow(indel)
  indel_GOI<-indel[,GOI]
  #coxph
  coefficient_CV <- matrix(0, sample,length(GOI))
  pvalue_CV <- matrix(0, sample,length(GOI))
  colnames(coefficient_CV)<-colnames(indel_GOI)
  rownames(coefficient_CV)<-rownames(indel_GOI)
  colnames(pvalue_CV)<-colnames(indel_GOI)
  rownames(pvalue_CV)<-rownames(indel_GOI)
  for (j in 1:sample) { 
    subdata <- indel_GOI[-j,] # LOO sample
    for (i in 1:length(GOI)) {
      df <-data.frame(OSMonth=indel$OSMonth[-j],Death=indel$Death[-j],CNV=subdata[,i],row.names = rownames(subdata))
      if ("Amp" %in% df$CNV & length(unique(subdata[,i]))==2)
      {
        df$CNV<- factor(df$CNV, levels = c("NC","Amp"))
        result_CV <- coxph(Surv(OSMonth, Death) ~ CNV, data=df)
      }
      if ("Del" %in% df$CNV & length(unique(subdata[,i]))==2) 
      {
        df$CNV<- factor(df$CNV, levels = c("NC","Del"))
        result_CV <- coxph(Surv(OSMonth, Death) ~ CNV, data=df)
      }
      coefficient_CV[j,i] <- result_CV$coefficient
      pvalue_CV[j,i] <- 1-pchisq(result_CV$score, df=1)
    }
  } 
  coefficient_GOI_CV_list[[l]]<-coefficient_CV
  #risk score
  risk.score.CV <- setNames(rep(0,sample), rownames(indel_GOI))
  weight <- 0
  for (i in 1:sample) {
    for (j in 1:length(GOI))  {
      if (indel_GOI[i,j]=="Del" | indel_GOI[i,j]=="Amp")      {weight <- weight + coefficient_CV[i,j]}
      if (indel_GOI[i,j]=="NC")  {weight <- weight}
    }
    risk.score.CV[i] <- weight
    weight <- 0
  }
  risk.score_GOI_CV_list[[l]]<-risk.score.CV
  ####FIGURE
  table(risk.score.CV)
  survival.data <- data.frame(OSMonth=indel$OSMonth,Death=indel$Death,row.names = rownames(indel))
  survival.data$risk <- ifelse (risk.score.CV > median(risk.score.CV) , "high" , "low")
  fit <- survfit(Surv(OSMonth, Death) ~ risk, data=survival.data)
  png(file=paste("GOI/",names(cohort)[l],"LOOCV median risk score.png"))
  plot <- ggsurvplot(fit, pval=TRUE, xlab = "Time in months", linetype = "strata",
                     title=paste(cohortname[l],"OS - 11 genes panel LOOCV risk score"),legend="bottom", legend.title = "Risk score", 
                     risk.table = T, tables.theme = theme_cleantable(), risk.table.fontsize=5,
                     font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14) 
  print(plot)
  dev.off()
}

####GOI Progression free survival####
for (l in 5:6){
  indel <- cohort_indel[[l-4]]
  sample <- nrow(indel)
  indel_GOI<-indel[,GOI]
  cox_result <-  data.frame(coefficients=rep(0,length(GOI)), pvalue=rep(NA,length(GOI)), row.names = colnames(indel_GOI))
  
  for (i in 1:length(GOI)) {
    df <-data.frame(PFSMonth=indel$PFSMonth,Progression=indel$Progression,CNV=indel_GOI[,i])
    if ("Amp" %in% df$CNV & length(unique(indel_GOI[,i]))==2) {
      df$CNV<- factor(df$CNV, levels = c("NC","Amp"))
      result <- coxph(Surv(PFSMonth, Progression) ~ CNV, data=df)
    }
    if ("Del" %in% df$CNV & length(unique(indel_GOI[,i]))==2) {
      df$CNV<- factor(df$CNV, levels = c("NC","Del"))
      result <- coxph(Surv(PFSMonth, Progression) ~ CNV, data=df)
    }
    cox_result[i,1] <- result$coefficients
    cox_result[i,2] <- 1-pchisq(result$score, df=1)
  }
  cox_result[is.na(cox_result)]<-0
  cox_result_GOI_list[[l]]<-cox_result
  
  risk.score <- setNames(rep(0,sample), rownames(indel_GOI))
  weight <- 0
  for (i in 1:sample) {
    for (j in 1:length(GOI))  {
      if (indel_GOI[i,j]=="Del" | indel_GOI[i,j]=="Amp") {weight <- weight + cox_result$coefficients[j]}
      if (indel_GOI[i,j]=="NC")                          {weight <- weight}
    }
    risk.score[i] <- weight
    weight <- 0
  }
  risk.score_GOI_list[[l]] <- risk.score
  ####FIGURE
  table(risk.score)
  survival.data <- data.frame(PFSMonth=indel$PFSMonth,Progression=indel$Progression,row.names = rownames(indel))
  survival.data$risk <- ifelse (risk.score > median(risk.score) , "high" , "low")
  fit <- survfit(Surv(PFSMonth, Progression) ~ risk, data=survival.data)
  table<-data.frame(summary(fit)$table)
  png(file=paste("GOI/",names(cohort)[l],"median risk score.png"))
  plot <- ggsurvplot(fit, pval=TRUE, xlab = "Time in months", linetype = "strata",
                     title=paste(cohortname[l],"PFS - 11 genes panel risk score"),legend="bottom", legend.title = "Risk score", 
                     risk.table = T, tables.theme = theme_cleantable(), risk.table.fontsize=5,
                     font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14) 
  print(plot)
  dev.off()
}
####GOI PFS LOOCV####
for (l in 5:6){
  indel <- cohort_indel[[l-4]]
  sample <- nrow(indel)
  indel_GOI<-indel[,GOI]
  #coxph
  coefficient_CV <- matrix(0, sample,length(GOI))
  pvalue_CV <- matrix(0, sample,length(GOI))
  colnames(coefficient_CV)<-colnames(indel_GOI)
  rownames(coefficient_CV)<-rownames(indel_GOI)
  colnames(pvalue_CV)<-colnames(indel_GOI)
  rownames(pvalue_CV)<-rownames(indel_GOI)
  for (j in 1:sample) { 
    subdata <- indel_GOI[-j,] # LOO sample
    for (i in 1:length(GOI)) {
      df <-data.frame(PFSMonth=indel$PFSMonth[-j],Progression=indel$Progression[-j],CNV=subdata[,i],row.names = rownames(subdata))
      if ("Amp" %in% df$CNV & length(unique(subdata[,i]))==2)
      {
        df$CNV<- factor(df$CNV, levels = c("NC","Amp"))
        result_CV <- coxph(Surv(PFSMonth, Progression) ~ CNV, data=df)
      }
      if ("Del" %in% df$CNV & length(unique(subdata[,i]))==2) 
      {
        df$CNV<- factor(df$CNV, levels = c("NC","Del"))
        result_CV <- coxph(Surv(PFSMonth, Progression) ~ CNV, data=df)
      }
      coefficient_CV[j,i] <- result_CV$coefficient
      pvalue_CV[j,i] <- 1-pchisq(result_CV$score, df=1)
    }
  } 
  coefficient_GOI_CV_list[[l]]<-coefficient_CV
  #risk score
  risk.score.CV <- setNames(rep(0,sample), rownames(indel_GOI))
  weight <- 0
  for (i in 1:sample) {
    for (j in 1:length(GOI))  {
      if (indel_GOI[i,j]=="Del" | indel_GOI[i,j]=="Amp")      {weight <- weight + coefficient_CV[i,j]}
      if (indel_GOI[i,j]=="NC")  {weight <- weight}
    }
    risk.score.CV[i] <- weight
    weight <- 0
  }
  risk.score_GOI_CV_list[[l]]<-risk.score.CV
  ####FIGURE
  table(risk.score.CV)
  survival.data <- data.frame(PFSMonth=indel$PFSMonth,Progression=indel$Progression,row.names = rownames(indel))
  survival.data$risk <- ifelse (risk.score.CV > median(risk.score.CV) , "high" , "low")
  fit <- survfit(Surv(PFSMonth, Progression) ~ risk, data=survival.data)
  png(file=paste("GOI/",names(cohort)[l],"LOOCV median risk score.png"))
  plot <- ggsurvplot(fit, pval=TRUE, xlab = "Time in months", linetype = "strata",
                     title=paste(cohortname[l],"PFS - 11 genes panel LOOCV risk score"),legend="bottom", legend.title = "Risk score", 
                     risk.table = T, tables.theme = theme_cleantable(), risk.table.fontsize=5,
                     font.title=c(16, "bold"),font.x=c(14, "bold"), font.y=c(14, "bold"),font.legend=14,font.tickslab=14) 
  print(plot)
  dev.off()
}

names(cox_result_GOI_list)=names(cohort)
names(risk.score_GOI_list)=names(cohort) 
names(coefficient_GOI_CV_list)<-names(cohort)
names(risk.score_GOI_CV_list)<-names(cohort)
####risk score####
risk.score_all <- unlist(risk.score_GOI_CV_list[c(1,3,5)])
risk.score_df <- data.frame(cohort=substr(names(risk.score_all),1,5),
                            patient=substr(names(risk.score_all),6,nchar(names(risk.score_all))),
                            score=round(risk.score_all,3))
risk.score_df$patient <- gsub("S.","",risk.score_df$patient)
risk.score_df$patient <- gsub("\\.","",risk.score_df$patient)
write.table(risk.score_df,"risk.score_df.txt",col.names = NA, sep="\t",quote = F)
#####risk score histogram####
data<-subset(risk.score_df,cohort=="pv1OS")
hist(data$score, main="PROMOTE pretreatment OS risk score", xlab="Risk score")
summary(data)
data<-subset(risk.score_df,cohort=="pv1PF")
hist(data$score, main="PROMOTE pretreatment PFS risk score", xlab="Risk score")
summary(data)
data<-subset(risk.score_df,cohort=="sv1OS")
hist(data$score, main="HR pretreatment OS risk score", xlab="Risk score",xlim = c(-0.5,1.6))
summary(data)
#####risk score boxplot####
data1<-subset(risk.score_df,cohort=="pv1OS")
data1<-merge(data1,promote_visit1,by.x = "patient", by.y = "PatientID")
data1$risk <-ifelse(data1$score>median(data1$score),"OS - High risk","OS - Low risk")
High<-subset(data1,risk=="OS - High risk")$OSMonth
Low<-subset(data1,risk=="OS - Low risk")$OSMonth
median(High)-median(Low)
ttest_p <- t.test(High, Low)
pvalue_p1 <- round(ttest_p$p.value,5)

data3<-subset(risk.score_df,cohort=="pv1PF")
data3<-merge(data3,promote_visit1,by.x = "patient", by.y = "PatientID")
data3$risk <-ifelse(data3$score>median(data3$score),"PFS - High risk","PFS - Low risk")
hist(data3$PFSMonth)
High<-subset(data3,risk=="PFS - High risk")$PFSMonth
Low<-subset(data3,risk=="PFS - Low risk")$PFSMonth
median(High)-median(Low)
ttest_p <- t.test(High, Low)
pvalue_p2 <- round(ttest_p$p.value,5)

colnames(data1)[9]="Months"
colnames(data3)[10]="Months"
data<-rbind(data1[,c(1,2,3,9,11)],data3[,c(1,2,3,10,11)])
png(file="PROMOTE OSPFS time and risk score.png")
plot <- ggplot(data, aes(x=risk, y=Months)) + geom_boxplot() +geom_jitter(shape=16, position=position_jitter(0.2))+ 
  labs(title="PROMOTE OS/PFS time and risk score",x="", y = "Months")+ theme_classic()+
  annotate("text", x = 1.5, y = 60, label = paste("p=",pvalue_p1,sep = ""),size=5)+
  annotate("text", x = 3.5, y = 60, label = paste("p=",pvalue_p2,sep = ""),size=5)
print(plot)
dev.off()

