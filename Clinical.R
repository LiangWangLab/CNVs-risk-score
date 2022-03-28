library(devtools)
library(ComplexHeatmap)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library("ggpmisc")
##############
promote <- read.table("promote_final.txt", header = T, sep = "\t")
promote$Response <- gsub("0","Responder", promote$Response)
promote$Response <- gsub("1","Non-responder", promote$Response)
####PSA####
PSA <- promote[,c("Patient","Response","Visit.1.PSA","Visit.2.PSA")]
PSA$Delta <-(PSA$Visit.2.PSA - PSA$Visit.1.PSA) / PSA$Visit.1.PSA *100
# PSA<-PSA[order(PSA$Delta, decreasing = TRUE),]
# PSA <- PSA %>% arrange(desc(Delta))
a=nrow(subset(PSA,Response=="Responder" & Delta<0))
b=nrow(subset(PSA,Response=="Non-responder" & Delta<0))
c=nrow(subset(PSA,Response=="Responder" & Delta>0))
d=nrow(subset(PSA,Response=="Non-responder" & Delta>0))
annotation <- matrix(c('PSA/Response','Responder',  'Non-responder',
                   'PSA-',     paste(a,"(",round(a/(a+b)*100,2),"%)"),            paste(b,"(",round(b/(a+b)*100,2),"%)"),         
                   'PSA+',     paste(c,"(",round(c/(c+d)*100,2),"%)"),            paste(d,"(",round(d/(c+d)*100,2),"%)")), 
                 ncol=3, byrow=TRUE)
colnames(annotation) <- annotation[1,];annotation <- annotation[-1, ] 
annotation <- as.data.frame(annotation)
test_df<-data.frame(
  "Responder" = c(a, c),
  "Non-responder" = c(b, d),
  row.names = c("PSA-", "PSA+"))
test <- fisher.test(test_df)
test
png(file="PROMOTE - PSA - Treatment response.png")
ggplot(data=PSA, aes(x=reorder(Patient, -Delta), y=Delta, fill=Response))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title="Prostate-Specific Antigen and Treatment response",x="Patients", y = "Change from baseline (%)")+ylim(-100,600)+
  theme(text = element_text(size=16),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x =element_blank())+
  annotate(geom = "table",x = 85,y = 150,label = list(annotation),size = 4)+
  annotate(geom = "text",x = 60,y = 100,label = paste("p =",format.pval(test$p.value)),size = 4)
dev.off()

####Chromogranin A (CgA)####
CGA<-promote[,c("Patient","Response","CGA.Visit.1","CGA.Visit.2")]
CGA$Delta <-(CGA$CGA.Visit.2 - CGA$CGA.Visit.1) / CGA$CGA.Visit.1 *100
a=nrow(subset(CGA,Response=="Responder" & Delta<0))
b=nrow(subset(CGA,Response=="Non-responder" & Delta<0))
c=nrow(subset(CGA,Response=="Responder" & Delta>0))
d=nrow(subset(CGA,Response=="Non-responder" & Delta>0))
annotation <- matrix(c('CGA/Response','Responder',  'Non-responder',
                       'CGA-',     paste(a,"(",round(a/(a+b)*100,2),"%)"),            paste(b,"(",round(b/(a+b)*100,2),"%)"),
                       'CGA+',     paste(c,"(",round(c/(c+d)*100,2),"%)"),            paste(d,"(",round(d/(c+d)*100,2),"%)")),
                     ncol=3, byrow=TRUE)
colnames(annotation) <- annotation[1,];annotation <- annotation[-1, ]
annotation <- as.data.frame(annotation)
test_df<-data.frame(
  "Responder" = c(a, c),
  "Non-responder" = c(b, d),
  row.names = c("CGA-", "CGA+"))
test <- fisher.test(test_df)
test
png(file="PROMOTE - CGA - Treatment response.png")
ggplot(data=CGA, aes(x=reorder(Patient, -Delta), y=Delta, fill=Response))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title="Chromogranin A and Treatment response",x="Patients", y = "Change from baseline (%)")+ylim(-100,150)+
  theme(text = element_text(size=16),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x =element_blank())+
  annotate(geom = "table",x = 90,y = 120,label = list(annotation),size = 4)+
  annotate(geom = "text",x = 60,y = 80,label = paste("p =",format.pval(test$p.value)),size = 4)
dev.off()
####LDH and Treatment response####
LDH <- promote[,c("Patient","Response","Serum.LDH.Visit.1","Serum.LDH.Visit.2")]
LDH$Delta <-(LDH$Serum.LDH.Visit.2 - LDH$Serum.LDH.Visit.1) / LDH$Serum.LDH.Visit.1 *100
a=nrow(subset(LDH,Response=="Responder" & Delta<0))
b=nrow(subset(LDH,Response=="Non-responder" & Delta<0))
c=nrow(subset(LDH,Response=="Responder" & Delta>0))
d=nrow(subset(LDH,Response=="Non-responder" & Delta>0))
annotation <- matrix(c('LDH/Response','Responder',  'Non-responder',
                       'LDH-',     paste(a,"(",round(a/(a+b)*100,2),"%)"),            paste(b,"(",round(b/(a+b)*100,2),"%)"),
                       'LDH+',     paste(c,"(",round(c/(c+d)*100,2),"%)"),            paste(d,"(",round(d/(c+d)*100,2),"%)")),
                     ncol=3, byrow=TRUE)
colnames(annotation) <- annotation[1,];annotation <- annotation[-1, ]
annotation <- as.data.frame(annotation)
test_df<-data.frame(
  "Responder" = c(a, c),
  "Non-responder" = c(b, d),
  row.names = c("LDH-", "LDH+"))
test <- fisher.test(test_df)
test
png(file="PROMOTE - LDH - Treatment response.png")
ggplot(data=LDH, aes(x=reorder(Patient, -Delta), y=Delta, fill=Response))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title="Lactate dehydrogenase and Treatment response",x="Patients", y = "Change from baseline (%)")+ylim(-50,130)+
  theme(text = element_text(size=16),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x =element_blank())+
  annotate(geom = "table",x = 80,y = 60,label = list(annotation),size = 4)+
  annotate(geom = "text",x = 60,y = 30,label = paste("p =",format.pval(test$p.value)),size = 4)
dev.off()

####PSA and risk score####
riskscore <- read.table("risk.score_df.txt", header = T, sep = "\t")
riskscore <- subset(riskscore, cohort=="pv1")
riskscore <- riskscore[riskscore$patient %in% promote$v1ID,]
riskscore$Risk <- ifelse(riskscore$socre > median(riskscore$socre),"High","Low")
PSArs <- merge(PSA,promote[,1:3],by = "Patient")
PSArs <- merge(PSArs,riskscore, by.x = "v1ID", by.y = "patient")

a=nrow(subset(PSArs,Risk=="Low" & Delta<0))
b=nrow(subset(PSArs,Risk=="High" & Delta<0))
c=nrow(subset(PSArs,Risk=="Low" & Delta>0))
d=nrow(subset(PSArs,Risk=="High" & Delta>0))
annotation <- matrix(c('PSA/Risk',  'Low',                                             'High',
                       'PSA-',     paste(a,"(",round(a/(a+b)*100,2),"%)"),            paste(b,"(",round(b/(a+b)*100,2),"%)"),         
                       'PSA+',     paste(c,"(",round(c/(c+d)*100,2),"%)"),            paste(d,"(",round(d/(c+d)*100,2),"%)")), 
                     ncol=3, byrow=TRUE)
colnames(annotation) <- annotation[1,];annotation <- annotation[-1, ] 
annotation <- as.data.frame(annotation)
test_df<-data.frame(
  "Low" = c(a, c),
  "High" = c(b, d),
  row.names = c("PSA-", "PSA+"))
test <- fisher.test(test_df)
test
png(file="PROMOTE - PSA - Risk score.png")
ggplot(data=PSArs, aes(x=reorder(Patient, -Delta), y=Delta, fill=Risk))+geom_bar(stat="identity",colour="black")+theme_classic()+
  labs(title="Prostate-Specific Antigen and Risk score",x="Patients", y = "Change from baseline (%)")+ylim(-100,600)+
  theme(text = element_text(size=16),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x =element_blank())+
  annotate(geom = "table",x = 80,y = 150,label = list(annotation),size = 4)+
  annotate(geom = "text",x = 60,y = 110,label = paste("p =",format.pval(test$p.value)),size = 5)+
  annotate(geom = "text",x = 60,y = 70,label = paste("Odd ratio =",round(test$estimate,2)),size = 5)
dev.off()
bar <- data.frame(t(test_df))
bar$Risk<-rownames(bar)
png(file="PROMOTE - PSA-Risk score-odd ratio.png")
ggplot(data=bar, aes(x=Risk, y=PSA..1)) +geom_bar(stat="identity",width =0.5) + labs(y="# of cases", title = "PSA+")+
  theme(text = element_text(size=16),axis.text = element_text(size=16))
dev.off()
####PSA and Single gene####
Indel <- read.table("pv1_indel.txt", header = T, sep = "\t")
PSA_Indel <- merge(PSArs,Indel, by.x = "v1ID", by.y = "X")
colnames(PSA_Indel)
for (i in 12:35){
  png(file=paste("PROMOTE - PSA -", colnames(PSA_Indel)[i],".png"))
  plot<-ggplot(data=PSA_Indel, aes(x=reorder(Patient, -Delta), y=Delta, fill=PSA_Indel[,i]))+geom_bar(stat="identity",colour="black")+theme_classic()+
    labs(title=paste("Prostate-Specific Antigen and", colnames(PSA_Indel)[i]),x="Patients", y = "Change from baseline (%)")+ylim(-100,600)+
    theme(text = element_text(size=16),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x =element_blank())+
    guides(fill=guide_legend(title=paste(colnames(PSA_Indel)[i])))
  print(plot)
  dev.off()
}
