library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library("ggpmisc")

Response <- read.table("promote_visit1_visit2.txt", header = T)
Response <- Response[27]#Response
pv1 <- read.table("promote_visit1_before_normalization.txt", header = T)
pv2 <- read.table("promote_visit2_before_normalization.txt", header = T)
coldata <- read.csv("colData_two_visit.txt", header = T, sep = "\t")
coldata <- coldata[,c(1,10)]
pv1 <- merge(pv1,coldata,by.x = 0, by.y = "PatientID")
row.names(pv1)=pv1$Patient
pv2 <- merge(pv2,coldata,by.x = 0, by.y = "PatientID")
row.names(pv2)=pv2$Patient
####Promote####
for (i in 2:25){
  ttest <- t.test(pv1[,i], pv2[,i], paired = TRUE, var.equal = TRUE)
  pvalue <- round(ttest$p.value,5)
  df <- rbind(pv1=pv1[i],pv2=pv2[i])
  df$ID <- rownames(df)
  df <- separate(df, "ID", c("Time", "ID"))
  df$Time <- gsub("pv","Visit ", df$Time)
  png(file=paste("promote/",colnames(pv1)[i],".png",sep = ""))
  plot <- ggplot(df, aes_string(x="Time", y=colnames(pv1)[i])) + geom_point() + geom_line(aes(group = ID)) + theme_classic() + labs(title=paste(colnames(pv1)[i]), x = "", y = "log2(copy number)") +
    theme(axis.text = element_text(colour = "black",size = 20),axis.title = element_text(size = 20),plot.title = element_text(size = 25,hjust = 0.5,face = "bold"),panel.border = element_rect(colour = "black", fill=NA, size=1.2)) +
    annotate("text", x = 1.5, y = max(df[,1])+0.1, label = paste("p=",pvalue,sep = ""),size=8)
  print(plot)
  dev.off()
}
####PROMOTE and response####
for (i in 2:25){
  df <- rbind(pv1=pv1[i],pv2=pv2[i])
  df$ID <- rownames(df)
  df <- separate(df, "ID", c("Time", "ID"))
  df$Time <- gsub("pv","Visit ", df$Time)
  df <- merge(df, Response, by.x = "ID", by.y = 0)
  df$Response <- gsub("0","Responder",df$Response)
  df$Response <- gsub("1","Nonresponder",df$Response)
  response1<-subset(df, Response=="Responder" & Time=="Visit 1")[,2]
  response2<-subset(df, Response=="Responder" & Time=="Visit 2")[,2]
  noresponse1<-subset(df, Response=="Nonresponder" & Time=="Visit 1")[,2]
  noresponse2<-subset(df, Response=="Nonresponder" & Time=="Visit 2")[,2]
  ttest0 <- t.test(response1, response2, paired = TRUE, var.equal = TRUE)
  pvalue0 <- round(ttest0$p.value,5)
  ttest1 <- t.test(noresponse1, noresponse2, paired = TRUE, var.equal = TRUE)
  pvalue1 <- round(ttest1$p.value,5)
  pvalue <- data.frame(label = c(paste("p=",pvalue0,sep = ""), paste("p=",pvalue1,sep = "")), 
                       Response=c("Responder","Nonresponder"), x=c(1.5,1.5), y=c(max(df[,2])+0.1, max(df[,2])+0.1))
  png(file=paste("promote_response/",colnames(pv1)[i],".png",sep = ""))
  plot <- ggplot(df, aes_string(x="Time", y=colnames(pv1)[i])) + geom_point() + geom_line(aes(group = ID)) + facet_grid(~Response) + 
    theme_classic() + labs(title=paste(colnames(pv1)[i]), x = "", y = "log2(copy number)") + geom_text(data = pvalue, mapping = aes(x = x, y = y, label = label), size =8)+
    theme(axis.text = element_text(colour = "black",size = 20),axis.title = element_text(size = 20),plot.title = element_text(size = 25,hjust = 0.5,face = "bold"),
          panel.border = element_rect(colour = "black", fill=NA, size=1.2),strip.text.x = element_text(size = 20)) 
  print(plot)
  dev.off()
}
####SDMS####
# sv1 <- read.table("SDMS_sharedpatients_visit1.txt", header = T)
# sv2 <- read.table("SDMS_sharedpatients_visit2.txt", header = T)

# for (i in 1:24){
#   ttest <- t.test(sv1[,i], sv2[,i], paired = TRUE, var.equal = TRUE)
#   pvalue <- round(ttest$p.value,5)
#   df <- rbind(sv1=sv1[i],sv2=sv2[i])
#   df$ID <- rownames(df)
#   df <- separate(df, "ID", c("Time", "ID"))
#   df$Time <- gsub("sv","Visit ", df$Time)
#   png(file=paste("SDMS/",colnames(sv1)[i],".png",sep = ""))
#   plot <- ggplot(df, aes_string(x="Time", y=colnames(sv1)[i])) + geom_point() + geom_line(aes(group = ID)) + theme_classic() + labs(title=paste(colnames(sv1)[i]), x = "", y = "Copy number ratio") +
#     theme(axis.text = element_text(colour = "black",size = 20),axis.title = element_text(size = 20),plot.title = element_text(size = 25,hjust = 0.5,face = "bold"),panel.border = element_rect(colour = "black", fill=NA, size=1.2)) +  
#     annotate("text", x = 1.5, y = max(df[,1])+0.1, label = paste("p=",pvalue,sep = ""),size=8) 
#   print(plot)
#   dev.off()
# }

