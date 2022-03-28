library(tidyverse)
library(data.table)
library(DESeq2)
####Promote####
Raw_promote<-read.csv("Promote_gene_counts.txt",header=T,sep = ",")
Raw_count_promote<-Raw_promote[,-(2:6)]
rownames(Raw_count_promote)<-Raw_count_promote[,1];Raw_count_promote<-Raw_count_promote[,-1]
colnames(Raw_count_promote)<-gsub("bams.","",colnames(Raw_count_promote))
colnames(Raw_count_promote)<-gsub(".bam","",colnames(Raw_count_promote))
Clinical_promote<-read.table("Promote clinical data.txt",header = T,sep = "\t")
colnames(Clinical_promote)[1]<-"PatientID"
Clinical_promote$PatientID<-gsub(" ","",Clinical_promote$PatientID)
Raw_count_promote<-Raw_count_promote[,colnames(Raw_count_promote) %in% Clinical_promote$PatientID]
#filter out gDNA
Clinical_promote<-subset(Clinical_promote,Type!="PromoteG")
Raw_count_promote<-Raw_count_promote[,colnames(Raw_count_promote) %in% Clinical_promote$PatientID]
#filter by mappable reads
Mappable_promote <- read.table("Promote_Read_summary.txt", sep="\t",header = T)
Clinical_promote<-merge(Clinical_promote,Mappable_promote,by.x="PatientID",by.y="Sample")
Clinical_promote<-subset(Clinical_promote,Mapped>2000000)
Raw_count_promote<-Raw_count_promote[,colnames(Raw_count_promote) %in% Clinical_promote$PatientID]
write.table(Clinical_promote[,1:11], file="Clinical_promote.txt", sep="\t", quote=F, row.names = F)
Clinical_promote$OStime<-Clinical_promote$OStime*30
Clinical_promote$PFStime<-Clinical_promote$PFStime*30
####SDMS####
Raw_SDMS<-read.csv("SDMS_gene_counts.txt",header=T,sep = "\t")
Raw_count_SDMS<-Raw_SDMS[,-(2:6)]
rownames(Raw_count_SDMS)<-Raw_count_SDMS[,1];Raw_count_SDMS<-Raw_count_SDMS[,-1]
colnames(Raw_count_SDMS)<-gsub(".bam","",colnames(Raw_count_SDMS))
colnames(Raw_count_SDMS) <- gsub("\\_.*", "", colnames(Raw_count_SDMS))
column_names <- unique(colnames(Raw_count_SDMS))
Raw_count_SDMS <- sapply(column_names, function(x) rowSums(Raw_count_SDMS[names(Raw_count_SDMS) %in% x])) %>% data.frame(check.names = F)
#sample with clinical information
Clinical_SDMS1<-read.table("SDMS clinical data.txt",header = T,sep = "\t")
Clinical_SDMS1$Patient.ID<-paste("X",Clinical_SDMS1$Patient.ID,sep = "")
#clinical information for 2nd-visit patients
Clinical_SDMS1$substring <- substr(Clinical_SDMS1$Patient.ID,3,5)
Clinical_SDMS3 <- data.frame(Patient.ID=colnames(Raw_count_SDMS)[121:162])
Clinical_SDMS3$substring <- substr(Clinical_SDMS3$Patient.ID,3,5)
Clinical_SDMS3<- merge(Clinical_SDMS1,Clinical_SDMS3,by="substring")
Clinical_SDMS4<- data.frame(PatientID=c(Clinical_SDMS3$Patient.ID.x,Clinical_SDMS3$Patient.ID.y),Death=rep(Clinical_SDMS3$X1..death..0.alive,2),OS_day=rep(Clinical_SDMS3$OS_day,2))
#data sample without clinical information
Clinical_SDMS2 <-data.frame(PatientID=subset(colnames(Raw_count_SDMS), !(colnames(Raw_count_SDMS) %in% Clinical_SDMS1$Patient.ID)))
Clinical_SDMS2<-subset(Clinical_SDMS2, !(Clinical_SDMS2$PatientID %in% Clinical_SDMS3$Patient.ID.y))
#clinical information final
Clinical_SDMS1<-subset(Clinical_SDMS1, !(Clinical_SDMS1$Patient.ID %in% Clinical_SDMS3$Patient.ID.x))
Clinical_SDMS <- data.frame(PatientID=c(Clinical_SDMS1$Patient.ID,Clinical_SDMS4$PatientID,Clinical_SDMS2$PatientID),
                             Death=c(Clinical_SDMS1$X1..death..0.alive,Clinical_SDMS4$Death,rep("NA",nrow(Clinical_SDMS2))),
                             OS_day=c(Clinical_SDMS1$OS_day,Clinical_SDMS4$OS_day,rep("NA",nrow(Clinical_SDMS2))))
Clinical_SDMS$Group <- "SDMS"
Clinical_SDMS$Visit <- ifelse(grepl("X1", Clinical_SDMS$PatientID), "1", "2")
Clinical_SDMS<-subset(Clinical_SDMS, Clinical_SDMS$PatientID %in% colnames(Raw_count_SDMS))
Clinical_SDMS$Type<-paste(Clinical_SDMS$Group,Clinical_SDMS$Visit)
#filter out NA clinical information 162 to 148
Clinical_SDMS<-subset(Clinical_SDMS,Death!="NA")
Raw_count_SDMS<-Raw_count_SDMS[,colnames(Raw_count_SDMS) %in% Clinical_SDMS$PatientID]
#filter by mappable reads 148 to 131
Mappable_SDMS <- read.table("SDMS_Read_summary.txt", sep="\t",header = T)
Mappable_SDMS$Sample<-paste("X",Mappable_SDMS$Sample,sep = "")
Clinical_SDMS<-merge(Clinical_SDMS,Mappable_SDMS,by.x="PatientID",by.y="Sample")
Clinical_SDMS<-subset(Clinical_SDMS,Mapped>2000000)
Raw_count_SDMS<-Raw_count_SDMS[,colnames(Raw_count_SDMS) %in% Clinical_SDMS$PatientID]
write.table(Clinical_SDMS[,1:6], file="Clinical_SDMS.txt", sep="\t", quote=F, row.names = F)
####Control####
Raw_Control <- read.csv("control_gene_counts.txt",header=T,sep = "\t")
Raw_count_Control<-Raw_Control[,-(2:6)]
rownames(Raw_count_Control)<-Raw_count_Control[,1];Raw_count_Control<-Raw_count_Control[,-1]
colnames(Raw_count_Control)<-gsub(".bam","",colnames(Raw_count_Control))
Raw_count_Control<-subset(Raw_count_Control,select=-(CPD_23))
Clinical_Control <- data.frame(PatientID=colnames(Raw_count_Control),Type="Control",Group="Control")
####Normalization####
colData <- data.frame(PatientID=c(Clinical_promote$PatientID,Clinical_SDMS$PatientID,Clinical_Control$PatientID),
                      Type=c(Clinical_promote$Type,Clinical_SDMS$Type,Clinical_Control$Type),
                      Group=c(Clinical_promote$Group,Clinical_SDMS$Group,Clinical_Control$Group),
                      Visit=c(Clinical_promote$Visit,Clinical_SDMS$Visit,rep("NA",15)),
                      Death=c(Clinical_promote$dead,Clinical_SDMS$Death,rep("NA",15)),
                      OSDay=c(Clinical_promote$OStime,Clinical_SDMS$OS_day,rep("NA",15)),
                      Response=c(Clinical_promote$Response,rep("NA",131),rep("NA",15)),
                      PFSDay=c(Clinical_promote$PFStime,rep("NA",131),rep("NA",15)),
                      Progression=c(Clinical_promote$Progression,rep("NA",131),rep("NA",15)),
                      Patient=c(Clinical_promote$EXN,substr(Clinical_SDMS$PatientID,3,5),rep("NA",15)))
write.table(colData, file="colData.txt", sep="\t", quote=F, row.names = F)
#colData$Type <-as.factor(colData$Type)
colData$Group <-as.factor(colData$Group)
Raw_count_All <-cbind(Raw_count_promote,Raw_count_SDMS,Raw_count_Control)
write.table(Raw_count_All, file="Raw_count_All.txt", sep="\t", quote=F, col.names = NA)
Raw_count_All_log <-log(Raw_count_All+1)
Raw_count_All_log <-round(log2(Raw_count_All+1),2)
write.table(Raw_count_All_log, file="Raw_count_All_log.txt", sep="\t", quote=F, col.names = NA)
Raw_count_All_median <-data.frame(apply(Raw_count_All,1,median))

dds <- DESeqDataSetFromMatrix(countData = Raw_count_All, colData = colData, design = ~ Group)
dds <- estimateSizeFactors(dds)
normalized_counts <- round(counts(dds, normalized=TRUE),2)
write.table(normalized_counts, file="Normalized_count_All.txt", sep="\t", quote=F, col.names = NA)
#genes of interest----
normalized_counts <- data.frame(t(normalized_counts))
GOI<-read.table("Genes of interest.txt")
normalized_counts <- normalized_counts[,GOI$V1]
write.table(normalized_counts, file="Normalized_count_24.txt", sep="\t", quote=F, col.names = NA)
normalized_counts_log <- round(log2(normalized_counts+1),2)
write.table(normalized_counts_log, file="Normalized_count_24_log.txt", sep="\t", quote=F, col.names = NA)

normalized_counts_control <- normalized_counts_log[grepl("CPD",rownames(normalized_counts_log)),]
write.table(normalized_counts_control, file="control_before_normalization.txt", sep="\t", quote=F, col.names = NA)
normalized_counts_control_median <- apply(normalized_counts_control,2,median)
normalized_counts_patients<-normalized_counts_log[!grepl("CPD",rownames(normalized_counts_log)),]

#before normalization
promote_v1_before <- normalized_counts_patients[promote_visit1$PatientID,]
write.table(promote_v1_before, file="promote_visit1_before_normalization.txt", sep="\t", quote=F, col.names = NA)
promote_v2_before <- normalized_counts_patients[promote_visit2$PatientID,]
write.table(promote_v2_before, file="promote_visit2_before_normalization.txt", sep="\t", quote=F, col.names = NA)
SDMS_v1_before <- normalized_counts_patients[SDMS_visit1$PatientID,]
write.table(SDMS_v1_before, file="SDMS_visit1_before_normalization.txt", sep="\t", quote=F, col.names = NA)
SDMS_v2_before <- normalized_counts_patients[SDMS_visit2$PatientID,]
write.table(SDMS_v2_before, file="SDMS_visit2_before_normalization.txt", sep="\t", quote=F, col.names = NA)

#normalize to control----
normalized_counts_patients_normalized <- normalized_counts_patients - normalized_counts_control_median

#patients with two visit
colData_two_visit<-colData[colData$Patient %in% colData$Patient[duplicated(colData$Patient)],]
colData_two_visit<-subset(colData_two_visit,Patient!="NA")
write.table(colData_two_visit, file="colData_two_visit.txt", sep="\t", quote=F, row.names = F)

####promote visit 1####
promote_visit1 <- subset(colData,Type=="Promote1")
promote_visit1_meta <- promote_visit1[,c(1,5:9)]
write.table(promote_visit1_meta, file="promote_visit1_meta.txt", sep="\t", quote=F, row.names = F)
promote_visit1_count <- normalized_counts_patients_normalized[promote_visit1$PatientID,]
write.table(promote_visit1_count, file="promote_visit1_count.txt", sep="\t", quote=F, col.names = NA)
####promote visit 2####
promote_visit2 <- subset(colData,Type=="Promote2")
promote_visit2_meta <- promote_visit2[,c(1,5:9)]
write.table(promote_visit2_meta, file="promote_visit2_meta.txt", sep="\t", quote=F, row.names = F)
promote_visit2_count <- normalized_counts_patients_normalized[promote_visit2$PatientID,]
write.table(promote_visit2_count, file="promote_visit2_count.txt", sep="\t", quote=F, col.names = NA)
####promote visit2 - visit1####
colData_promote <- subset(colData_two_visit,Group=="Promote")
promote_visit1_visit2_count <- normalized_counts_patients_normalized[colData_promote$PatientID,]
promote_visit1_visit2_gene <- merge(promote_visit1_visit2_count,colData_promote[,c(1,4,10)],by.x=0,by.y="PatientID")
promote_visit1_visit2_gene1 <- subset(promote_visit1_visit2_gene,Visit==1)
write.table(promote_visit1_visit2_gene1, file="promote_sharedpatients_visit1.txt", sep="\t", quote=F, row.names = F)
rownames(promote_visit1_visit2_gene1) <- promote_visit1_visit2_gene1[,"Patient"]
promote_visit1_visit2_gene1 <- promote_visit1_visit2_gene1[,2:25]
promote_visit1_visit2_gene2 <- subset(promote_visit1_visit2_gene,Visit==2)
write.table(promote_visit1_visit2_gene2, file="promote_sharedpatients_visit2.txt", sep="\t", quote=F, row.names = F)
rownames(promote_visit1_visit2_gene2) <- promote_visit1_visit2_gene2[,"Patient"]
promote_visit1_visit2_gene2 <- promote_visit1_visit2_gene2[,2:25]
promote_visit1_visit2_gene <- promote_visit1_visit2_gene2 - promote_visit1_visit2_gene1
colData_promote_unique <- colData_promote[,5:10]
colData_promote_unique <- colData_promote_unique[!duplicated(colData_promote_unique$Patient), ]
promote_visit1_visit2 <- merge(promote_visit1_visit2_gene,colData_promote_unique,by.x=0,by.y="Patient")
rownames(promote_visit1_visit2) <- promote_visit1_visit2[,1];promote_visit1_visit2 <- promote_visit1_visit2[,-1]
write.table(promote_visit1_visit2, file="promote_visit1_visit2.txt", sep="\t", quote=F, col.names = NA)
#write.table(promote_visit1_visit2_gene1, file="promote_sharedpatients_visit1.txt", sep="\t", quote=F, col.names = NA)
#write.table(promote_visit1_visit2_gene2, file="promote_sharedpatients_visit2.txt", sep="\t", quote=F, col.names = NA)
####SDMS visit 1####
SDMS_visit1 <- subset(colData,Type=="SDMS 1")
SDMS_visit1_meta <- SDMS_visit1[,c(1,5:6)]
write.table(SDMS_visit1_meta, file="SDMS_visit1_meta.txt", sep="\t", quote=F, row.names = F)
SDMS_visit1_count <- normalized_counts_patients_normalized[SDMS_visit1$PatientID,]
write.table(SDMS_visit1_count, file="SDMS_visit1_count.txt", sep="\t", quote=F, col.names = NA)
####SDMS visit 2####
SDMS_visit2 <- subset(colData,Type=="SDMS 2")
SDMS_visit2_meta <- SDMS_visit2[,c(1,5:6)]
write.table(SDMS_visit2_meta, file="SDMS_visit2_meta.txt", sep="\t", quote=F, row.names = F)
SDMS_visit2_count <- normalized_counts_patients_normalized[SDMS_visit2$PatientID,]
write.table(SDMS_visit2_count, file="SDMS_visit2_count.txt", sep="\t", quote=F, col.names = NA)
####SDMS visit2 - visit1####
colData_SDMS <- subset(colData_two_visit,Group=="SDMS")
SDMS_visit1_visit2_count <- normalized_counts_patients_normalized[colData_SDMS$PatientID,]
SDMS_visit1_visit2_gene <- merge(SDMS_visit1_visit2_count,colData_SDMS[,c(1,4,10)],by.x=0,by.y="PatientID")
SDMS_visit1_visit2_gene1 <- subset(SDMS_visit1_visit2_gene,Visit==1)
rownames(SDMS_visit1_visit2_gene1) <- SDMS_visit1_visit2_gene1[,"Patient"]
SDMS_visit1_visit2_gene1 <- SDMS_visit1_visit2_gene1[,2:25]
SDMS_visit1_visit2_gene2 <- subset(SDMS_visit1_visit2_gene,Visit==2)
rownames(SDMS_visit1_visit2_gene2) <- SDMS_visit1_visit2_gene2[,"Patient"]
SDMS_visit1_visit2_gene2 <- SDMS_visit1_visit2_gene2[,2:25]
SDMS_visit1_visit2_gene <- SDMS_visit1_visit2_gene2 - SDMS_visit1_visit2_gene1
colData_SDMS_unique <- colData_SDMS[,c(5,6,10)]
colData_SDMS_unique <- colData_SDMS_unique[!duplicated(colData_SDMS_unique$Patient), ]
SDMS_visit1_visit2 <- merge(SDMS_visit1_visit2_gene,colData_SDMS_unique,by.x=0,by.y="Patient")
rownames(SDMS_visit1_visit2) <- SDMS_visit1_visit2[,1];SDMS_visit1_visit2 <- SDMS_visit1_visit2[,-1]
write.table(SDMS_visit1_visit2, file="SDMS_visit1_visit2.txt", sep="\t", quote=F, col.names = NA)
write.table(SDMS_visit1_visit2_gene1, file="SDMS_sharedpatients_visit1.txt", sep="\t", quote=F, col.names = NA)
write.table(SDMS_visit1_visit2_gene2, file="SDMS_sharedpatients_visit2.txt", sep="\t", quote=F, col.names = NA)
