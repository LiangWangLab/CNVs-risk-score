library(devtools)
library(ComplexHeatmap)
library(dplyr)
x=1
y=1
w=1
h=1
####Oncoprint - visit 1####
col = c(Amp="indianred1", Del="deepskyblue1")
alter_fun = list(
  Amp = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                       gp = gpar(fill = col["Amp"], col = "black")),
  Del = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                       gp = gpar(fill = col["Del"], col = "black"))
)
####promote####
pv1<-read.table("pv1_indel.txt",header = T)
pv1<-pv1[,1:24]
pv1[]<-sapply(pv1, function(x) gsub("NC","",x))
pv1 <-t(pv1)
sample_order <- names(sort(apply(pv1, 2, function(x) sum(x==""))))
oncoPrint(pv1, alter_fun = alter_fun, col = col, top_annotation =NULL, right_annotation =NULL,
          column_title = "PROMOTE visit 1 CNVs - Genes of interest", column_order = sample_order)

####SDMS####
sv1<-read.table("sv1_indel.txt",header = T)
sv1<-sv1[,1:24]
sv1[]<-sapply(sv1, function(x) gsub("NC","",x))
sv1 <-t(sv1)
sample_order <- names(sort(apply(sv1, 2, function(x) sum(x==""))))
oncoPrint(sv1, alter_fun = alter_fun, col = col, top_annotation =NULL, right_annotation =NULL,
          column_title = "SDMS visit 1 CNVs - Genes of interest", column_order = sample_order)

####Oncoprint - visit 2####
col = c(Amp="indianred1", Del="deepskyblue1")
alter_fun = list(
  Amp = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                       gp = gpar(fill = col["Amp"], col = "black")),
  Del = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                       gp = gpar(fill = col["Del"], col = "black"))
)
####promote####
pv2<-read.table("pv2_indel.txt",header = T)
pv2<-pv2[,1:24]
pv2[]<-sapply(pv2, function(x) gsub("NC","",x))
pv2 <-t(pv2)
sample_order <- names(sort(apply(pv2, 2, function(x) sum(x==""))))
oncoPrint(pv2, alter_fun = alter_fun, col = col, top_annotation =NULL, right_annotation =NULL,
          column_title = "PROMOTE visit 2 CNVs - Genes of interest", column_order = sample_order)

####SDMS####
sv2<-read.table("sv2_indel.txt",header = T)
sv2<-sv2[,1:24]
sv2[]<-sapply(sv2, function(x) gsub("NC","",x))
sv2 <-t(sv2)
sample_order <- names(sort(apply(sv2, 2, function(x) sum(x==""))))
oncoPrint(sv2, alter_fun = alter_fun, col = col, top_annotation =NULL, right_annotation =NULL,
          column_title = "SDMS visit 2 CNVs - Genes of interest", column_order = sample_order)

####Oncoprint - two visits patients####
col = c(Amp1="indianred1", Amp2="indianred1", Del1="deepskyblue1", Del2="deepskyblue1")
alter_fun = list(
  Amp1 = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = col["Amp1"], col = "black"))},
  Amp2 = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col["Amp2"], col = "black"))},
  Del1 = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = col["Del1"], col = "black"))},
  Del2 = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col["Del2"], col = "black"))}
)
####promote####
promote <- read.table("promote_final.txt", header = T, sep = "\t")
promote <-promote[,c(1:3)]
pv1<-read.table("pv1_indel.txt",header = T)
pv1<-merge(pv1,promote,by.x = 0,by.y="v1ID")
row.names(pv1)=pv1$Patient
pv1<-pv1[,c(2:25)]
pv2<-read.table("pv2_indel.txt",header = T)
pv2<-merge(pv2,promote,by.x = 0,by.y="v2ID")
row.names(pv2)=pv2$Patient
pv2<-pv2[,c(2:25)]

#oncoprint matrix
pv1[]<-sapply(pv1, function(x) gsub("Amp","Amp1",x))
pv1[]<-sapply(pv1, function(x) gsub("Del","Del1",x))
pv2[]<-sapply(pv2, function(x) gsub("Amp","Amp2",x))
pv2[]<-sapply(pv2, function(x) gsub("Del","Del2",x))

promote_indel <- matrix("",82,24)
colnames(promote_indel)<-gsub("v1.","",colnames(pv1))
row.names(promote_indel)<-rownames(pv1)
for (i in 1:24) {
  promote_indel[,i]<-paste0(pv1[,i],";", pv2[,i],sep="")
}
promote_indel[]<-sapply(promote_indel, function(x) gsub("NC;NC","",x))
promote_indel[]<-sapply(promote_indel, function(x) gsub("NC;","",x))
promote_indel[]<-sapply(promote_indel, function(x) gsub(";NC","",x))

promote_indel <-t(promote_indel)

sample_order <- names(sort(apply(promote_indel, 2, function(x) sum(x==""))))
oncoPrint(promote_indel, alter_fun = alter_fun, col = col, top_annotation =NULL, right_annotation =NULL,
          column_title = "PROMOTE two visits CNVs - Genes of interest", column_order = sample_order)

####SDMS####
sv1<-read.table("sv1_indel.txt",header = T)
sv2<-read.table("sv2_indel.txt",header = T)
coldata<-read.table("colData_two_visit.txt", header = T, sep = "\t")
coldata<-subset(coldata,Group=="SDMS")
sv1<-sv1[rownames(sv1) %in% coldata$PatientID,]
sv2<-sv2[rownames(sv2) %in% coldata$PatientID,]
#oncoprint matrix
sv1[]<-sapply(sv1, function(x) gsub("Amp","Amp1",x))
sv1[]<-sapply(sv1, function(x) gsub("Del","Del1",x))
sv2[]<-sapply(sv2, function(x) gsub("Amp","Amp2",x))
sv2[]<-sapply(sv2, function(x) gsub("Del","Del2",x))
SDMS_indel <- matrix("",33,24)
colnames(SDMS_indel)<-colnames(sv1)[1:24]
row.names(SDMS_indel)<-rownames(sv1)
for (i in 1:24) {
  SDMS_indel[,i]<-paste0(sv1[,i],";", sv2[,i],sep="")
}
SDMS_indel[]<-sapply(SDMS_indel, function(x) gsub("NC;NC","",x))
SDMS_indel[]<-sapply(SDMS_indel, function(x) gsub("NC;","",x))
SDMS_indel[]<-sapply(SDMS_indel, function(x) gsub(";NC","",x))

SDMS_indel <-t(SDMS_indel)


sample_order <- names(sort(apply(SDMS_indel, 2, function(x) sum(x==""))))
oncoPrint(SDMS_indel, alter_fun = alter_fun, col = col, top_annotation =NULL, right_annotation =NULL,
          column_title = "SDMS two visits CNVs - Genes of interest", column_order = sample_order)
