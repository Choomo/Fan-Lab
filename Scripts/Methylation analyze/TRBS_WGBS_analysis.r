setwd("D:/DATA/TRBS_WGBS/")
####################################################
###1. Generate input files for DSS 
for (f in 1:length(list.files(pattern="[.]bismark.cov"))){
  cov=read.table(list.files(pattern="[.]bismark.cov")[f])
  cov=cov[,-3]
  cov[,3]=cov[,4]+cov[,5] #total reads
  cov=cov[,-5] ## X:methyl reads
  colnames(cov)=c("chr","pos","N","X")
  filenames=gsub("GSM\\d+[_]WGBS[_]","",list.files(pattern="[.]bismark.cov")[f])
  filenames=gsub("bismark[.]cov","csv",filenames)
  write.csv(cov,file=filenames)
  }
 
 ##2. Call DML & DMR
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DSS")

library(DSS)
require(bsseq)
#path <- file.path("D:/DATA/TRBS_WGBS")
path <- file.path("/mnt/d/DATA/TRBS_WGBS/")
GS.NSC.1 <- read.csv(file.path(path, "GS-NSC1.csv"), header=TRUE)
GS.NSC.2 <- read.csv(file.path(path, "GS-NSC2.csv"), header=TRUE)
#IN.NSC.1 <- read.csv(file.path(path, "IN-NSC1.csv"), header=TRUE)
#IN.NSC.2 <- read.csv(file.path(path, "IN-NSC2.csv"), header=TRUE)
#J1.ESC <- read.csv(file.path(path, "J1-ESC.csv"), header=TRUE)
J1.NSC.1 <- read.csv(file.path(path, "J1-NSC1.csv"), header=TRUE)
J1.NSC.2 <- read.csv(file.path(path, "J1-NSC2.csv"), header=TRUE)
#KO.ESC <- read.csv(file.path(path, "KO-ESC.csv"), header=TRUE)
KO.NSC <- read.csv(file.path(path, "KO-NSC.csv"), header=TRUE)
MK.NSC.1 <- read.csv(file.path(path, "MK-NSC1.csv"), header=TRUE)
MK.NSC.2 <- read.csv(file.path(path, "MK-NSC2.csv"), header=TRUE)
PL.NSC.1 <- read.csv(file.path(path, "PL-NSC1.csv"), header=TRUE)
PL.NSC.2 <- read.csv(file.path(path, "PL-NSC2.csv"), header=TRUE)

GS.NSC.1=GS.NSC.1[,-1]
GS.NSC.2=GS.NSC.2[,-1]
MK.NSC.1 =MK.NSC.1 [,-1]
MK.NSC.2 =MK.NSC.2 [,-1]
J1.NSC.1=J1.NSC.1[,-1]
J1.NSC.2=J1.NSC.2[,-1]
PL.NSC.1=PL.NSC.1[,-1]
PL.NSC.2=PL.NSC.2[,-1]
KO.NSC=KO.NSC[,-1]


BSobj1<- makeBSseqData(list(J1.NSC.1,J1.NSC.2,GS.NSC.1, GS.NSC.2), c("J1.NSC.1","J1.NSC.2","GS.NSC.1","GS.NSC.2"))
BSobj2<- makeBSseqData(list(J1.NSC.1,J1.NSC.2,MK.NSC.1, MK.NSC.2), c("J1.NSC.1","J1.NSC.2","MK.NSC.1","MK.NSC.2"))
BSobj4<- makeBSseqData(list(J1.NSC.1,J1.NSC.2,PL.NSC.1, PL.NSC.2), c("J1.NSC.1","J1.NSC.2","PL.NSC.1","PL.NSC.2"))
BSobj5<- makeBSseqData(list(J1.NSC.1,J1.NSC.2,KO.NSC), c("J1.NSC.1","J1.NSC.2","KO.NSC"))
BSobj3<- makeBSseqData(list(J1.NSC.1,J1.NSC.2,IN.NSC.1, IN.NSC.2), c("J1.NSC.1","J1.NSC.2","IN.NSC.1","IN.NSC.2"))


dmlTest1 <- DMLtest(BSobj1, group1=c("J1.NSC.1", "J1.NSC.2"), group2=c("GS.NSC.1", "GS.NSC.2"))
dmlTest2 <- DMLtest(BSobj2, group1=c("J1.NSC.1", "J1.NSC.2"), group2=c("MK.NSC.1", "MK.NSC.2"))
dmlTest3 <- DMLtest(BSobj3, group1=c("J1.NSC.1", "J1.NSC.2"), group2=c("IN.NSC.1", "IN.NSC.2"))
dmlTest4 <- DMLtest(BSobj4, group1=c("J1.NSC.1", "J1.NSC.2"), group2=c("PL.NSC.1", "PL.NSC.2"))
dmlTest5 <- DMLtest(BSobj5, group1=c("J1.NSC.1", "J1.NSC.2"), group2="KO.NSC")


dmls1 <- callDML(dmlTest1,delta=0.1, p.threshold=0.001)
dmls2 <- callDML(dmlTest2,delta=0.1, p.threshold=0.001)
dmls3 <- callDML(dmlTest3,delta=0.1, p.threshold=0.001)
dmls4 <- callDML  (dmlTest4,delta=0.1, p.threshold=0.001)
dmls5 <- callDML(dmlTest5,delta=0.1, p.threshold=0.001)

dmrs1 <- callDMR(dmlTest1, p.threshold=0.01)
dmrs2 <- callDMR(dmlTest2, p.threshold=0.01)
dmrs3 <- callDMR(dmlTest3, p.threshold=0.01)
dmrs4 <- callDMR(dmlTest4, p.threshold=0.01)
dmrs5 <- callDMR(dmlTest5, p.threshold=0.01)

saveRDS(BSobj1,"")

save.image()

###################
gtf=read.table("mm10.refGene.gene_position.txt",header=T)

dml1=read.csv("DML_GS.csv",header=T,row.names=1)
dml2=read.csv("DML_MK.csv",header=T,row.names=1)
dml3=read.csv("dmls_KO.csv",header=T,row.names=1)
dml4=read.csv("dmls_PL.csv",header=T,row.names=1)

dmr1=read.csv("DMR_GS.csv",header=T,row.names=1)
dmr2=read.csv("DMR_MK.csv",header=T,row.names=1)
dmr3=read.csv("dmrs_KO.csv",header=T)
dmr4=read.csv("dmrs_PL.csv",header=T,row.names=1)


gs1=read.table("GSM3231360_WGBS_GS-NSC1.bismark.cov")
gs1[,3]=paste(gs1[,1],gs1[,2],sep="-")
gs2=read.table("GSM3231361_WGBS_GS-NSC2.bismark.cov")
gs2[,3]=paste(gs2[,1],gs2[,2],sep="-")

wt1=read.table("GSM3231355_WGBS_J1-NSC1.bismark.cov")
wt1[,3]=paste(wt1[,1],wt1[,2],sep="-")
wt2=read.table("GSM3231356_WGBS_J1-NSC2.bismark.cov")
wt2[,3]=paste(wt2[,1],wt2[,2],sep="-")


mk1=read.table("GSM3231362_WGBS_MK-NSC1.bismark.cov")
mk1[,3]=paste(mk1[,1],mk1[,2],sep="-")
MK2=read.table("GSM3231363_WGBS_MK-NSC2.bismark.cov")
MK2[,3]=paste(MK2[,1],MK2[,2],sep="-")
 mk= merge(x=mk1,y=MK2,by="V3")
 mk=mk[,-c(7:8)]
colnames(mk)[c(4,7)]=c("mk1_meth","mk2_meth")

library(dplyr)
gs = merge(x=gs1,y=gs2,by="V3")
gs=gs[,-c(7:8)]
colnames(gs)[c(4,7)]=c("gs1_meth","gs2_meth")
wt = merge(x=wt1,y=wt2,by="V3")
wt=wt[,-c(7:8)]
wt=wt[,-c(2:3)]
colnames(wt)[c(4,7)]=c("wt1_meth","wt2_meth")

all = merge(x=wt,y=mk,by="V3") ##14646847  lines    
apply(all[,c(2,5,10,13)],2,mean)

boxplot(all[,c(2,5,10,13)])

all2=all[apply(all[,c(3,6,11,14)],1,min)>4,] # 1789413  lines with each C site detected reads >=5 

boxplot(all2[,c(2,5,10,13)])
apply(all2[,c(2,5,10,13)],2,mean)

dml1=read.csv("DML_MK.csv",header=T,row.names=1)
dml1[,13]=paste(dml1$chr,dml1$pos,sep="-")
colnames(dml1)[13]="V3"

df = merge(x=all,y=dml1,by="V3",all.x=TRUE) ##left join
saveRDS(df,"WT-MK-ALL methylation and DMLS.rds")
#df2=na.omit(df2)
#apply(df2[,c(2,5,10,13)],2,mean)

dim(df)
df[,28]=apply(df[,c(2,5)],1,mean)
df[,29]=apply(df[,c(10,13)],1,mean)
colnames(df)[28:29]=c("wt_mean","mk_mean")
df2=na.omit(df) ##10692    29


png("WT-MK global methylation comparision.png",width=5,height=5,units="in",res=350)
par(mar = c(5,4,4,5) + .1)
smoothScatter(x=df$wt_mean,y=df$mk_mean, nrpoints = 0,xlab="WT (%methylation)",ylab="M548K (%methylation)")
points(x=df2$wt_mean,y=df2$gs_mean,pch=46,col="black")
dev.off()
library(ggplot2)
library(reshape2)
library(ggpubr)
df3=df2[,c(1,2,5,10,13)]
colnames(df3)[2:3]=c("wt1_meth","wt2_meth")
rownames(df3)=df3[,1]
df3=melt(df3,id="V3")
group=c(rep("WT",21384),rep("M548K",21384))
df3=cbind(group,df3)

compare_means(value~ group, data =df3)

my_comparisons <- list( c("WT", "M548K"))
df3$group=factor(df3$group)

png("WT-MK-difDMLs_boxplot.png",width=6,height=4.5,units="in",res=350)
ggboxplot(df3,x="variable",y="value",fill="group")+
  ylab("DMLs %methylation")+
  xlab("")+
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y =0.02)  +
  theme( strip.text.x = element_text(size=18, face="bold"),
         axis.text.y=element_text(size=14, face="bold"),
         axis.text.x=element_text(size=14, face="bold"))+
  theme_bw()
dev.off()
