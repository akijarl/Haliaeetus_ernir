setwd("/home/aki/Documents/Rannsóknir/Haaliaeetus/Haliaeetus_ernir/ADMIXTURE/")
require(ggplot2)

cv<-read.csv("logout_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.cv.error",header=F,sep=" " )
colnames(cv)<-c("K","CV")

ggplot(cv)+
  geom_line(aes(K,CV))+
  geom_vline(xintercept = 4,linetype="dashed")+
  scale_x_continuous(breaks = c(1:15))+
  theme_classic()

name <- read.table("logout_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.list")
tbl2=read.table("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.2.Q")
row.names(tbl2)<-unlist(name)
tbl2$pop <- c(rep("IS_C",25),rep("NO_C",12),rep("DK_H",5),rep("DK_C",11),rep("EE_C",3),rep("GL_C",12),rep("GL_H",8),rep("IS_H",2),rep("NO_H",13),"TU_H")
tbl2<-tbl2[c(which(tbl2$pop=="IS_C"),which(tbl2$pop=="IS_H"),which(tbl2$pop=="GL_C"),which(tbl2$pop=="GL_H"),which(tbl2$pop=="NO_C"),which(tbl2$pop=="NO_H"),which(tbl2$pop=="DK_C"),which(tbl2$pop=="DK_H"),which(tbl2$pop=="EE_C"),which(tbl2$pop=="TU_H")),]

barplot(t(as.matrix(tbl2)), col=c("blue","red"),xlab=NA, ylab="Ancestry", border=NA,cex.names = 0.75,las=2)

tbl3=read.table("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.3.Q")
row.names(tbl3)<-unlist(name)
tbl3$pop <- c(rep("IS_C",25),rep("NO_C",12),rep("DK_H",5),rep("DK_C",11),rep("EE_C",3),rep("GL_C",12),rep("GL_H",8),rep("IS_H",2),rep("NO_H",13),"TU_H")
tbl3<-tbl3[c(which(tbl3$pop=="IS_C"),which(tbl3$pop=="IS_H"),which(tbl3$pop=="GL_C"),which(tbl3$pop=="GL_H"),which(tbl3$pop=="NO_C"),which(tbl3$pop=="NO_H"),which(tbl3$pop=="DK_C"),which(tbl3$pop=="DK_H"),which(tbl3$pop=="EE_C"),which(tbl3$pop=="TU_H")),]
barplot(t(as.matrix(tbl3)), col=c("red","blue","orange"),xlab=NA, yaxt="n", border=NA,cex.names = 0.75,las=2)

tbl4=read.table("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.4.Q")
row.names(tbl4)<-unlist(name)
tbl4$pop <- c(rep("IS_C",25),rep("NO_C",12),rep("DK_H",5),rep("DK_C",11),rep("EE_C",3),rep("GL_C",12),rep("GL_H",8),rep("IS_H",2),rep("NO_H",13),"TU_H")
tbl4<-tbl4[c(which(tbl4$pop=="IS_C"),which(tbl4$pop=="IS_H"),which(tbl4$pop=="GL_C"),which(tbl4$pop=="GL_H"),which(tbl4$pop=="NO_C"),which(tbl4$pop=="NO_H"),which(tbl4$pop=="DK_C"),which(tbl4$pop=="DK_H"),which(tbl4$pop=="EE_C"),which(tbl4$pop=="TU_H")),]
barplot(t(as.matrix(tbl4)), col=c("green","orange","red","blue"),xlab=NA,ylab="Ancestry", border=NA,cex.names = 0.75,las=2,cex.axis = 1.5,cex.lab=1.5)

tbl5=read.table("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.5.Q")
row.names(tbl5)<-unlist(name)
tbl5$pop <- c(rep("IS_C",25),rep("NO_C",12),rep("DK_H",5),rep("DK_C",11),rep("EE_C",3),rep("GL_C",12),rep("GL_H",8),rep("IS_H",2),rep("NO_H",13),"TU_H")
tbl5<-tbl5[c(which(tbl5$pop=="IS_C"),which(tbl5$pop=="IS_H"),which(tbl5$pop=="GL_C"),which(tbl5$pop=="GL_H"),which(tbl5$pop=="NO_C"),which(tbl5$pop=="NO_H"),which(tbl5$pop=="DK_C"),which(tbl5$pop=="DK_H"),which(tbl5$pop=="EE_C"),which(tbl5$pop=="TU_H")),]
barplot(t(as.matrix(tbl5)), col=c("blue","red","purple","green","orange"),xlab=NA, ylab="Ancestry", border=NA,cex.names = 0.75,las=2,cex.axis = 1.5,cex.lab=1.5)

tbl6=read.table("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.6.Q")
row.names(tbl6)<-unlist(name)
tbl6$pop <- c(rep("IS_C",25),rep("NO_C",12),rep("DK_H",5),rep("DK_C",11),rep("EE_C",3),rep("GL_C",12),rep("GL_H",8),rep("IS_H",2),rep("NO_H",13),"TU_H")
tbl6<-tbl6[c(which(tbl6$pop=="IS_C"),which(tbl6$pop=="IS_H"),which(tbl6$pop=="GL_C"),which(tbl6$pop=="GL_H"),which(tbl6$pop=="NO_C"),which(tbl6$pop=="NO_H"),which(tbl6$pop=="DK_C"),which(tbl6$pop=="DK_H"),which(tbl6$pop=="EE_C"),which(tbl6$pop=="TU_H")),]
barplot(t(as.matrix(tbl6)), col=c("blue","purple","red","pink","orange","green"),xlab=NA, yaxt="n", border=NA,cex.names = 0.75,las=2)


K6_IS<-tbl6[tbl6$pop=="IS_C"|tbl6$pop=="IS_H",]
colnames(K6_IS)<-c("Denmark","Estonia","Iceland_1","Iceland_2","Norway","Greenland","Pop")
str(K6_IS)
summary(K6_IS)

write.table(K6_IS,"K6_IS.csv",sep="\t")
