setwd("/home/aki/Documents/Rannsóknir/Haaliaeetus/Haliaeetus_ernir/ADMIXTURE/")

cv<-read.csv("logout_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.cv.error",header=F,sep=" " )
colnames(cv)<-c("K","CV")

ggplot(cv)+
  geom_line(aes(K,CV))+
  geom_vline(xintercept = 4)+
  theme_classic()

name <- read.table("logout_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.list")
tbl=read.table("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.6.Q")
row.names(tbl)<-unlist(name)
barplot(t(as.matrix(tbl)), col=rainbow(6),xlab=NA, ylab="Ancestry", border=NA,cex.names = 0.75,las=2)
