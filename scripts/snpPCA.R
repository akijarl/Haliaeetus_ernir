setwd("~/Documents/UI/rannsoknir/Haliaeetus/")
library(gdsfmt)
require(SNPRelate)
require(adegenet)
require(pegas)
require(hierfstat)
require(ggplot2)

vcf.fn<-"all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.vcf"
snpgdsVCF2GDS(vcf.fn, "ernir.gds",method = "copy.num.of.ref")

ernir<-snpgdsOpen("ernir.gds")
#snpgdsClose(ernir)

pcaC <- snpgdsPCA(ernir,autosome.only=F)

pc.percent <- pcaC$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pcaC$sample.id,
                  EV1 = pcaC$eigenvect[,1],
                  EV2 = pcaC$eigenvect[,2],
                  EV3 = pcaC$eigenvect[,3],
                  EV4 = pcaC$eigenvect[,4],
                  #EV5 = pcaC$eigenvect[,5],
                  #EV6 = pcaC$eigenvect[,6],
                  stringsAsFactors = FALSE)
head(tab)

Pop <- c(rep("IS_c",25),rep("NO_c",12),rep("DK_h",5),rep("DK_c",11),rep("EE_c",3),rep("GL_c",12),rep("GL_h",8),rep("IS_h",2),rep("NO_h",13),"TU_h")
Per_exp<-head(round(pc.percent, 2))

tab$pop<-Pop
cl<-c("blue", "cyan","orange", "darkgreen","green","red","magenta", "black","grey","brown")
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E15759","#F28E2B","#4E79A7")
shp<-c(0,1,2,3,4,5,6,7,8,10)
ggplot(tab, aes(EV1,EV2,color=Pop,shape=Pop,label=sample.id)) +
  xlab(paste("PC1 (",Per_exp[1],"%)",sep="")) + 
  ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + 
  geom_point(size=3) + 
  stat_ellipse(level=0.75,size=1)+
  scale_shape_manual(name="Pop", labels=unique(Pop)[order(unique(Pop))], values=shp)+
  scale_color_manual(name="Pop", labels=unique(Pop)[order(unique(Pop))], values=cl)+
  labs(color="") + 
  theme_classic()

#smartPCA
spca_tab<-read.table("../VCF/ernir_LD.eigenvec")
spca_Per_exp<-read.table("../VCF/ernir_LD.eigenval")
spca_tab$V12<-Pop
colnames(spca_tab)<-c("sample.id","EV1","EV2","EV3", "EV4","EV5","EV6","EV7","EV8","EV9","EV10","PopS")
spca_tab$PopS<-factor(spca_tab$PopS,levels = c("TU_h","EE_c","NO_h","DK_h","NO_c", "DK_c","IS_c","IS_h","GL_c","GL_h"))
ggplot(spca_tab, aes(EV1,EV2,color=Pop,shape=Pop,label=sample.id)) +
  xlab(paste("PC1 (",round(spca_Per_exp[1,],2),"%)",sep="")) + 
  ylab(paste("PC2 (",round(spca_Per_exp[2,],2),"%)",sep="")) + 
  geom_point(size=3) + 
  #geom_point(size=3, subset = .(label == 'point'))+
  stat_ellipse(level=0.75,size=1)+
  scale_shape_manual(name="Pop", labels=unique(Pop)[order(unique(Pop))], values=shp)+
  scale_color_manual(name="Pop", labels=unique(Pop)[order(unique(Pop))], values=cl)+
  labs(color="") + 
  theme_classic()

spcaM_tab<-read.table("Mainland_LD.eigenvec")
spcaM_Per_exp<-read.table("Mainland_LD.eigenval")
colnames(spcaM_tab)<-c("sample.id","EV1","EV2","EV3", "EV4","EV5","EV6","EV7","EV8","EV9","EV10","Pop")
#spca_tab$PopS<-factor(spca_tab$PopS,levels = c("TU_h","EE_c","NO_h","DK_h","NO_c", "DK_c","IS_c","IS_h","GL_c","GL_h"))
clM<-c("blue", "cyan","orange", "black","grey","brown")
shpM<-c(0,1,2,7,8,10)
ggplot(spcaM_tab, aes(EV1,EV2,color=Pop,shape=Pop)) +
  xlab(paste("PC1 (",round(spcaM_Per_exp[1,],2),"%)",sep="")) + 
  ylab(paste("PC2 (",round(spcaM_Per_exp[2,],2),"%)",sep="")) + 
  geom_point(size=3) + 
  #geom_point(size=3, subset = .(label == 'point'))+
  stat_ellipse(level=0.75,size=1)+
  scale_shape_manual(name="Pop", labels=unique(spcaM_tab$Pop)[order(unique(spcaM_tab$Pop))], values=shpM)+
  scale_color_manual(name="Pop", labels=unique(spcaM_tab$Pop)[order(unique(spcaM_tab$Pop))], values=clM)+
  labs(color="") + 
  theme_classic()

ibs <- snpgdsIBS(ernir, autosome.only = F, remove.monosnp = F)

loc<-cmdscale(1-ibs$ibs,k=2)
x<-loc[,1];y<-loc[,2]
ls.gdsn(ernir,recursive = T)
Pop<-as.factor(Pop)
Sample<-as.factor(read.gdsn(index.gdsn(ernir, "sample.id")))

summary(Pop)
cl2<-c(rep("blue",11), rep("cyan",5),rep("orange",3), rep("darkgreen",12),rep("green",8),rep("red",25),rep("magenta",2),rep("black",12),rep("grey",13),rep("brown",1))
shp2<-c(rep(0,11), rep(1,5),rep(2,3), rep(3,12),rep(4,8),rep(5,25),rep(6,2),rep(7,12),rep(8,13),rep(10,1))

par(mar=c(5.1,4.1,4.1,8.1),xpd=T)
plot(x,y,col=cl2,pch=shp2,xlab="",cex=1.5,ylab="",main="cmdscale(IBS Distance)")
legend("topright",bty="n",legend=levels(Pop),pch=shp,text.col=cl,inset=c(-0.15,0))

ibs.hc<-snpgdsHCluster(ibs)
set.seed(100)
rv <- snpgdsCutTree(ibs.hc, label.H=TRUE, label.Z=TRUE)
plot(rv$dendrogram, main="Ernir",hang=-1)
table(rv$samp.group)


#Estimating IBD Using Maximum Likelihood Estimation (MLE)
# Estimate IBD coefficients
set.seed(100)
snp.id <- sample(snpset.id, 1500)  # random 1500 SNPs
ibd <- snpgdsIBDMLE(genofile, sample.id=YRI.id, snp.id=snp.id,
                    maf=0.05, missing.rate=0.05, num.thread=2)

# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="YRI samples (MLE)")
lines(c(0,1), c(1,0), col="red", lty=2)