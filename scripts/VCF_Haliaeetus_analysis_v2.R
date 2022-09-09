#####################################################
#Code for Haliaeetus LD fitlered data
#####################################################

#setwd("/home/aki/Documents/Rannsóknir/Haaliaeetus/VCF")
setwd("/Users/akijarl/Documents/UI/rannsoknir/Haliaeetus/")

require(cowplot)
require(vcfR)
require(ggplot2)
require(ggrepel)
require(seqinr)

#require(BiocManager)
#BiocManager::install("qvalue")
#require(devtools)
#devtools::install_github("whitlock/OutFLANK")
require(OutFLANK)

#BiocManager::install("SNPRelate")
require(SNPRelate)

#remotes::install_github("privefl/bigsnpr")
require(bigsnpr)

#install.packages("robust")
require(robust)

#install.packages("bigstatsr")
require(bigstatsr)

X<-read.vcfR("Haliaeetus_LDfilt_snpsOnly.vcf")
#Y<-maf(X, element=1)
#y<-maf(X, element=2)

#y_o<-y[order(y[,"Frequency"],decreasing = T),]
#barplot(y_o[,"Frequency"], names.arg="", main="Minor allele frequency")

#Y_o<-Y[order(Y[,"Frequency"],decreasing = T),]
#barplot(Y_o[,"Frequency"], names.arg="", main="Major allele frequency")

queryMETA(X)
queryMETA(X, element = 'DP')

dp <- extract.gt(X, element = "DP", as.numeric=TRUE)
#barplot(dp, ylab = "Read depth (DP)")
#abline(h=min(dp),col="red")
#abline(h=summary(dp)[[2]],col="orange")
#abline(h=mean(dp),col="orange")
#abline(h=median(dp),col="cyan")
#abline(h=summary(dp)[[5]],col="orange")
#abline(h=max(dp),col="red")

colMeans(dp)
summary(colMeans(dp))

#contemporary depth
summary(colMeans(dp[,c(1:25,26:37,43:53,54:56,57:68)]))

#historical depth
summary(colMeans(dp[,c(-1:-25,-26:-37,-43:-53,-54:-56,-57:-68)]))
#gt <- extract.gt(X)
#hets <- is_het(gt)

#Annotated contig names in TS_minE_best
#Y<-X[which(X@fix[,1]%in%Annot_minE_best$TR),]
#Y<-X[which(X@fix[,1]%in%Annot_minE_best_TSL$TR),]
#Y<-X[which(X@fix[,1]%in%Annot_minE_best_NRecip_TSL$TR),]

#dp <- extract.info(Y, element = "DP", as.numeric=TRUE)
#barplot(dp, ylab = "Read depth (DP)")
#summary(dp)

geno <- extract.gt(X) # Character matrix containing the genotypes
position <- getPOS(X) # Positions in bp
chromosome <- getCHROM(X) # Chromosome information
chromosome_n<-as.numeric(factor(chromosome))
pos_loc<-paste(chromosome_n,position,sep="_")

tab<-read.table("ernir_LD.eigenvec")
pop<-tab$V12

#test <- geno[!geno %in% c("0/0","0/1", "1/0","1/1", "0/2", "2/0","2/2","1/2", "2/1","2/3", "0/3", "1/3")]
#test[!is.na(test)]

mult<-NULL
for(i in 1:nrow(geno)){
  if(any(geno[i,]%in%c("0/2","2/2","1/2","0/3","1/3","2/3"))){
    mult<-c(mult,i)
  }
}

geno_mul<-geno[mult,]

colnames(geno)
table(as.vector(geno[,1:25])) #Iceland
table(as.vector(geno[,26:37])) #Norway
table(as.vector(geno[,38:42])) #Denmark_H
table(as.vector(geno[,43:53])) #Denmark
table(as.vector(geno[,54:56])) #Estonia
table(as.vector(geno[,57:68])) #Greenland
table(as.vector(geno[,69:76])) #Greenland_H
table(as.vector(geno[,77:78])) #Iceland_H
table(as.vector(geno[,79:91])) #Norway_H
table(as.vector(geno[,92])) #Turkey_H

table(as.vector(geno))
sum(is.na(geno))/(nrow(geno)*ncol(geno))

geno_c<-geno[,c(1:25,26:37,43:53,54:56,57:68)]
colnames(geno_c)

sum(is.na(geno_c))/(nrow(geno_c)*ncol(geno_c))

nll_c<-NULL
for(i in 1:nrow(geno_c)){
  if(length(unique(geno_c[i,]))==1){
    nll_c<-c(nll_c,i)
  }
}

table(as.vector(geno_c[-nll_c,]))

geno_h<-geno[,c(-1:-25,-26:-37,-43:-53,-54:-56,-57:-68)]
colnames(geno_h)

sum(is.na(geno_h))/(nrow(geno_h)*ncol(geno_h))

nll_h<-NULL
for(i in 1:nrow(geno_h)){
  if(length(unique(geno_h[i,]))==1){
    nll_h<-c(nll_h,i)
  }
}

table(as.vector(geno_h[-nll_h,]))


geno_sing<-geno[-mult,]
table(as.vector(geno_sing))

nll<-NULL
for(i in 1:nrow(geno_sing)){
  if(length(unique(geno_sing[i,]))==1){
    nll<-c(nll,i)
  }
}

table(as.vector(geno[-nll]))

geno_sing<-geno_sing[-nll,]

sum(is.na(geno_sing))/(nrow(geno_sing)*ncol(geno_sing))

geno_sing_c<-geno_sing[,c(1:25,26:37,43:53,54:56,57:68)]
colnames(geno_sing_c)
length(row.names(geno_sing_c))
length(pos_loc)
pos_loc<-pos_loc[-mult]
pos_loc<-pos_loc[-nll]

G <- matrix(NA, nrow = nrow(geno_sing_c), ncol = ncol(geno_sing_c),dimnames = list(pos_loc,colnames(geno_sing_c)) )

G[geno_sing_c %in% c("0/0")] <- 0
G[geno_sing_c %in% c("0/1", "1/0")] <- 1
G[geno_sing_c %in% c("1/1")] <- 2
# G[geno %in% c("0/2", "2/0", "2|0", "0|2")] <- 3
# G[geno %in% c("1/2", "2/1", "2|1", "1|2")] <- 4
# G[geno %in% c("2/2", "2|2")] <- 5
# G[geno %in% c("0/3", "3/0", "3|0", "0|3")] <- 6
# G[geno %in% c("1/3", "3/1", "3|1", "1|3")] <- 7
#G[geno %in% c("2/3", "3/2", "3|2", "2|3")] <- 8
#G[geno %in% c("3/3", "3|3")] <- 9

table(as.vector(G))

#test2 <- G[!G %in% c(0:9)]
#test2[!is.na(test2)]

#(sum(G %in% c(0:9))+sum(is.na(G[!G %in% c(0:9)]))) / (nrow(G)*ncol(G))

#G[is.na(G)]<-9



#table(as.vector(geno[nll,]))


#dat<-list(chromosome,position,G,pop,ind)

SNPmat<-(t(G))

#colnames(SNPmat)<-pos_loc
G_nM<-SNPmat[ , apply(SNPmat, 2, function(x) !any(is.na(x)))]
chrom_fil<-colnames(G_nM)[colnames(G_nM)%in%pos_loc]

chrom_filt<-as.numeric(sub("_.*", "", chrom_fil))
pos_filt<-as.numeric(sub(".*_", "", chrom_fil))

#modified function MakeDiploidFSTMat
MakeDiploidFSTMat<-function(SNPmat,locusNames,popNames){
  locusname <- unlist(locusNames)
  popname <- unlist(popNames)
  snplevs <- levels(as.factor(unlist(SNPmat)))
  if(any(!(snplevs%in%c(0,1,2,9)))==TRUE) {
    print("Error: Your snp matrix has a character other than 0,1,2 or 9")
    break
  }
  if (dim(SNPmat)[1] != length(popname)) {
    print("Error: your population names do not match your SNP matrix")
    break
  }
  if (dim(SNPmat)[2] != length(locusname)) {
    print("Error:  your locus names do not match your SNP matrix")
    break
  }
  writeLines("Calculating FSTs, may take a few minutes...")
  nloci <- length(locusname)
  FSTmat <- matrix(NA, nrow = nloci, ncol = 8)
  for (i in 1:nloci) {
    FSTmat[i, ] = unlist(getFSTs_diploids(popname, SNPmat[,i]))
    if (i%%10000 == 0) {
      print(paste(i, "done of", nloci))
    }
  }
  outTemp = as.data.frame(FSTmat)
  outTemp = cbind(locusname, outTemp)
  colnames(outTemp) = c("LocusName", "He", "FST", "T1", "T2", 
                        "FSTNoCorr", "T1NoCorr", "T2NoCorr", "meanAlleleFreq")
  return(outTemp)
}

getFSTs_diploids = function(popNameList, SNPDataColumn){  
  #eliminating the missing data for this locus
  popnames=unlist(as.character(popNameList))
  popNameTemp=popnames[which(SNPDataColumn!=9)]
  snpDataTemp=SNPDataColumn[SNPDataColumn!=9]
  
  HetCounts <- tapply(snpDataTemp, list(popNameTemp,snpDataTemp), length)
  HetCounts[is.na(HetCounts)] = 0
  
  #Case: all individuals are genetically identical at this locus
  if(dim(HetCounts)[2]==1){
    return (list(He=NA,FST=NA, T1=NA, T2=NA,FSTNoCorr=NA, T1NoCorr=NA, T2NoCorr=NA,meanAlleleFreq = NA))
  }
  
  if(dim(HetCounts)[2]==2){
    if(paste(colnames(HetCounts),collapse="")=="01"){HetCounts=cbind(HetCounts,"2"=0)}
    if(paste(colnames(HetCounts),collapse="")=="12"){HetCounts=cbind("0"=0,HetCounts)} 
    if(paste(colnames(HetCounts),collapse="")=="02"){HetCounts=cbind(HetCounts[,1],"1"=0, HetCounts[,2])}
  }
  
  out = WC_FST_Diploids_2Alleles(HetCounts)	
  return(out)
}

pop_c <- pop[c(1:25,26:37,43:53,54:56,57:68)]
pop_h <- pop[c(38:42,69:76)]
#my_fst <- MakeDiploidFSTMat(t(G), locusNames = position, popNames = pop)
my_fst <- MakeDiploidFSTMat(G_nM, locusNames = pos_filt, popNames = pop_c)

plot(my_fst$He, my_fst$FST, xlab="Heterozygosity", ylab=expression(paste("F"[ST], " across all populations", sep="")), main=expression(paste("Per locus F"[ST], " vs Heterozygosity", sep="")))

# If a locus has a much lower sample size compared to the rest, it could have a broader error distribution (and therefore incorrectly inferred as an outlier).
plot(my_fst$FST, my_fst$FSTNoCorr)
abline(0,1)

#### SNP trimming from OutFLANK vignette ####
# identify a quasi-independent set of SNPs to calculate FSTbar and df. A common way of obtaining these SNPs is to thin for linkage disequilibrium (SNP thinning), which typically moves along a genome in a sliding window and thins SNPs based on linkage disequilibrium with each other. 
#This may be based on a combination of 
#(i) “pruning,” which sequentially scans the genome and performs pairwise thinning based on a given threshold of correlation, 
#(ii) “clumping,” which may incorporate some information about the importance of SNPs based on summary statistics, and 
#(iii) removing SNPs in long-range LD regions (Prive et al. 2017).

Ga<-add_code256(big_copy(G_nM,type="raw"),code=bigsnpr:::CODE_012)
ind.keep<-snp_clumping(G=Ga,infos.chr=chrom_filt ,infos.pos=pos_filt)
m <- ncol(Ga)
length(ind.keep) / m

out_trim <- OutFLANK(my_fst[ind.keep,], NumberOfSamples=length(unique(pop)), qthreshold = 0.001, Hmin = 0.1)
str(out_trim)
head(out_trim$results)
summary(out_trim$results$OutlierFlag)
summary(out_trim$results$pvaluesRightTail)

OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.01, binwidth = 0.01, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
## Zoom in on right tail
OutFLANKResultsPlotter(out_trim , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.01, binwidth = 0.01, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)

hist(out_trim$results$pvaluesRightTail)

P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)

(my_out <- P1$OutlierFlag==TRUE)
plot(P1$He, P1$FST, pch=19, col=rgb(0,0,0,0.1),xlab="Heterozygosity", ylab=expression(paste("F"[ST], " across all populations", sep="")), main=expression(paste("Per locus F"[ST], " vs Heterozygosity", sep="")))
points(P1$He[my_out], P1$FST[my_out], col="blue")

hist(P1$pvaluesRightTail)

OutLoc<-P1[P1$OutlierFlag==TRUE,]
OutLoc<-OutLoc[!is.na(OutLoc$OutlierFlag),]

chrom_out<-as.numeric(sub("_.*", "", OutLoc$LocusName))
pos_out<-as.numeric(sub(".*_", "", OutLoc$LocusName))
Annot_out<-NULL
for(i in chrom_out){
  Annot_out<-c(Annot_out,levels(factor(chromosome))[i])
}

# OutLoc_annot<-data.frame(Annot_out,chrom_out,pos_out,OutLoc[,-1])
# OutAnnot<-Annot_minE_best_NRecip_TSL[Annot_minE_best_NRecip_TSL$TR%in%OutLoc_annot$Annot_out]


# setwd("/home/aki/Documents/Rannsóknir/Haaliaeetus/Haliaeetus_ernir/")
# 
# het_pops<-read.table("all_IGND_noSMALL")

#setwd("/home/aki/Documents/Rannsóknir/Haaliaeetus/VCF/")
setwd("/home/aki/Documents/Rannsóknir/Haaliaeetus/Haliaeetus_ernir/")

require(tidyverse)
require(vcfR)


hwe <- fread("LDfilt_hardy.hwe")

het <- fread("LDfilt_HET_all.het")

X<-read.vcfR("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.vcf")
all_pops<-colnames(X@gt)[-1]

all_pops[!all_pops%in%het$INDV]

Pop <- c(rep("IS_C",25),rep("NO_C",12),rep("DK_H",5),rep("DK_C",11),rep("EE_C",3),rep("GL_C",12),rep("GL_H",8),rep("IS_H",2),rep("NO_H",13),"TU_H")

het$POP<-Pop

head(het)
het$O.HET<-het$N_SITES-het$`O(HOM)`
het$E.HET<-het$N_SITES-het$`E(HOM)`
het$hetO_frac<-het$O.HET/het$N_SITES
het$POP<-factor(het$POP, levels = c("GL_C","GL_H","IS_C","IS_H","NO_C","NO_H","DK_C","DK_H","EE_C","TU_H"))
het[order(het$O.HET),]
het[order(het$het_frac),]

par(mar=c(5.2,5.2,0.2,0.2))
boxplot(het$O.HET/het$N_SITES~het$POP, ylab = "Observed heterozygosity", xlab = "", cex.lab=3, cex.axis=2)
mtext("Population_time", side = 1, line = 3.5, cex = 3 )

par(mar=c(5.2,5.2,0.2,0.2))
boxplot(het$E.HET/het$N_SITES~het$POP, ylab = "Expected heterozygosity", xlab = "", cex.lab=3, cex.axis=2)
mtext("Population_time", side = 1, line = 3.5, cex = 3 )

par(mar=c(5.2,5.2,0.2,0.2))
boxplot(het$F~het$POP, ylab = "Inbreeding coefficient", xlab = "", cex.lab=3, cex.axis=2)
mtext("Population_time", side = 1, line = 3.5, cex = 3 )

#het_obs_tmp<-het[,c(1,2,3,4,6,7)]
#het_data_obs_exp_for_ggplot<-het_obs_tmp[order(het_obs_tmp$ORDER),]
het_data_obs_exp_for_ggplot<-het[,c(1,2,3,4,6,7)]
het_data_obs_exp_for_ggplot$obs_or_exp<-rep("obs", 92)
colnames(het_data_obs_exp_for_ggplot)[6]<-c("HET")

#het_exp_tmp<-het_data[,c(1,2,3,4,6,8)]
het_exp_tmp2<-het[,c(1,2,3,4,6,8)]
het_exp_tmp2$obs_or_exp<-rep("exp", 92)
colnames(het_exp_tmp2)[6]<-c("HET")

het_inb_tmp2<-het[,c(1,2,3,4,6,5)]
het_inb_tmp2$obs_or_exp<-rep("inb", 92)
colnames(het_inb_tmp2)[6]<-c("HET")

het_data_obs_exp_for_ggplot_true<-rbind(het_data_obs_exp_for_ggplot, het_exp_tmp2)

ggplot_het_bothexpandobsandinb_090421<-
  ggplot(het_data_obs_exp_for_ggplot_true, aes(x=POP, y=HET/N_SITES, color = obs_or_exp)) +
  geom_boxplot() +
  geom_boxplot(data = het_inb_tmp2, aes(x=POP, y=HET, color = obs_or_exp), fill="#636363") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        title = element_text(size = 16), legend.text = element_text(size = 14)) +
  #annotate(geom= "text", x=seq_len(unique(het_data_obs_exp_for_ggplot_true$CO_TI)), y=10, label=het_data_obs_exp_for_ggplot_true$obs_or_exp) +
  scale_color_manual(values = c( "#aaaaaa", "#000000", "#000000")) + 
  labs(x = "Population_time", 
       y = "Heterozygosity & Inbreeding coefficient",
       color="") + 
  theme(legend.position = "none")

ggsave(plot = ggplot_het_bothexpandobsandinb_090421, filename = "het_plot_bothhetandobsandinb_130421.png", width = 30, height = 20, units = "cm", device = "png")

het_gg_observed_090421<-
  ggplot(het, aes(x=POP, y=O.HET/N_SITES)) + 
  geom_boxplot() + 
  labs(title = "", 
       x = "Population_time", 
       y = "Observed heterozygosity") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        title = element_text(size = 16), legend.text = element_text(size = 14)) +
  scale_y_continuous(position = "left", limits = c(0.125, 0.3)) 

het_gg_expected_090421<-
  ggplot(het, aes(x=POP, y=E.HET/N_SITES), ) + 
  geom_boxplot() + 
  labs(title = "", 
       x = "Population_time", 
       y = "Expected heterozygosity") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        title = element_text(size = 16), legend.text = element_text(size = 14)) + 
  scale_y_continuous(position = "right", limits = c(0.125, 0.3))

het_gg_2figures_090421<-
  ggarrange(het_gg_observed_090421, het_gg_expected_090421,
            labels = c("A", "B"),
            ncol = 2, nrow = 1, common.legend = TRUE)

het[, .(meanFIS=mean(F)), by=POP]
het[, .(meanFIS=var(F)), by=POP]

het[, .(meanE.HET=mean(E.HET)), by=POP]
het[, .(meanE.HET=var(E.HET)), by=POP]

het[, .(meanO.HET=mean(O.HET)), by=POP]
het[, .(meanO.HET=var(O.HET)), by=POP]

het[, .(meanHETfrac=mean(het_frac)), by=POP]

het_data_obs_exp_for_ggplot_true
het[, .(meanE.HET=mean(E.HET)), by=POP]


ernir<-snpgdsOpen("ernir.gds")
#snpgdsClose(ernir)
ibs <- snpgdsIBS(ernir, autosome.only = F, remove.monosnp = F)

IBS<-data.frame(ibs$ibs)
colnames(IBS)<-all_pops
rownames(IBS)<-all_pops

for(i in 1:ncol(IBS)){
    IBS[which(IBS[,i]==1.0000000),i]<-NA
}

POP<-as.factor(Pop)
Pop_col<-c(rep("red",25),rep("orange",12),rep("cyan",5),rep("blue",11),rep("black",3),rep("green",12),rep("lightgreen",8),rep("pink",2),rep("yellow",13),"grey")
loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
race <- as.factor(read.gdsn(index.gdsn(genofile, "sample.annot/pop.group")))
plot(x, y, col=Pop_col, xlab = "", ylab = "", main = "cmdscale(IBS Distance)",pch=16,cex=2)
legend("topleft", legend=unique(Pop), text.col=c("red", "orange", "cyan","blue","black","green","lightgreen","pink","yellow","grey"))


summary(as.factor(Pop))
mean(colMeans(IBS[Pop=="GL_C",Pop=="GL_C"],na.rm = T))
mean(colMeans(IBS[Pop=="GL_H",Pop=="GL_H"],na.rm = T))
mean(colMeans(IBS[Pop=="IS_C",Pop=="IS_C"],na.rm = T))
mean(colMeans(IBS[Pop=="IS_H",Pop=="IS_H"],na.rm = T))
mean(colMeans(IBS[Pop=="NO_C",Pop=="NO_C"],na.rm = T))
mean(colMeans(IBS[Pop=="NO_H",Pop=="NO_H"],na.rm = T))
mean(colMeans(IBS[Pop=="DK_C",Pop=="DK_C"],na.rm = T))
mean(colMeans(IBS[Pop=="DK_H",Pop=="DK_H"],na.rm = T))
mean(colMeans(IBS[Pop=="EE_C",Pop=="EE_C"],na.rm = T))
mean(colMeans(IBS[Pop=="TU_H",Pop=="TU_H"],na.rm = T))


hweGC <- read.table("LDfilt_hardy_GC.hwe",header = T)

hardy_GL_C<-read.table(text = gsub("/", "\t", readLines("LDfilt_hardy_GC.hwe")), header = T)
hardy_GL_C$CHR <- as.numeric(revalue(hardy_GL_C$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_GL_C$bin <- hardy_GL_C$POS%/%100000
hardy_GL_C$hetero<-hardy_GL_C$HET/(hardy_GL_C$OBS.HOM1+hardy_GL_C$HET+hardy_GL_C$HOM2.)
hardy_GL_C_tapp_min<-tapply(hardy_GL_C$hetero, list(hardy_GL_C$CHR, hardy_GL_C$bin), quantile, p=0.05)
hardy_GL_C_tapp_q1<-tapply(hardy_GL_C$hetero, list(hardy_GL_C$CHR, hardy_GL_C$bin), quantile, p=0.25)
hardy_GL_C_tapp_q2<-tapply(hardy_GL_C$hetero, list(hardy_GL_C$CHR, hardy_GL_C$bin), quantile, p=0.5)
hardy_GL_C_tapp_q3<-tapply(hardy_GL_C$hetero, list(hardy_GL_C$CHR, hardy_GL_C$bin), quantile, p=0.75)
hardy_GL_C_tapp_max<-tapply(hardy_GL_C$hetero, list(hardy_GL_C$CHR, hardy_GL_C$bin), quantile, p=0.95)
new.hardy_GL_C_tapp_min=as.data.frame(cbind(rep(row.names(hardy_GL_C_tapp_min),each=855),rep(0:854,each=nrow(hardy_GL_C_tapp_min)),as.numeric(t(hardy_GL_C_tapp_min))))
new.hardy_GL_C_tapp_q1=as.data.frame(cbind(rep(row.names(hardy_GL_C_tapp_q1),each=855),rep(0:854,each=nrow(hardy_GL_C_tapp_q1)),as.numeric(t(hardy_GL_C_tapp_q1))))
new.hardy_GL_C_tapp_q2=as.data.frame(cbind(rep(row.names(hardy_GL_C_tapp_q2),each=855),rep(0:854,each=nrow(hardy_GL_C_tapp_q2)),as.numeric(t(hardy_GL_C_tapp_q2))))
new.hardy_GL_C_tapp_q3=as.data.frame(cbind(rep(row.names(hardy_GL_C_tapp_q3),each=855),rep(0:854,each=nrow(hardy_GL_C_tapp_q3)),as.numeric(t(hardy_GL_C_tapp_q3))))
new.hardy_GL_C_tapp_max=as.data.frame(cbind(rep(row.names(hardy_GL_C_tapp_max),each=855),rep(0:854,each=nrow(hardy_GL_C_tapp_max)),as.numeric(t(hardy_GL_C_tapp_max))))
bins<-rep(0:854,length.out=22230)
table_het_per_scaf_GL_C<-new.hardy_GL_C_tapp_min
table_het_per_scaf_GL_C$V2<-bins
table_het_per_scaf_GL_C$V4<-new.hardy_GL_C_tapp_q1$V3
table_het_per_scaf_GL_C$V5<-new.hardy_GL_C_tapp_q2$V3
table_het_per_scaf_GL_C$V6<-new.hardy_GL_C_tapp_q3$V3
table_het_per_scaf_GL_C$V7<-new.hardy_GL_C_tapp_max$V3
colnames(table_het_per_scaf_GL_C)<-c("scaffold", "window", "5%", "25%", "mean", "75%", "95%")
table_het_per_scaf_GL_C[table_het_per_scaf_GL_C == NA] <- NA
#plot(table_het_per_scaf_GL_C$mean)
table_het_per_scaf_GL_C_sub<-table_het_per_scaf_GL_C[complete.cases(table_het_per_scaf_GL_C),]
SNP<-c(1:(nrow(table_het_per_scaf_GL_C_sub)))
df_table_het_per_scaf_GL_C<-data.frame(SNP,table_het_per_scaf_GL_C_sub)
df_table_het_per_scaf_GL_C$scaffold <- as.numeric(df_table_het_per_scaf_GL_C$scaffold)
df_table_het_per_scaf_GL_C$window <- as.numeric(df_table_het_per_scaf_GL_C$window)
df_table_het_per_scaf_GL_C$mean <- as.numeric(df_table_het_per_scaf_GL_C$mean)
manhattan(df_table_het_per_scaf_GL_C, chr="scaffold", bp="window", p="mean", snp="SNP", logp=FALSE, ylab="Heterozygosity in windows of 100K SNPs")

hardy_GL_C<-read.table(text = gsub("/", "\t", readLines("LDfilt_hardy_GC.hwe")), header = T)
hardy_GL_C$hetero<-hardy_GL_C$HET/(hardy_GL_C$OBS.HOM1+hardy_GL_C$HET+hardy_GL_C$HOM2.)
hardy_GL_C$hetero_exp<-hardy_GL_C$HET.1/(hardy_GL_C$E.HOM1+hardy_GL_C$HET.1+hardy_GL_C$HOM2..1)
hardy_GL_C$CHR <- as.numeric(revalue(hardy_GL_C$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_GL_H<-read.table(text = gsub("/", "\t", readLines("LDfilt_hardy_GH.hwe")), header = T)
hardy_GL_H$hetero<-hardy_GL_H$HET/(hardy_GL_H$OBS.HOM1+hardy_GL_H$HET+hardy_GL_H$HOM2.)
hardy_GL_H$hetero_exp<-hardy_GL_H$HET.1/(hardy_GL_H$E.HOM1+hardy_GL_H$HET.1+hardy_GL_H$HOM2..1)
hardy_GL_H$CHR <- as.numeric(revalue(hardy_GL_H$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_IS_C<-read.table(text = gsub("/", "\t", readLines("LDfilt_hardy_IC.hwe")), header = T)
hardy_IS_C$hetero<-hardy_IS_C$HET/(hardy_IS_C$OBS.HOM1+hardy_IS_C$HET+hardy_IS_C$HOM2.)
hardy_IS_C$hetero_exp<-hardy_IS_C$HET.1/(hardy_IS_C$E.HOM1+hardy_IS_C$HET.1+hardy_IS_C$HOM2..1)
hardy_IS_C$CHR <- as.numeric(revalue(hardy_IS_C$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_IS_H<-read.table(text = gsub("/", "\t", readLines("LDfilt_hardy_IH.hwe")), header = T)
hardy_IS_H$hetero<-hardy_IS_H$HET/(hardy_IS_H$OBS.HOM1+hardy_IS_H$HET+hardy_IS_H$HOM2.)
hardy_IS_H$hetero_exp<-hardy_IS_H$HET.1/(hardy_IS_H$E.HOM1+hardy_IS_H$HET.1+hardy_IS_H$HOM2..1)
hardy_IS_H$CHR <- as.numeric(revalue(hardy_IS_H$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_NO_C<-read.table(text = gsub("/", "\t", readLines("LDfilt_hardy_NC.hwe")), header = T)
hardy_NO_C$hetero<-hardy_NO_C$HET/(hardy_NO_C$OBS.HOM1+hardy_NO_C$HET+hardy_NO_C$HOM2.)
hardy_NO_C$hetero_exp<-hardy_NO_C$HET.1/(hardy_NO_C$E.HOM1+hardy_NO_C$HET.1+hardy_NO_C$HOM2..1)
hardy_NO_C$CHR <- as.numeric(revalue(hardy_NO_C$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_NO_H<-read.table(text = gsub("/", "\t", readLines("LDfilt_hardy_NH.hwe")), header = T)
hardy_NO_H$hetero<-hardy_NO_H$HET/(hardy_NO_H$OBS.HOM1+hardy_NO_H$HET+hardy_NO_H$HOM2.)
hardy_NO_H$hetero_exp<-hardy_NO_H$HET.1/(hardy_NO_H$E.HOM1+hardy_NO_H$HET.1+hardy_NO_H$HOM2..1)
hardy_NO_H$CHR <- as.numeric(revalue(hardy_NO_H$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_DK_C<-read.table(text = gsub("/", "\t", readLines("LDfilt_hardy_DC.hwe")), header = T)
hardy_DK_C$hetero<-hardy_DK_C$HET/(hardy_DK_C$OBS.HOM1+hardy_DK_C$HET+hardy_DK_C$HOM2.)
hardy_DK_C$hetero_exp<-hardy_DK_C$HET.1/(hardy_DK_C$E.HOM1+hardy_DK_C$HET.1+hardy_DK_C$HOM2..1)
hardy_DK_C$CHR <- as.numeric(revalue(hardy_DK_C$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_DK_H<-read.table(text = gsub("/", "\t", readLines("LDfilt_hardy_DH.hwe")), header = T)
hardy_DK_H$hetero<-hardy_DK_H$HET/(hardy_DK_H$OBS.HOM1+hardy_DK_H$HET+hardy_DK_H$HOM2.)
hardy_DK_H$hetero_exp<-hardy_DK_H$HET.1/(hardy_DK_H$E.HOM1+hardy_DK_H$HET.1+hardy_DK_H$HOM2..1)
hardy_DK_H$CHR <- as.numeric(revalue(hardy_DK_H$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_EE_C<-read.table(text = gsub("/", "\t", readLines("LDfilt_hardy_EC.hwe")), header = T)
hardy_EE_C$hetero<-hardy_EE_C$HET/(hardy_EE_C$OBS.HOM1+hardy_EE_C$HET+hardy_EE_C$HOM2.)
hardy_EE_C$hetero_exp<-hardy_EE_C$HET.1/(hardy_EE_C$E.HOM1+hardy_EE_C$HET.1+hardy_EE_C$HOM2..1)
hardy_EE_C$CHR <- as.numeric(revalue(hardy_EE_C$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_TU_H<-read.table(text = gsub("/", "\t", readLines("LDfilt_hardy_TH.hwe")), header = T)
hardy_TU_H$hetero<-hardy_TU_H$HET/(hardy_TU_H$OBS.HOM1+hardy_TU_H$HET+hardy_TU_H$HOM2.)
hardy_TU_H$hetero_exp<-hardy_TU_H$HET.1/(hardy_TU_H$E.HOM1+hardy_TU_H$HET.1+hardy_TU_H$HOM2..1)
hardy_TU_H$CHR <- as.numeric(revalue(hardy_TU_H$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))


#pdf("Heterozygosity_mean_along_chromosome_10pops_030621.pdf", width=13, height=8)
par(mfrow=c(5,2))
par(mar=c(2.5,2.5,1.5,0.5))
boxplot(hetero~CHR,data=hardy_GL_C, ylab="", xlab="", cex.axis=1.5)
title("A", adj=0.02, line = 0.5)
boxplot(hetero~CHR,data=hardy_GL_H, ylab="", xlab="", cex.axis=1.5)
title("B", adj=0.02, line = 0.5)
boxplot(hetero~CHR,data=hardy_IS_C, ylab="", xlab="", cex.axis=1.5)
title("C", adj=0.02, line = 0.5)
boxplot(hetero~CHR,data=hardy_IS_H, ylab="", xlab="", cex.axis=1.5)
title("D", adj=0.02, line = 0.5)
boxplot(hetero~CHR,data=hardy_NO_C, ylab="", xlab="", cex.axis=1.5)
title("E", adj=0.02, line = 0.5)
boxplot(hetero~CHR,data=hardy_NO_H, ylab="", xlab="", cex.axis=1.5)
title("F", adj=0.02, line = 0.5)
boxplot(hetero~CHR,data=hardy_DK_C, ylab="", xlab="", cex.axis=1.5)
title("G", adj=0.02, line = 0.5)
boxplot(hetero~CHR,data=hardy_DK_H, ylab="", xlab="", cex.axis=1.5)
title("H", adj=0.02, line = 0.5)
boxplot(hetero~CHR,data=hardy_EE_C, ylab="", xlab="", cex.axis=1.5)
title("I", adj=0.02, line = 0.5)
boxplot(hetero~CHR,data=hardy_TU_H, ylab="", xlab="", cex.axis=1.5)
title("J", adj=0.02, line = 0.5)
dev.off()

cumsum_hardy_GL_C<-cumsum(sort(na.omit(hardy_GL_C$P_HWE)))
cumsum_hardy_GL_H<-cumsum(sort(na.omit(hardy_GL_H$P_HWE)))
cumsum_hardy_IS_C<-cumsum(sort(na.omit(hardy_IS_C$P_HWE)))
cumsum_hardy_IS_H<-cumsum(sort(na.omit(hardy_IS_H$P_HWE)))
cumsum_hardy_NO_C<-cumsum(sort(na.omit(hardy_NO_C$P_HWE)))
cumsum_hardy_NO_H<-cumsum(sort(na.omit(hardy_NO_H$P_HWE)))
cumsum_hardy_DK_C<-cumsum(sort(na.omit(hardy_DK_C$P_HWE)))
cumsum_hardy_DK_H<-cumsum(sort(na.omit(hardy_DK_H$P_HWE)))
cumsum_hardy_EE_C<-cumsum(sort(na.omit(hardy_EE_C$P_HWE)))
cumsum_hardy_TU_H<-cumsum(sort(na.omit(hardy_TU_H$P_HWE)))

mean(na.omit(hardy_GL_C$hetero))
mean(na.omit(hardy_GL_H$hetero))
mean(na.omit(hardy_IS_C$hetero))
mean(na.omit(hardy_IS_H$hetero))
mean(na.omit(hardy_NO_C$hetero))
mean(na.omit(hardy_NO_H$hetero))
mean(na.omit(hardy_DK_C$hetero))
mean(na.omit(hardy_DK_H$hetero))
mean(na.omit(hardy_EE_C$hetero))
mean(na.omit(hardy_TU_H$hetero))

var(na.omit(hardy_GL_C$hetero))
var(na.omit(hardy_GL_H$hetero))
var(na.omit(hardy_IS_C$hetero))
var(na.omit(hardy_IS_H$hetero))
var(na.omit(hardy_NO_C$hetero))
var(na.omit(hardy_NO_H$hetero))
var(na.omit(hardy_DK_C$hetero))
var(na.omit(hardy_DK_H$hetero))
var(na.omit(hardy_EE_C$hetero))
var(na.omit(hardy_TU_H$hetero))

mean(na.omit(hardy_GL_C$hetero_exp))
mean(na.omit(hardy_GL_H$hetero_exp))
mean(na.omit(hardy_IS_C$hetero_exp))
mean(na.omit(hardy_IS_H$hetero_exp))
mean(na.omit(hardy_NO_C$hetero_exp))
mean(na.omit(hardy_NO_H$hetero_exp))
mean(na.omit(hardy_DK_C$hetero_exp))
mean(na.omit(hardy_DK_H$hetero_exp))
mean(na.omit(hardy_EE_C$hetero_exp))
mean(na.omit(hardy_TU_H$hetero_exp))

var(na.omit(hardy_GL_C$hetero_exp))
var(na.omit(hardy_GL_H$hetero_exp))
var(na.omit(hardy_IS_C$hetero_exp))
var(na.omit(hardy_IS_H$hetero_exp))
var(na.omit(hardy_NO_C$hetero_exp))
var(na.omit(hardy_NO_H$hetero_exp))
var(na.omit(hardy_DK_C$hetero_exp))
var(na.omit(hardy_DK_H$hetero_exp))
var(na.omit(hardy_EE_C$hetero_exp))
var(na.omit(hardy_TU_H$hetero_exp))

