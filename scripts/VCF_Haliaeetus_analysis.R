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

#write.csv(OutLoc_annot,"OutLoc_annot.csv",quote = F,row.names = F)
#write.csv(OutAnnot,"OutAnnot.csv",quote = F,row.names = F)

OutLocAnnot<-merge(OutLoc_annot,OutAnnot,by.x="Annot_out",by.y="TR")
#write.csv(OutLocAnnot,"OutLocAnnot.csv",quote = F,row.names = F)
#OutLocAnnot<-read.csv("OutLocAnnot.csv")

OutLocAnnot<-read.csv("OutLocAnnot.csv")

Y<-read.vcfR("All.calls.FFF.recode.Annot_NRecip_TSL.vcf")

#write.vcf(Z, file = "All.calls.FFF.recode.Annot_NRecip_TSL_outliers.vcf", mask = FALSE, APPEND = FALSE)
###################################################
# Try the same as above but with no T. ventricosus
###################################################

G_nT<-G[,-c(68:77)] # Filter T.v.
SNPmat_nT<-(t(G_nT))
#G[is.na(G)]<-9

chromosome_n<-as.numeric(factor(chromosome))
pos_loc<-paste(chromosome,position,sep="_")
colnames(SNPmat_nT)<-pos_loc
G_nMT<-SNPmat_nT[ , apply(SNPmat_nT, 2, function(x) !any(is.na(x)))]
chrom_fil<-colnames(G_nMT)[colnames(G_nMT)%in%pos_loc]

#chrom_filt<-as.numeric(sub("_.*", "", chrom_fil))
#pos_filt<-as.numeric(sub(".*_", "", chrom_fil))
pop_nT<-pop[-c(68:77)]
my_fst_nT <- MakeDiploidFSTMat(G_nMT, locusNames = chrom_fil, popNames = pop_nT)

plot(my_fst_nT$He, my_fst_nT$FST, xlab="Heterozygosity", ylab="F_ST across all populations", main="Per locus F_ST vs Heterozygosity")

# If a locus has a much lower sample size compared to the rest, it could have a broader error distribution (and therefore incorrectly inferred as an outlier).
plot(my_fst_nT$FST, my_fst_nT$FSTNoCorr)
abline(0,1)

Ga_nT<-add_code256(big_copy(G_nMT,type="raw"),code=bigsnpr:::CODE_012)
ind.keep<-snp_clumping(G=Ga_nT,infos.chr=chrom_filt ,infos.pos=pos_filt)
m <- ncol(Ga_nT)
length(ind.keep) / m

out_trim <- OutFLANK(my_fst_nT[ind.keep,], NumberOfSamples=length(unique(pop_nT)), qthreshold = 0.01, Hmin = 0.1)
str(out_trim)

OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.01, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
## Zoom in on right tail
OutFLANKResultsPlotter(out_trim , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.01, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)

hist(out_trim$results$pvaluesRightTail)

P1 <- pOutlierFinderChiSqNoCorr(my_fst_nT, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)

(my_out <- P1$OutlierFlag==TRUE)
plot(P1$He, P1$FST, pch=19, col=rgb(0,0,0,0.1))
points(P1$He[my_out], P1$FST[my_out], col="blue")

hist(P1$pvaluesRightTail)

OutLoc<-P1[P1$OutlierFlag==TRUE,]
OutLoc<-OutLoc[!is.na(OutLoc$OutlierFlag),]

#######################################################
# PCA with only significant outliers
#######################################################

#Z<-Y[apply(Y@fix[,c(1,2)],1,paste,collapse="") %in% apply(cbind(as.character(OutLocAnnot$Annot_out),OutLocAnnot$pos_out),1,paste,collapse="")]

#write.vcf(Z, file = "All.calls.FFF.recode.Annot_NRecip_TSL_outliers.vcf", mask = FALSE, APPEND = FALSE)
#vcf.fn<-"~/Desktop/Data_Files/Exon_cap_seq-UseThisOne/All.calls.FFF.recode.Annot_NRecip_TSL_outliers.vcf"
#snpgdsVCF2GDS(vcf.fn, "All_calls_annot_outliers.gds",method = "copy.num.of.ref")

#snpgdsClose(AllCallAnnot)
AllCallAnnot <- snpgdsOpen("All_calls_annot_outliers.gds")

pcaC <- snpgdsPCA(AllCallAnnot,autosome.only=F)

pc.percent <- pcaC$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pcaC$sample.id,
                  EV1 = pcaC$eigenvect[,1],
                  EV2 = pcaC$eigenvect[,2],
                  EV3 = pcaC$eigenvect[,3],
                  EV4 = pcaC$eigenvect[,4],
                  EV5 = pcaC$eigenvect[,5],
                  EV6 = pcaC$eigenvect[,6],
                  stringsAsFactors = FALSE)

FiltSampNames<-read.csv("Names_filt.csv")

tab$sample.id<-FiltSampNames$Pop_ind

Per_exp<-head(round(pc.percent, 2))

P1_2 <- ggplot(tab, aes(EV1,EV2,color=FiltSampNames$Pop,label=sample.id)) +xlab(paste("PC1 (",Per_exp[1],"%)",sep="")) + ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + geom_point() + stat_ellipse(level=0.75)+labs(color="") + theme_classic()
P3_2 <- ggplot(tab, aes(EV3,EV2,color=FiltSampNames$Pop)) +xlab(paste("PC3 (",Per_exp[3],"%)",sep="")) + ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + geom_point() + stat_ellipse(level=0.75) +geom_text_repel(data=tab, aes(label=sample.id), size=3)+ guides(col=FALSE)  +theme_classic()
leg<-get_legend(P1_2+theme_classic()+theme(legend.margin=margin(t=5, r=0, b=0, l=0, unit="cm")))
P1_2_noL <- ggplot(tab, aes(EV1,EV2,color=FiltSampNames$Pop,label=sample.id)) +xlab(paste("PC1 (",Per_exp[1],"%)",sep="")) +ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + geom_point() + stat_ellipse(level=0.75)+labs(color="")+ theme_classic() + guides(col=FALSE)
plot_grid(leg,P1_2_noL, NULL, P3_2,ncol = 2, nrow = 2,rel_widths = c(1,3,3))



No_Tv<-pcaC$sample.id[-c(68:77)]
pcaC_noTv <- snpgdsPCA(AllCallAnnot,autosome.only=F,sample.id=No_Tv)

pc_noTv.percent <- pcaC_noTv$varprop*100
head(round(pc_noTv.percent, 2))

tab_noTv <- data.frame(sample.id = pcaC_noTv$sample.id,
                       EV1 = pcaC_noTv$eigenvect[,1],
                       EV2 = pcaC_noTv$eigenvect[,2],
                       EV3 = pcaC_noTv$eigenvect[,3],
                       EV4 = pcaC_noTv$eigenvect[,4],
                       EV5 = pcaC_noTv$eigenvect[,5],
                       EV6 = pcaC_noTv$eigenvect[,6],
                       stringsAsFactors = FALSE)
head(tab_noTv)
tail(tab_noTv)

P1_2 <- ggplot(tab_noTv, aes(EV1,EV2,color = FiltSampNames$Pop[-c(68:77)],label=sample.id)) + xlab(paste("PC1 (",Per_exp[1],"%)",sep="")) + ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + geom_point() + stat_ellipse(level=0.75)+labs(color="") + theme_classic()
P3_2 <- ggplot(tab_noTv, aes(EV3,EV2,color = FiltSampNames$Pop[-c(68:77)])) + xlab(paste("PC3 (",Per_exp[3],"%)",sep="")) + ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + geom_point() + stat_ellipse(level=0.75) + geom_text_repel(data=tab_noTv, aes(label=sample.id), size=3) + guides(col=FALSE) + theme_classic()
leg <- get_legend(P1_2 + theme_classic() + theme(legend.margin=margin(t=5, r=0, b=0, l=0, unit="cm")))
P1_2_noL <- ggplot(tab_noTv, aes(EV1,EV2,color=FiltSampNames$Pop[-c(68:77)],label=sample.id)) +xlab(paste("PC1 (",Per_exp[1],"%)",sep="")) +ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + geom_point() + stat_ellipse(level=0.75)+labs(color="")+ theme_classic() + guides(col=FALSE)
plot_grid(leg,P1_2_noL, NULL, P3_2,ncol = 2, nrow = 2,rel_widths = c(1,3,3))

#######################################################
# PCA with only non-significant outliers
#######################################################

#outlier SNPs (387 across 133 sequences) and non-outlier SNPs (2189 across 247 sequences)

#Z<-Y[which(!Y@fix[,1]%in%OutLocAnnot$Annot_out),]

#write.vcf(Z, file = "All.calls.FFF.recode.Annot_NRecip_TSL_NONoutliers.vcf", mask = FALSE, APPEND = FALSE)
#vcf.fn<-"All.calls.FFF.recode.Annot_NRecip_TSL_NONoutliers.vcf"
#snpgdsVCF2GDS(vcf.fn, "All_calls_annot_NONoutliers.gds",method = "copy.num.of.ref")

#snpgdsClose(AllCallAnnot)
AllCallAnnot <- snpgdsOpen("All_calls_annot_NONoutliers.gds")

pcaC <- snpgdsPCA(AllCallAnnot,autosome.only=F)

pc.percent <- pcaC$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pcaC$sample.id,
                  EV1 = pcaC$eigenvect[,1],
                  EV2 = pcaC$eigenvect[,2],
                  EV3 = pcaC$eigenvect[,3],
                  EV4 = pcaC$eigenvect[,4],
                  EV5 = pcaC$eigenvect[,5],
                  EV6 = pcaC$eigenvect[,6],
                  stringsAsFactors = FALSE)

FiltSampNames<-read.csv("Names_filt.csv")

tab$sample.id<-FiltSampNames$Pop_ind

Per_exp<-head(round(pc.percent, 2))

P1_2 <- ggplot(tab, aes(EV1,EV2,color=FiltSampNames$Pop,label=sample.id)) +xlab(paste("PC1 (",Per_exp[1],"%)",sep="")) + ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + geom_point() + stat_ellipse(level=0.75)+labs(color="") + theme_classic()
P3_2 <- ggplot(tab, aes(EV3,EV2,color=FiltSampNames$Pop)) +xlab(paste("PC3 (",Per_exp[3],"%)",sep="")) + ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + geom_point() + stat_ellipse(level=0.75) +geom_text_repel(data=tab, aes(label=sample.id), size=3)+ guides(col=FALSE)  +theme_classic()
leg<-get_legend(P1_2+theme_classic()+theme(legend.margin=margin(t=5, r=0, b=0, l=0, unit="cm")))
P1_2_noL <- ggplot(tab, aes(EV1,EV2,color=FiltSampNames$Pop,label=sample.id)) +xlab(paste("PC1 (",Per_exp[1],"%)",sep="")) +ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + geom_point() + stat_ellipse(level=0.75)+labs(color="")+ theme_classic() + guides(col=FALSE)
plot_grid(leg,P1_2_noL, NULL, P3_2,ncol = 2, nrow = 2,rel_widths = c(1,3,3))

################################################################################
Y<-read.vcfR("All.calls.FFF.recode.Annot_NRecip_TSL_outliers.vcf.gz")
mito <- read.csv("Tg_outlier_sequences/Mito.csv", header=F)
mito <- mito$V4
M<-X[which(X@fix[,1]%in%mito),]

write.vcf(M, file = "All.calls.FFF.recode.Annot_NRecip_TSL_mito.vcf", mask = FALSE, APPEND = FALSE)
vcf.fn<-"All.calls.FFF.recode.Annot_NRecip_TSL_mito.vcf"
snpgdsVCF2GDS(vcf.fn, "All_calls_annot_mito.gds",method = "copy.num.of.ref")

AllCallAnnot <- snpgdsOpen("All_calls_annot_mito.gds")
#snpgdsClose(AllCallAnnot)

pcaC <- snpgdsPCA(AllCallAnnot,autosome.only=F)

pc.percent <- pcaC$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pcaC$sample.id,
                  EV1 = pcaC$eigenvect[,1],
                  EV2 = pcaC$eigenvect[,2],
                  EV3 = pcaC$eigenvect[,3],
                  EV4 = pcaC$eigenvect[,4],
                  EV5 = pcaC$eigenvect[,5],
                  EV6 = pcaC$eigenvect[,6],
                  stringsAsFactors = FALSE)

FiltSampNames<-read.csv("Names_filt.csv")

tab$sample.id<-FiltSampNames$Pop_ind

Per_exp<-head(round(pc.percent, 2))

P1_2 <- ggplot(tab, aes(EV1,EV2,color=FiltSampNames$Pop,label=sample.id)) +xlab(paste("PC1 (",Per_exp[1],"%)",sep="")) + ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + geom_point() + stat_ellipse(level=0.75)+labs(color="") + theme_classic()
P3_2 <- ggplot(tab, aes(EV3,EV2,color=FiltSampNames$Pop)) +xlab(paste("PC3 (",Per_exp[3],"%)",sep="")) + ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + geom_point() + stat_ellipse(level=0.75) +geom_text_repel(data=tab, aes(label=sample.id), size=3)+ guides(col=FALSE)  +theme_classic()
leg<-get_legend(P1_2+theme_classic()+theme(legend.margin=margin(t=5, r=0, b=0, l=0, unit="cm")))
P1_2_noL <- ggplot(tab, aes(EV1,EV2,color=FiltSampNames$Pop,label=sample.id)) +xlab(paste("PC1 (",Per_exp[1],"%)",sep="")) +ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + geom_point() + stat_ellipse(level=0.75)+labs(color="")+ theme_classic() + guides(col=FALSE)
plot_grid(leg,P1_2_noL, NULL, P3_2,ncol = 2, nrow = 2,rel_widths = c(1,3,3))

#####################################################################################################

Z<-read.vcfR("All.calls.FFF.recode.Annot_NRecip_TSL_NONoutliers.vcf")

length(getCHROM(Y))
length(unique(getCHROM(Y)))
length(getCHROM(Z))
length(unique(getCHROM(Z)))

#refs<-list.files()[grep("_ref_reform.fasta",list.files())]

refs<-list.files()[grep("CNS.fa",list.files())]
out <- read.table("Outlier_fragment_names.txt")

for(i in refs){
  print(i)
  S <- read.fasta(i)
  Sf<-S[names(S)%in%out$V1]
  print(paste("Lengths of filtered",i,"=",length(Sf)))
  write.fasta(sequences=Sf, names=names(Sf), file.out=paste(substr(i,1,3),"_CNS_outliers.fas",sep=""))
}                         

write.fasta(sequences=Sf, names=names(Sf), file.out="Tg_ExoCap_orfs_outliers.fas")


##################################################################################
# Environmental data
##################################################################################

TP<-read.csv("SOCATv2021_TropicalPacfic.csv")
IO<-read.csv("SOCATv2021_Indian.csv",sep="\t")

Env<-merge(TP,IO,all = T)

