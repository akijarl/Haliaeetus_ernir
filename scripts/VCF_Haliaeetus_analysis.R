#####################################################
#Code for Haliaeetus LD fitlered data
#####################################################

setwd("/home/aki/Documents/Rannsóknir/Haliaeetus/VCF")
#setwd("/Users/akijarl/Documents/UI/rannsoknir/Haliaeetus/")

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

X<-read.vcfR("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.vcf")
#Y<-maf(X, element=1)
#y<-maf(X, element=2)

#y_o<-y[order(y[,"Frequency"],decreasing = T),]
#barplot(y_o[,"Frequency"], names.arg="", main="Minor allele frequency")

#Y_o<-Y[order(Y[,"Frequency"],decreasing = T),]
#barplot(Y_o[,"Frequency"], names.arg="", main="Major allele frequency")

#remove indels
X <- extract.indels(X)

queryMETA(X)
queryMETA(X, element = 'DP')

dp <- extract.gt(X, element = "DP", as.numeric=TRUE)
# barplot(dp, ylab = "Read depth (DP)")
# abline(h=min(dp),col="red")
# abline(h=summary(dp)[[2]],col="orange")
# abline(h=mean(dp),col="orange")
# abline(h=median(dp),col="cyan")
# abline(h=summary(dp)[[5]],col="orange")
# abline(h=max(dp),col="red")

summary(dp)
colMeans(dp)
summary(colMeans(dp))
barplot(colMeans(dp), ylab = "Mean read depth (DP)")
#Contemporary depth
summary(colMeans(dp[,c(1:25,26:37,43:53,54:56,57:68)]))

#historical depth
summary(colMeans(dp[,c(-1:-25,-26:-37,-43:-53,-54:-56,-57:-68)]))
#gt <- extract.gt(X)
#hets <- is_het(gt)

#Annotated Contig names in TS_minE_best
#Y<-X[which(X@fix[,1]%in%Annot_minE_best$TR),]
#Y<-X[which(X@fix[,1]%in%Annot_minE_best_TSL$TR),]
#Y<-X[which(X@fix[,1]%in%Annot_minE_best_NRecip_TSL$TR),]

#dp <- extract.info(Y, element = "DP", as.numeric=TRUE)
#barplot(dp, ylab = "Read depth (DP)")
#summary(dp)

geno <- extract.gt(X) # Character matrix Containing the genotypes
position <- getPOS(X) # Positions in bp
chromosome <- getCHROM(X) # Chromosome information
chromosome_n<-as.numeric(factor(chromosome))
pos_loc<-paste(chromosome_n,position,sep="_")

pos_loc_all<-pos_loc

tab<-read.table("../VCF/ernir_LD.eigenvec")
pop<-tab$V12

#test <- geno[!geno %in% c("0/0","0/1", "1/0","1/1", "0/2", "2/0","2/2","1/2", "2/1","2/3", "0/3", "1/3")]
#test[!is.na(test)]

#identify variants that have unusual genotypes (beyond 0/0, 0/1, 1/1)
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

#separate out the Contemporary individuals
geno_c<-geno[,c(1:25,26:37,43:53,54:56,57:68)]
colnames(geno_c)

sum(is.na(geno_c))/(nrow(geno_c)*ncol(geno_c))

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

table(as.vector(geno_sing[-nll]))

geno_sing<-geno_sing[-nll,]

sum(is.na(geno_sing))/(nrow(geno_sing)*ncol(geno_sing))

pos_loc_sing<-pos_loc[-mult]
pos_loc_sing<-pos_loc_sing[-nll]

G <- matrix(NA, nrow = nrow(geno_sing), ncol = ncol(geno_sing),dimnames = list(pos_loc_sing,colnames(geno_sing)) )

G[geno_sing %in% c("0/0")] <- 0
G[geno_sing %in% c("0/1", "1/0")] <- 1
G[geno_sing %in% c("1/1")] <- 2

# G[geno %in% c("0/2", "2/0", "2|0", "0|2")] <- 1
# G[geno %in% c("1/2", "2/1", "2|1", "1|2")] <- 4
# G[geno %in% c("2/2", "2|2")] <- 2
# G[geno %in% c("0/3", "3/0", "3|0", "0|3")] <- 6
# G[geno %in% c("1/3", "3/1", "3|1", "1|3")] <- 7
#G[geno %in% c("2/3", "3/2", "3|2", "2|3")] <- 8
#G[geno %in% c("3/3", "3|3")] <- 9

Gall <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno),dimnames = list(pos_loc,colnames(geno)) )

Gall[geno %in% c("0/0")] <- 0
Gall[geno %in% c("0/1", "1/0")] <- 1
Gall[geno %in% c("1/1")] <- 2
Gall[geno %in% c("0/2", "2/0")] <- 1
Gall[geno %in% c("1/2", "2/1")] <- 1
Gall[geno %in% c("2/2")] <- 2
Gall[geno %in% c("0/3", "3/0")] <- 1
Gall[geno %in% c("1/3", "3/1")] <- 1
Gall[geno %in% c("2/3", "3/2")] <- 1
Gall[geno %in% c("3/3")] <- 2

table(as.vector(Gall))

geno_sing_c<-geno_sing[,c(1:25,26:37,43:53,54:56,57:68)]
colnames(geno_sing_c)
row.names(geno_sing_c)
length(row.names(geno_sing_c))

#idenify sites that are invariant in the sample (i.e. all )
nll_c<-NULL
for(i in 1:nrow(geno_sing_c)){
  if(length(unique(geno_sing_c[i,]))==1){
    nll_c<-c(nll_c,i)
  }
}

table(as.vector(geno_sing_c[-nll_c,]))

geno_sing_c<-geno_sing_c[-nll_c,]

pos_loc_c<-pos_loc_sing[-nll_c]

Gc <- matrix(NA, nrow = nrow(geno_sing_c), ncol = ncol(geno_sing_c),dimnames = list(pos_loc_c,colnames(geno_sing_c)) )

Gc[geno_sing_c %in% c("0/0")] <- 0
Gc[geno_sing_c %in% c("0/1", "1/0")] <- 1
Gc[geno_sing_c %in% c("1/1")] <- 2

geno_sing_h<-geno_sing[,c(-1:-25,-26:-37,-43:-53,-54:-56,-57:-68)]

nll_h<-NULL
for(i in 1:nrow(geno_sing_h)){
  if(length(unique(geno_sing_h[i,]))==1){
    nll_h<-c(nll_h,i)
  }
}

geno_sing_h<-geno_sing_h[-nll_h,]

pos_loc_h<-pos_loc_sing[-nll_h]

Gh <- matrix(NA, nrow = nrow(geno_sing_h), ncol = ncol(geno_sing_h),dimnames = list(pos_loc_h,colnames(geno_sing_h)) )

Gh[geno_sing_h %in% c("0/0")] <- 0
Gh[geno_sing_h %in% c("0/1", "1/0")] <- 1
Gh[geno_sing_h %in% c("1/1")] <- 2

#test2 <- G[!G %in% c(0:9)]
#test2[!is.na(test2)]

#(sum(G %in% c(0:9))+sum(is.na(G[!G %in% c(0:9)]))) / (nrow(G)*ncol(G))

#table(as.vector(geno[nll,]))

#dat<-list(chromosome,position,G,pop,ind)

SNPmat<-(t(G))

SNPmatC<-(t(Gc))

SNPmatA<-(t(Gall))

#colnames(SNPmat)<-pos_loc
# G_nM<-SNPmat[ , apply(SNPmat, 2, function(x) !any(is.na(x)))]
# chrom_fil<-colnames(G_nM)[colnames(G_nM)%in%pos_loc]
# 
# chrom_filt<-as.numeric(sub("_.*", "", chrom_fil))
# pos_filt<-as.numeric(sub(".*_", "", chrom_fil))

G[is.na(G)]<-9

Gc[is.na(Gc)]<-9

Gh[is.na(Gh)]<-9

Gall[is.na(Gall)]<-9


summary(as.factor(Gc))/(nrow(Gc)*ncol(Gc))

missGc<-NULL
for(i in 1:nrow(Gc)){
  if(sum(Gc[i,]==9)/sum(summary(as.factor(Gc[i,])))>0.49){
    missGc<-c(missGc,i)
  }
}

missGcInd<-NULL
for(i in 1:ncol(Gc)){
  if(sum(Gc[,i]==9)/sum(summary(as.factor(Gc[,i])))>0.49){
    missGcInd<-c(missGcInd,i)
  }
}

colnames(Gc)[missGcInd]

#pos_loc_c<-pos_loc_c[-missGc]
#Gc<-Gc[-missGc,]
Gc<-Gc[,-missGcInd]


summary(as.factor(Gh))/(nrow(Gh)*ncol(Gh))

missGh<-NULL
for(i in 1:nrow(Gh)){
  if(sum(Gh[i,]==9)/sum(summary(as.factor(Gh[i,])))>0.49){
    missGh<-c(missGh,i)
  }
}

missGhInd<-NULL
for(i in 1:ncol(Gh)){
  if(sum(Gh[,i]==9)/sum(summary(as.factor(Gh[,i])))>0.49){
    missGhInd<-c(missGhInd,i)
  }
}

colnames(Gh)[missGhInd]

sum(Gh[,14]==9)/sum(summary(as.factor(Gh[,14])))
sum(Gh[,15]==9)/sum(summary(as.factor(Gh[,15])))
sum(Gh[,29]==9)/sum(summary(as.factor(Gh[,29])))

pos_loc_h<-pos_loc_h[-missGh]
Gh<-Gh[-missGh,]

#Gh<-Gh[,-c(14,15,29)]
Gh<-Gh[,-missGhInd]

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
pop_c <- pop_c[-missGcInd]
pop_h <- pop[c(-1:-25,-26:-37,-43:-53,-54:-56,-57:-68)]
pop_h <- pop_h[-missGhInd]
my_fst <- MakeDiploidFSTMat(SNPmat, locusNames = pos_loc, popNames = pop)

my_fst_c <- MakeDiploidFSTMat(t(Gc), locusNames = pos_loc_c, popNames = pop_c)

my_fst_h <- MakeDiploidFSTMat(t(Gh), locusNames = pos_loc_h, popNames = pop_h)

my_fst_all <- MakeDiploidFSTMat(t(Gall), locusNames = pos_loc_all, popNames = pop)

plot(my_fst$He, my_fst$FST, xlab="Heterozygosity", ylab=expression(paste("F"[ST], " across all populations", sep="")), main=expression(paste("Per locus F"[ST], " vs Heterozygosity", sep="")))

plot(my_fst_c$He, my_fst_c$FST, xlab="Heterozygosity", ylab=expression(paste("F"[ST], " across all populations", sep="")), main=expression(paste("Per locus F"[ST], " vs Heterozygosity", sep="")))

# If a locus has a much lower sample size compared to the rest, it could have a broader error distribution (and therefore incorrectly inferred as an outlier).
plot(my_fst$FST, my_fst$FSTNoCorr)
abline(a=0,b=1,col="red")

plot(my_fst_c$FST, my_fst_c$FSTNoCorr)
abline(a=0,b=1,col="red")

#### SNP trimming from OutFLANK vignette ####
# identify a quasi-independent set of SNPs to calculate FSTbar and df. A common way of obtaining these SNPs is to thin for linkage disequilibrium (SNP thinning), which typically moves along a genome in a sliding window and thins SNPs based on linkage disequilibrium with each other. 
#This may be based on a combination of 
#(i) “pruning,” which sequentially scans the genome and performs pairwise thinning based on a given threshold of correlation, 
#(ii) “clumping,” which may incorporate some information about the importance of SNPs based on summary statistics, and 
#(iii) removing SNPs in long-range LD regions (Prive et al. 2017).

# Ga<-add_code256(big_copy(G_nM,type="raw"),code=bigsnpr:::CODE_012)
# ind.keep<-snp_clumping(G=Ga,infos.chr=chrom_filt ,infos.pos=pos_filt)
# m <- ncol(Ga)
# length(ind.keep) / m

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
sum(my_out,na.rm = T)
plot(P1$He, P1$FST, pch=19, col=rgb(0,0,0,0.1),xlab="Heterozygosity", ylab=expression(paste("F"[ST], " across all populations", sep="")), main=expression(paste("Per locus F"[ST], " vs Heterozygosity", sep="")))
points(P1$He[my_out], P1$FST[my_out], col="blue")

hist(P1$pvaluesRightTail)

plot(as.factor(P1[!is.na(P1$He),]$LocusName[P1[!is.na(P1$He),]$He>0.1]), P1[!is.na(P1$He),]$FST[P1[!is.na(P1$He),]$He>0.1],
     xlab="Position", ylab="FST", col=rgb(0,0,0,0.2))
points(P1$LocusName[my_out], P1$FST[my_out], col="magenta", pch=20)  

OutLoc<-P1[P1$OutlierFlag==TRUE,]
OutLoc<-OutLoc[!is.na(OutLoc$OutlierFlag),]

#Contemporary
ind.keep_c = 1:nrow(Gc)
out_trim_c <- OutFLANK(my_fst_c[ind.keep_c,], NumberOfSamples=length(unique(pop_c)), qthreshold = 0.001, Hmin = 0.1)
str(out_trim_c)
head(out_trim_c$results)
summary(out_trim_c$results$OutlierFlag)
summary(out_trim_c$results$pvaluesRightTail)

OutFLANKResultsPlotter(out_trim_c, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.01, binwidth = 0.01, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
## Zoom in on right tail
OutFLANKResultsPlotter(out_trim_c , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.01, binwidth = 0.01, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)

hist(out_trim_c$results$pvaluesRightTail)

P1_c <- pOutlierFinderChiSqNoCorr(my_fst_c, Fstbar = out_trim_c$FSTNoCorrbar, dfInferred = out_trim_c$dfInferred, qthreshold = 0.05, Hmin=0.1)

(my_out_c <- P1_c$OutlierFlag==TRUE)
sum(my_out_c,na.rm = T)
plot(P1_c$He, P1_c$FST, pch=19, col=rgb(0,0,0,0.1),xlab="Heterozygosity", ylab=expression(paste("F"[ST], " across all populations", sep="")), main=expression(paste("Per locus F"[ST], " vs Heterozygosity", sep="")))
points(P1_c$He[my_out], P1_c$FST[my_out], col="blue")

hist(P1_c$pvaluesRightTail)

plot(as.factor(P1_c[!is.na(P1_c$He),]$LocusName[P1_c[!is.na(P1_c$He),]$He>0.1]), P1_c[!is.na(P1_c$He),]$FST[P1_c[!is.na(P1_c$He),]$He>0.1],
     xlab="Position", ylab="FST", col=rgb(0,0,0,0.2))
points(P1_c$LocusName[my_out_c], P1_c$FST[my_out_c], col="magenta", pch=20)  

OutLoc_c<-P1_c[P1_c$OutlierFlag==TRUE,]
OutLoc_c<-OutLoc_c[!is.na(OutLoc_c$OutlierFlag),]

#historic
ind.keep_h<-1:nrow(Gh)

out_trim_h <- OutFLANK(my_fst_h[ind.keep_h,], NumberOfSamples=length(unique(pop_h)), qthreshold = 0.001, Hmin = 0.1)
str(out_trim_h)
head(out_trim_h$results)
summary(out_trim_h$results$OutlierFlag)
summary(out_trim_h$results$pvaluesRightTail)

OutFLANKResultsPlotter(out_trim_h, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.01, binwidth = 0.01, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
## Zoom in on right tail
OutFLANKResultsPlotter(out_trim_h , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.01, binwidth = 0.01, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)

hist(out_trim_h$results$pvaluesRightTail)

P1_h <- pOutlierFinderChiSqNoCorr(my_fst_h, Fstbar = out_trim_h$FSTNoCorrbar, dfInferred = out_trim_h$dfInferred, qthreshold = 0.05, Hmin=0.1)

P1_h_FST1<-P1_h[P1_h$FST==1,]
P1_h_FST1<-P1_h_FST1[!is.na(P1_h_FST1$FST),]


(my_out_h <- P1_h$OutlierFlag==TRUE)
sum(my_out_h,na.rm = T)
plot(P1_h$He, P1_h$FST, pch=19, col=rgb(0,0,0,0.1),xlab="Heterozygosity", ylab=expression(paste("F"[ST], " across all populations", sep="")), main=expression(paste("Per locus F"[ST], " vs Heterozygosity", sep="")))
points(P1_h$He[my_out], P1_h$FST[my_out], col="blue")

hist(P1_h$pvaluesRightTail)

plot(as.factor(P1_h[!is.na(P1_h$He),]$LocusName[P1_h[!is.na(P1_h$He),]$He>0.1]), P1_h[!is.na(P1_h$He),]$FST[P1_h[!is.na(P1_h$He),]$He>0.1],
     xlab="Position", ylab="FST", col=rgb(0,0,0,0.2))
points(P1_h$LocusName[my_out_h], P1_h$FST[my_out_h], col="magenta", pch=20)  

OutLoc_h<-P1_h[P1_h$OutlierFlag==TRUE,]
OutLoc_h<-OutLoc_h[!is.na(OutLoc_h$OutlierFlag),]


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

#setwd("/home/aki/Documents/Rannsóknir/Haliaeetus/VCF/")
setwd("/home/aki/Documents/Rannsóknir/Haliaeetus/Haliaeetus_ernir/")

require(tidyverse)
require(data.table)
require(vcfR)

hwe <- fread("LDfilt_hardy.hwe")

het <- fread("../Haliaeetus_ernir/LDfilt_HET_all.het")

#X<-read.vcfR("../VCF/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.vcf")
all_pops<-colnames(X@gt)[-1]

all_pops[!all_pops%in%het$INDV]

Pop <- c(rep("IS_C",25),rep("NO_C",12),rep("DK_H",5),rep("DK_C",11),rep("EE_C",3),rep("GL_C",12),rep("GL_H",8),rep("IS_H",2),rep("NO_H",13),"TU_H")

het$POP<-Pop

head(het)
het$O.HET<-het$N_SITES-het$`O(HOM)`
het$E.HET<-het$N_SITES-het$`E(HOM)`
het$hetO_frac<-het$O.HET/het$N_SITES
het$hetE_frac<-het$E.HET/het$N_SITES
het$POP<-factor(het$POP, levels = c("GL_C","GL_H","IS_C","IS_H","NO_C","NO_H","DK_C","DK_H","EE_C","TU_H"))
het[order(het$O.HET),]
het[order(het$hetO_frac),]

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
  #annotate(geom= "text", x=seq_len(unique(het_data_obs_exp_for_ggplot_true$POP)), y=10, label=het_data_obs_exp_for_ggplot_true$obs_or_exp) +
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

require(ggpubr)
het_gg_2figures_090421<-
  ggarrange(het_gg_observed_090421, het_gg_expected_090421,
            labels = c("A", "B"),
            ncol = 2, nrow = 1, common.legend = TRUE)

het[, .(meanFIS=mean(F)), by=POP]
het[, .(meanFIS=sd(F)), by=POP]

het[, .(meanFIS=mean(F))]
het[, .(meanFIS=sd(F))]

het[, .(meanE.HET=mean(E.HET)), by=POP]
het[, .(meanE.HET=sd(E.HET)), by=POP]

het[, .(meanO.HET=mean(O.HET)), by=POP]
het[, .(meanO.HET=sd(O.HET)), by=POP]

het[, .(meanHETefrac=mean(hetE_frac)), by=POP]
het[, .(meanHETefrac=sd(hetE_frac)), by=POP]

het[, .(meanHETofrac=mean(hetO_frac)), by=POP]
het[, .(meanHETofrac=sd(hetO_frac)), by=POP]

het$hetO_frac

het_data_obs_exp_for_ggplot_true
het[, .(meanE.HET=mean(E.HET)), by=POP]

wilcox.test(het$hetO_frac[het$POP=="GL_C"], het$hetO_frac[het$POP=="GL_H"], exact = T)
wilcox.test(het$hetO_frac[het$POP=="IS_C"], het$hetO_frac[het$POP=="IS_H"], exact = T)
wilcox.test(het$hetO_frac[het$POP=="NO_C"], het$hetO_frac[het$POP=="NO_H"], exact = T)
wilcox.test(het$hetO_frac[het$POP=="DK_C"], het$hetO_frac[het$POP=="DK_H"], exact = T)

wilcox.test(het$hetO_frac[het$POP=="NO_C" |het$POP=="NO_H" |
                                het$POP=="DK_C" |het$POP=="NO_H" | 
                                het$POP=="EE_C"| het$POP=="TU"], 
            het$hetO_frac[het$POP=="GL_C" | het$POP=="GL_H" |
                                het$POP=="IS_C" | het$POP=="IS_H"
            ], exact = T)

## IBS
ernir<-snpgdsOpen("../VCF/ernir.gds")
#snpgdsClose(ernir)
ibs <- snpgdsIBS(ernir, autosome.only = F, remove.monosnp = F)

IBS<-data.frame(ibs$ibs)
colnames(IBS)<-all_pops
rownames(IBS)<-all_pops

for(i in 1:ncol(IBS)){
    IBS[which(IBS[,i]==1.0000000),i]<-NA
}
#mean(as.numeric(as.dist(1-IBS[Pop=="GL_C",Pop=="GL_C"])))
IBS_dist<-1-IBS

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

mean(colMeans(IBS_dist[Pop=="GL_C",Pop=="GL_C"],na.rm = T))
mean(colMeans(IBS_dist[Pop=="GL_H",Pop=="GL_H"],na.rm = T))
mean(colMeans(IBS_dist[Pop=="IS_C",Pop=="IS_C"],na.rm = T))
mean(colMeans(IBS_dist[Pop=="IS_H",Pop=="IS_H"],na.rm = T))
mean(colMeans(IBS_dist[Pop=="NO_C",Pop=="NO_C"],na.rm = T))
mean(colMeans(IBS_dist[Pop=="NO_H",Pop=="NO_H"],na.rm = T))
mean(colMeans(IBS_dist[Pop=="DK_C",Pop=="DK_C"],na.rm = T))
mean(colMeans(IBS_dist[Pop=="DK_H",Pop=="DK_H"],na.rm = T))
mean(colMeans(IBS_dist[Pop=="EE_C",Pop=="EE_C"],na.rm = T))
mean(colMeans(IBS_dist[Pop=="TU_H",Pop=="TU_H"],na.rm = T))

mean(as.numeric(as.dist(IBS_dist)))
sd(as.numeric(as.dist(IBS_dist)))

sd(as.numeric(as.dist(IBS_dist[Pop=="GL_C",Pop=="GL_C"])))
sd(as.numeric(as.dist(IBS_dist[Pop=="GL_H",Pop=="GL_H"])))
sd(as.numeric(as.dist(IBS_dist[Pop=="IS_C",Pop=="IS_C"])))
sd(as.numeric(as.dist(IBS_dist[Pop=="IS_H",Pop=="IS_H"])))
sd(as.numeric(as.dist(IBS_dist[Pop=="NO_C",Pop=="NO_C"])))
sd(as.numeric(as.dist(IBS_dist[Pop=="NO_H",Pop=="NO_H"])))
sd(as.numeric(as.dist(IBS_dist[Pop=="DK_C",Pop=="DK_C"])))
sd(as.numeric(as.dist(IBS_dist[Pop=="DK_H",Pop=="DK_H"])))
sd(as.numeric(as.dist(IBS_dist[Pop=="EE_C",Pop=="EE_C"])))
sd(as.numeric(as.dist(IBS_dist[Pop=="TU_H",Pop=="TU_H"])))

wilcox.test(as.numeric(as.dist(IBS_dist[Pop=="GL_C",Pop=="GL_C"])),as.numeric(as.dist(IBS_dist[Pop=="GL_H",Pop=="GL_H"])),alternative = "greater")
wilcox.test(as.numeric(as.dist(IBS_dist[Pop=="GL_H",Pop=="GL_H"])),as.numeric(as.dist(IBS_dist[Pop=="GL_C",Pop=="GL_C"])),alternative = "greater")
wilcox.test(as.numeric(as.dist(IBS_dist[Pop=="IS_C",Pop=="IS_C"])),as.numeric(as.dist(IBS_dist[Pop=="IS_H",Pop=="IS_H"])),alternative = "greater")
wilcox.test(as.numeric(as.dist(IBS_dist[Pop=="IS_H",Pop=="IS_H"])),as.numeric(as.dist(IBS_dist[Pop=="IS_C",Pop=="IS_C"])),alternative = "greater")
wilcox.test(as.numeric(as.dist(IBS_dist[Pop=="NO_C",Pop=="NO_C"])),as.numeric(as.dist(IBS_dist[Pop=="NO_H",Pop=="NO_H"])),alternative = "greater")
wilcox.test(as.numeric(as.dist(IBS_dist[Pop=="NO_H",Pop=="NO_H"])),as.numeric(as.dist(IBS_dist[Pop=="NO_C",Pop=="NO_C"])),alternative = "greater")
wilcox.test(as.numeric(as.dist(IBS_dist[Pop=="DK_C",Pop=="DK_C"])),as.numeric(as.dist(IBS_dist[Pop=="DK_H",Pop=="DK_H"])),alternative = "greater")
wilcox.test(as.numeric(as.dist(IBS_dist[Pop=="DK_H",Pop=="DK_H"])),as.numeric(as.dist(IBS_dist[Pop=="DK_C",Pop=="DK_C"])),alternative = "greater")

require(plyr)
require(qqman)

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

hardy_all<-read.table(text = gsub("/", "\t", readLines("LDfilt_hardy.hwe")), header = T)
hardy_all$hetero<-hardy_all$HET/(hardy_all$OBS.HOM1+hardy_all$HET+hardy_all$HOM2.)
hardy_all$hetero_exp<-hardy_all$HET.1/(hardy_all$E.HOM1+hardy_all$HET.1+hardy_all$HOM2..1)
hardy_all$CHR <- as.numeric(revalue(hardy_all$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))


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

mean(na.omit(hardy_all$hetero))

sd(na.omit(hardy_GL_C$hetero))
sd(na.omit(hardy_GL_H$hetero))
sd(na.omit(hardy_IS_C$hetero))
sd(na.omit(hardy_IS_H$hetero))
sd(na.omit(hardy_NO_C$hetero))
sd(na.omit(hardy_NO_H$hetero))
sd(na.omit(hardy_DK_C$hetero))
sd(na.omit(hardy_DK_H$hetero))
sd(na.omit(hardy_EE_C$hetero))
sd(na.omit(hardy_TU_H$hetero))

sd(na.omit(hardy_all$hetero))

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

mean(na.omit(hardy_all$hetero_exp))

sd(na.omit(hardy_GL_C$hetero_exp))
sd(na.omit(hardy_GL_H$hetero_exp))
sd(na.omit(hardy_IS_C$hetero_exp))
sd(na.omit(hardy_IS_H$hetero_exp))
sd(na.omit(hardy_NO_C$hetero_exp))
sd(na.omit(hardy_NO_H$hetero_exp))
sd(na.omit(hardy_DK_C$hetero_exp))
sd(na.omit(hardy_DK_H$hetero_exp))
sd(na.omit(hardy_EE_C$hetero_exp))
sd(na.omit(hardy_TU_H$hetero_exp))

sd(na.omit(hardy_all$hetero_exp))

###ROH
ROH_92_relaxed_310321<-read.table("../Haliaeetus_ernir/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_LD_prune0.5_w134.indiv", header = T)
# ROH_92_relaxed_130921_plot<-
ggplot(data=ROH_92_relaxed_310321,
       aes(x=KB/1000, y=NSEG, color=Pop,  size=Pop, fill=Pop, shape=Pop)) + # for shape insert shape=Pop
  scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9,10),name="Pop_Time") + # for shape uncomments 
  geom_point() +
  scale_size_manual(values=c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5), element_blank(), name="") +
  scale_color_manual(values = c("#1F78B4", "#A6CEE3", "#6A3D9A", "#33A02C", "#B2DF8A", "#E31A1C", "#FB9A99", "#FF7F00", "#FDBF6F", "#CAB2D6" ), name="Pop_Time") +
  scale_fill_manual(values = c("#1F78B4", "#A6CEE3", "#6A3D9A", "#33A02C", "#B2DF8A", "#E31A1C", "#FB9A99", "#FF7F00", "#FDBF6F", "#CAB2D6" ),  name="Pop_Time") +
  guides(size = F, shape=F) + 
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  labs(color  = "Pop_Time", shape = "Pop_Time")+
  guides(color = guide_legend(override.aes = list(size = 5)),shape = guide_legend(override.aes = list(size = 5))) +
  #ylim(0,500) +
  labs(x ="Length of ROH (Mb)", y = "Number of ROH segments") +
  #geom_hline(yintercept = 0.000147102) +
  #geom_hline(yintercept = 0.0001613868, linetype='dashed') +
  #geom_hline(yintercept = 7.567803e-05, linetype='dotted') +
  theme(
    #axis.text.x = element_blank(),
    #axis.title.x = element_blank(),
    #axis.ticks.x=element_blank(),
    axis.title.x = element_text(size = 16),
    #axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    text = element_text(size = 14))+
  theme_classic()

autSize <- 1103368343 #length of the autosome

ROH_92_relaxed_310321$FROH<-(1000*ROH_92_relaxed_310321$KB)/autSize
ROH_92_relaxed_310321$POP<-Pop

ROH<-aggregate(ROH_92_relaxed_310321$FROH, list(ROH_92_relaxed_310321$POP), FUN=mean)
miss <- c("BB38","DKH2","DKH14","GLH1","GLH15","ISH4","ISH6","K9","K8","K6","K2","K4","NOH6","NOH7")
ROH_78<-ROH_92_relaxed_310321[!ROH_92_relaxed_310321$FID%in%miss,]
ROH_78<-aggregate(ROH_78$FROH, list(ROH_78$POP), FUN=mean)

colnames(ROH)<-c("Group","F_ROH")
colnames(ROH_78)<-c("Group","F_ROH")

FH <- data.frame(Group=c("GL_C","GL_H","IS_C","IS_H","NO_C","NO_H","DK_C","DK_H","EE_C","TU_H"),F_H=c(0.362,0.386,0.456,0.4096,-0.2204,-0.206,-0.0889,-0.227,-0.186,-0.139))
F_comp<-merge(ROH,FH,by = "Group")

plot(F_comp$F_ROH,F_comp$F_H)

F_comp_filt<-merge(ROH_78,FH[!FH$Group=="IS_H",],by = "Group")

F_comp<-data.frame("FH"=het$F,"FROH"=ROH_92_relaxed_310321$FROH)

ggplot(data=F_comp)+
  geom_point(aes(FROH,FH))+
  theme_classic()

ggplot()+
  geom_point(aes(het$missingness,(1000*ROH_92_relaxed_310321$KB)))+
  theme_classic()

cor.test(F_comp$FH,F_comp$FROH)

ROH_92_relaxed_310321$FROH


setwd("/home/aki/Documents/Rannsóknir/Haaliaeetus/VCF/")
IS_GL <- read.table("Iceland_Greenland_c_c_comp.weir.fst",header = T)
IS_NO <- read.table("Iceland_Norway_c_c_comp.weir.fst",header = T)
IS_DK <- read.table("Iceland_Denmark_c_c_comp.weir.fst",header = T)
IS_EE <- read.table("Iceland_Estonia_c_c_comp.weir.fst",header = T)
GL_NO <- read.table("Greenland_Norway_c_c_comp.weir.fst",header = T)
GL_DK <- read.table("Greenland_Denmark_c_c_comp.weir.fst",header = T)
GL_EE <- read.table("Greenland_Estonia_c_c_comp.weir.fst",header = T)
NO_DK <- read.table("Norway_Denmark_c_c_comp.weir.fst",header = T)
NO_EE <- read.table("Norway_Estonia_c_c_comp.weir.fst",header = T)
DK_EE <- read.table("Denmark_Estonia_c_c_comp.weir.fst",header = T)

ISh_GLh <- read.table("Iceland_Greenland_h_h_comp.weir.fst",header = T)
ISh_NOh <- read.table("Iceland_Norway_h_h_comp.weir.fst",header = T)
ISh_DKh <- read.table("Iceland_Denmark_h_h_comp.weir.fst",header = T)
ISh_TUh <- read.table("Iceland_Turkey_h_h_comp.weir.fst",header = T)
GLh_NOh <- read.table("Greenland_Norway_h_h_comp.weir.fst",header = T)
GLh_DKh <- read.table("Greenland_Denmark_h_h_comp.weir.fst",header = T)
GLh_TUh <- read.table("Greenland_Turkey_h_h_comp.weir.fst",header = T)
NOh_DKh <- read.table("Norway_Denmark_h_h_comp.weir.fst",header = T)
NOh_TUh <- read.table("Norway_Turkey_h_h_comp.weir.fst",header = T)
DKh_TUh <- read.table("Denmark_Turkey_h_h_comp.weir.fst",header = T)

IS_ch <- read.table("Iceland_c_h_comp.weir.fst",header = T)
GL_ch <- read.table("Greenland_c_h_comp.weir.fst",header = T)
NO_ch <- read.table("Norway_c_h_comp.weir.fst",header = T)
DK_ch <- read.table("Denmark_c_h_comp.weir.fst",header = T)

mean(IS_GL$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(IS_NO$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(IS_DK$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(IS_EE$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(GL_NO$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(GL_DK$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(GL_EE$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(NO_DK$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(NO_EE$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(DK_EE$WEIR_AND_COCKERHAM_FST,na.rm = T)

mean(ISh_GLh$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(ISh_NOh$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(ISh_DKh$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(ISh_TUh$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(GLh_NOh$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(GLh_DKh$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(GLh_TUh$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(NOh_DKh$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(NOh_TUh$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(DKh_TUh$WEIR_AND_COCKERHAM_FST,na.rm = T)

mean(IS_ch$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(GL_ch$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(NO_ch$WEIR_AND_COCKERHAM_FST,na.rm = T)
mean(DK_ch$WEIR_AND_COCKERHAM_FST,na.rm = T)

IS_ch$CHROM <- revalue(IS_ch$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
GL_ch$CHROM <- revalue(GL_ch$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
NO_ch$CHROM <- revalue(NO_ch$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
DK_ch$CHROM <- revalue(DK_ch$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))

IS_GL$CHROM <- revalue(IS_GL$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))

IS_ch$CHROM <- as.numeric(IS_ch$CHROM)
GL_ch$CHROM <- as.numeric(GL_ch$CHROM)
NO_ch$CHROM <- as.numeric(NO_ch$CHROM)
DK_ch$CHROM <- as.numeric(DK_ch$CHROM)

IS_GL$CHROM <- as.numeric(IS_GL$CHROM)

IS_ch$SNP<-seq.int(nrow(IS_ch))
GL_ch$SNP<-seq.int(nrow(GL_ch))
NO_ch$SNP<-seq.int(nrow(NO_ch))
DK_ch$SNP<-seq.int(nrow(DK_ch))

IS_GL$SNP<-seq.int(nrow(IS_GL))

IS_ch<-IS_ch[!is.na(IS_ch$WEIR_AND_COCKERHAM_FST),]
GL_ch<-GL_ch[!is.na(GL_ch$WEIR_AND_COCKERHAM_FST),]
NO_ch<-NO_ch[!is.na(NO_ch$WEIR_AND_COCKERHAM_FST),]
DK_ch<-DK_ch[!is.na(DK_ch$WEIR_AND_COCKERHAM_FST),]

IS_GL<-IS_GL[!is.na(IS_GL$WEIR_AND_COCKERHAM_FST),]

manhattan(IS_ch, chr="CHROM", bp="POS", p="WEIR_AND_COCKERHAM_FST", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5, ylim=c(0.0,1.1))
manhattan(GL_ch, chr="CHROM", bp="POS", p="WEIR_AND_COCKERHAM_FST", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5, ylim=c(0.0,1.1))
manhattan(NO_ch, chr="CHROM", bp="POS", p="WEIR_AND_COCKERHAM_FST", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5, ylim=c(0.0,1.1))
manhattan(DK_ch, chr="CHROM", bp="POS", p="WEIR_AND_COCKERHAM_FST", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5, ylim=c(0.0,1.1))
manhattan(IS_GL, chr="CHROM", bp="POS", p="WEIR_AND_COCKERHAM_FST", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5, ylim=c(0.0,1.1))




setwd("/home/aki/Documents/Rannsóknir/Haaliaeetus/VCF/FST_full/")

IS_ch <- read.table("fst_weir_ISM_ISH_110321.weir.fst",header = T)
IS_ch$CHROM <- revalue(IS_ch$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
IS_ch$CHROM <- as.numeric(IS_ch$CHROM)
IS_ch$SNP<-seq.int(nrow(IS_ch))
IS_ch<-IS_ch[!is.na(IS_ch$WEIR_AND_COCKERHAM_FST),]
manhattan(IS_ch, chr="CHROM", bp="POS", p="WEIR_AND_COCKERHAM_FST", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5, ylim=c(0.0,1.1))

GL_ch <- read.table("fst_weir_GLH_GLM_110321.weir.fst",header = T)
GL_ch$CHROM <- revalue(GL_ch$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
GL_ch$CHROM <- as.numeric(GL_ch$CHROM)
GL_ch$SNP<-seq.int(nrow(GL_ch))
GL_ch<-GL_ch[!is.na(GL_ch$WEIR_AND_COCKERHAM_FST),]
manhattan(GL_ch, chr="CHROM", bp="POS", p="WEIR_AND_COCKERHAM_FST", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5, ylim=c(0.0,1.1))

NO_ch <- read.table("fst_weir_NOM_NOH_110321.weir.fst",header = T)
NO_ch$CHROM <- revalue(NO_ch$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
NO_ch$CHROM <- as.numeric(NO_ch$CHROM)
NO_ch$SNP<-seq.int(nrow(NO_ch))
NO_ch<-NO_ch[!is.na(NO_ch$WEIR_AND_COCKERHAM_FST),]
manhattan(NO_ch, chr="CHROM", bp="POS", p="WEIR_AND_COCKERHAM_FST", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5, ylim=c(0.0,1.1))

DK_ch <- read.table("fst_weir_DKM_DKH_110321.weir.fst",header = T)
DK_ch$CHROM <- revalue(DK_ch$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
DK_ch$CHROM <- as.numeric(DK_ch$CHROM)
DK_ch$SNP<-seq.int(nrow(DK_ch))
DK_ch<-DK_ch[!is.na(DK_ch$WEIR_AND_COCKERHAM_FST),]
manhattan(DK_ch, chr="CHROM", bp="POS", p="WEIR_AND_COCKERHAM_FST", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5, ylim=c(0.0,1.1))


### All in one 3.2*10-9 - 100K years
setwd("E:/Research_AJL/Haliaeetus/Haliaeetus_ernir/stairway/dblcheck")
setwd("/home/aki/Documents/Rannsóknir/Haliaeetus/Haliaeetus_ernir/stairway/dblcheck")

GL_C_2.3_109_ang_auto <- read.table("GL_C_angsd_120721_10boot_100miosites_T15.6_mut3.29_10-9.final.summary", head=T)
IS_C_2.3_109_ang_auto <- read.table("IS_C_angsd_120721_10boot_100miosites_T15.6_mut3.29_10-9.final.summary", head=T)
NO_C_2.3_109_ang_auto <- read.table("NO_C_11_angsd_120721_10boot_100miosites_T15.6_mut3.29_10-9.final.summary", head=T)
DK_C_2.3_109_ang_auto <- read.table("DK_C_angsd_120721_10boot_100miosites_T15.6_mut3.29_10-9.final.summary", head=T)

GL_C_2.3_109_ang_auto$label <- "Greenland"
IS_C_2.3_109_ang_auto$label <- "Iceland"
NO_C_2.3_109_ang_auto$label <- "Norway"
DK_C_2.3_109_ang_auto$label <- "Denmark"

Con_SP <- rbind(GL_C_2.3_109_ang_auto,IS_C_2.3_109_ang_auto,NO_C_2.3_109_ang_auto,DK_C_2.3_109_ang_auto)

options(scipen=5)

#pdf("plot_stairway_4pops_raw_95_150721_T15.6_mut2.3_10-9_100Kyears_angsd_auto.pdf", height = 8, width = 12)
par(mar=c(6,4.5,0.5,0.5))
plot(GL_C_2.3_109_ang_auto$year, GL_C_2.3_109_ang_auto$Ne_median/1000, col="green", type="l", xlim=c(0,1600000), ylim=c(0,1000), cex.axis=1.5, cex.lab=1.5, ylab="Ne/1000", xlab = "Years ago")
lines(GL_C_2.3_109_ang_auto$year, GL_C_2.3_109_ang_auto$Ne_2.5./1000, col="green", lty=2, lwd=1)
lines(GL_C_2.3_109_ang_auto$year, GL_C_2.3_109_ang_auto$Ne_97.5./1000, col="green", lty=2, lwd=1)

lines(IS_C_2.3_109_ang_auto$year, IS_C_2.3_109_ang_auto$Ne_median/1000, col="red", lwd=1)
lines(IS_C_2.3_109_ang_auto$year, IS_C_2.3_109_ang_auto$Ne_2.5./1000, col="red", lty=2, lwd=1)
lines(IS_C_2.3_109_ang_auto$year, IS_C_2.3_109_ang_auto$Ne_97.5./1000, col="red", lty=2, lwd=1)

lines(NO_C_2.3_109_ang_auto$year, NO_C_2.3_109_ang_auto$Ne_median/1000, col="orange", lwd=1)
lines(NO_C_2.3_109_ang_auto$year, NO_C_2.3_109_ang_auto$Ne_2.5./1000, col="orange", lty=2, lwd=1)
lines(NO_C_2.3_109_ang_auto$year, NO_C_2.3_109_ang_auto$Ne_97.5./1000, col="orange", lty=2, lwd=1)

lines(DK_C_2.3_109_ang_auto$year, DK_C_2.3_109_ang_auto$Ne_median/1000, col="blue", lwd=1)
lines(DK_C_2.3_109_ang_auto$year, DK_C_2.3_109_ang_auto$Ne_2.5./1000, col="blue", lty=2, lwd=1)
lines(DK_C_2.3_109_ang_auto$year, DK_C_2.3_109_ang_auto$Ne_97.5./1000, col="blue", lty=2, lwd=1)
# legend("topleft",        # Add legend to plot
#        legend = c("GL_C", "IS_C", "NO_C", "DK_C"),
#        col = c("green", "red", "orange", "blue"),
#        pch = 16,
#        cex = 1.5)
#dev.off()

require(ggplot2)

ggplot(Con_SP)+
  geom_line(aes(year/1.95,(Ne_median/1000),colour=label),lwd=1.5)+
  geom_line(aes(year/1.95,(Ne_2.5./1000),colour=label),lty=3,lwd=1.5, alpha=0.3)+
  geom_line(aes(year/1.95,(Ne_97.5./1000),colour=label),lty=3,lwd=1.5, alpha=0.3)+
  xlab("Years ago")+
  ylab("Ne/1000")+
  coord_cartesian(xlim=c(0,200000),ylim=c(0,600))+
  #scale_linetype_manual(values = rep("solid",10),guide="none")+
  scale_color_manual(name="Sample",values = c("blue","green","red","orange"), labels=c("Denmark","Greenland","Iceland","Norway"))+
  #guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  scale_x_continuous(breaks = c(1000,seq(10000,200000,10000)))+
  scale_y_continuous(breaks = c(seq(0,600,10)))+
  geom_vline(xintercept=10000,lty=2,lwd=0.75)+
  geom_vline(xintercept=25000,lty=3,lwd=0.75)+
  geom_vline(xintercept=110000,lty=4,lwd=0.75)+
  geom_vline(xintercept=150000,lty=4,lwd=0.75)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))


ggplot(Con_SP)+
  geom_line(aes(year/1.95,(Ne_median),colour=label),lwd=1.5)+
  geom_line(aes(year/1.95,(Ne_2.5.),colour=label),lty=3,lwd=1.5, alpha=0.3)+
  geom_line(aes(year/1.95,(Ne_97.5.),colour=label),lty=3,lwd=1.5, alpha=0.3)+
  xlab("Years ago")+
  ylab("Ne")+
  #coord_cartesian(xlim=c(0,200),ylim=c(0,3000))+
  coord_cartesian(xlim=c(0,35000),ylim=c(0,150000))+
  #scale_linetype_manual(values = rep("solid",10),guide="none")+
  scale_color_manual(name="Sample",values = c("blue","green","red","orange"), labels=c("Denmark","Greenland","Iceland","Norway"))+
  #guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  scale_x_continuous(breaks = c(seq(0,35000,1000)))+
  scale_y_continuous(breaks = c(seq(0,150000,10000)))+
  geom_vline(xintercept=10000,lty=2,lwd=0.75)+
  geom_vline(xintercept=20000,lty=3,lwd=0.75)+
  geom_vline(xintercept=110000,lty=4,lwd=0.75)+
  geom_vline(xintercept=150000,lty=4,lwd=0.75)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))


ggplot(Con_SP)+
  geom_line(aes(year/2.223,(Ne_median/1000),colour=label),lwd=1.5)+
  geom_line(aes(year/2.223,(Ne_2.5./1000),colour=label),lty=3,lwd=1.5, alpha=0.3)+
  geom_line(aes(year/2.223,(Ne_97.5./1000),colour=label),lty=3,lwd=1.5, alpha=0.3)+
  xlab("Years ago")+
  ylab("Ne/1000")+
  coord_cartesian(xlim=c(0,200000),ylim=c(0,600))+
  #scale_linetype_manual(values = rep("solid",10),guide="none")+
  scale_color_manual(name="Sample",values = c("blue","green","red","orange"), labels=c("Denmark","Greenland","Iceland","Norway"))+
  #guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  scale_x_continuous(breaks = c(1000,seq(10000,200000,10000)))+
  geom_vline(xintercept=10000,lty=2,lwd=0.75)+
  geom_vline(xintercept=25000,lty=3,lwd=0.75)+
  geom_vline(xintercept=110000,lty=4,lwd=0.75)+
  geom_vline(xintercept=150000,lty=4,lwd=0.75)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))

#### het vs depth
mean_depth_true<-read.table("Truemeandepth_onlysitespresentinindividual.txt", header = T)
missingness<-read.table("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.imiss", header = T)
het_data_withdepth<-het_data
het_data_withdepth$mean_depth<-mean_depth_true$MEAN_DEPTH
het_data_withdepth$F_MISS<-missingness$F_MISS
