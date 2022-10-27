setwd("C:/Users/ccrha/Dropbox/PhD eagle Island/Analyser/ThirdPaper_fullgenome_modernandhistoric/")
#install.packages("ggplot2")
#install.packages("rlang")
#install.packages("broom")
library(ggplot2)
library(RcppCNPy)
library(MASS)
library(RColorBrewer)
library(tidyr)
library(tidyverse)
library(gridExtra) 
library(plotrix)
library(ggpubr)
library(dplyr)
library(gghighlight)
library(viridis)
library(htmltools)
#install.packages("PopGenome")
library(vcfR)
#library(PopGenome)
library(devtools)
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("gdsfmt")
#  BiocManager::install("SNPRelate")
#install_github("zhengxwen/gdsfmt")
#install_github("zhengxwen/SNPRelate")
#library(gdsfmt)
library(SNPRelate)
#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)
# install_github("https://github.com/shenglin-liu/vcf2sfs.git")
library(vcf2gt)
source("vcf2sfs.R")
#install.packages("qqman")
library(qqman)
library(plyr)
library(hierfstat)

####vcf2gt("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05.recode.vcf")

## Heterozygosity and inbreeding
het_data<-read.csv("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.het.csv", sep = "," , header = T)
colnames(het_data)[1]<-c("ORDER")
head(het_data)
het_data$O.HET<-het_data$N_SITES-het_data$O.HOM.
het_data$E.HET<-het_data$N_SITES-het_data$E.HOM.
het_data$het_frac<-het_data$O.HET/het_data$N_SITES
het_data$CO_TI<-factor(het_data$CO_TI, levels = c("GL_C","GL_H","IS_C","IS_H","NO_C","NO_H","DK_C","DK_H","EE_C","TU_H"))
het_data[order(het_data$O.HET),]
#pdf("Rplot_Heterozygosity_92ind_miss75_100321.pdf", height = 8, width = 13)
#png("Rplot_Heterozygosity_92ind_miss75_100321.png", height = 800, width = 1300)
par(mar=c(5.2,5.2,0.2,0.2))
boxplot(het_data$O.HET/het_data$N_SITES~het_data$CO_TI, ylab = "Observed heterozygosity", xlab = "", cex.lab=3, cex.axis=2)
mtext("Population_time", side = 1, line = 3.5, cex = 3 )

#pdf("Rplot_Heterozygosity_expected_92ind_miss75_100321.pdf", height = 8, width = 13)
#png("Rplot_Heterozygosity_expected_92ind_miss75_100321.png", height = 800, width = 1300)
par(mar=c(5.2,5.2,0.2,0.2))
boxplot(het_data$E.HET/het_data$N_SITES~het_data$CO_TI, ylab = "Expected heterozygosity", xlab = "", cex.lab=3, cex.axis=2)
mtext("Population_time", side = 1, line = 3.5, cex = 3 )
#dev.off()

het_obs_tmp<-het_data[,c(1,2,3,4,7,9)]
het_data_obs_exp_for_ggplot<-het_obs_tmp[order(het_obs_tmp$ORDER),]
het_data_obs_exp_for_ggplot$obs_or_exp<-rep("obs", 92)
colnames(het_data_obs_exp_for_ggplot)[6]<-c("HET")

het_exp_tmp<-het_data[,c(1,2,3,4,7,10)]
het_exp_tmp2<-het_exp_tmp[order(het_exp_tmp$ORDER),]
het_exp_tmp2$obs_or_exp<-rep("exp", 92)
colnames(het_exp_tmp2)[6]<-c("HET")

het_inb_tmp<-het_data[,c(1,2,3,4,7,8)]
het_inb_tmp2<-het_inb_tmp[order(het_inb_tmp$ORDER),]
het_inb_tmp2$obs_or_exp<-rep("inb", 92)
colnames(het_inb_tmp2)[6]<-c("HET")

het_data_obs_exp_for_ggplot_true<-rbind(het_data_obs_exp_for_ggplot, het_exp_tmp2)

ggplot_het_bothexpandobsandinb_090421<-
  ggplot(het_data_obs_exp_for_ggplot_true, aes(x=CO_TI, y=HET/N_SITES, color = obs_or_exp)) +
  geom_boxplot() +
  geom_boxplot(data = het_inb_tmp2, aes(x=CO_TI, y=HET, color = obs_or_exp), fill="#636363") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        title = element_text(size = 16), legend.text = element_text(size = 14)) +
  #annotate(geom= "text", x=seq_len(unique(het_data_obs_exp_for_ggplot_true$CO_TI)), y=10, label=het_data_obs_exp_for_ggplot_true$obs_or_exp) +
  scale_color_manual(values = c( "#aaaaaa", "#000000", "#000000")) + 
  labs(x = "Population_time", 
       y = "Heterozygosity & Inbreeding coefficient",
       color="") + 
  theme(legend.position = "none")
  
#ggsave(plot = ggplot_het_bothexpandobsandinb_090421, filename = "het_plot_bothhetandobsandinb_130421.png", width = 30, height = 20, units = "cm", device = "png")

  # het_data_obs_exp_for_ggplot_true %>% 
  #   mutate(CO_TI = fct_relevel(CO_TI, 
  #                             "GL_C", "GL_H", "IS_C", "IS_H", "NO_C", "NO_H", "DK_C", "DK_H", "EE_C", "TU_H")) %>%
  # ggplot( aes(x=interaction(CO_TI, obs_or_exp, lex.order = T), y=HET/N_SITES, color = obs_or_exp)) +
  #   geom_boxplot() +
  #   geom_boxplot(data = het_inb_tmp2, aes(x=interaction(CO_TI, obs_or_exp, lex.order = T), y=HET, color = obs_or_exp), fill="#636363") +
  #   theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
  #         title = element_text(size = 16), legend.text = element_text(size = 14)) +
  #   #annotate(geom= "text", x=seq_len(unique(het_data_obs_exp_for_ggplot_true$CO_TI)), y=10, label=het_data_obs_exp_for_ggplot_true$obs_or_exp) +
  #   scale_color_manual(values = c( "#aaaaaa", "#000000", "#000000")) + 
  #   labs(x = "Population_time", 
  #        y = "Heterozygosity & Inbreeding coefficient",
  #        color="") + 
  #   theme(legend.position = "none")

het_gg_observed_090421<-
  ggplot(het_data, aes(x=CO_TI, y=O.HET/N_SITES)) + 
       geom_boxplot() + 
  labs(title = "", 
      x = "Population_time", 
      y = "Observed heterozygosity") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        title = element_text(size = 16), legend.text = element_text(size = 14)) +
  scale_y_continuous(position = "left", limits = c(0.125, 0.3)) 

het_gg_expected_090421<-
  ggplot(het_data, aes(x=CO_TI, y=E.HET/N_SITES), ) + 
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

#png("Rplot_Inbreeding_92ind_miss75_100321.png", height = 800, width = 1300)
par(mar=c(5.2,5.2,0.2,0.2))
boxplot(het_data$F~het_data$CO_TI, ylab = "Inbreeding coefficient", xlab = "", cex.lab=3, cex.axis=2)
mtext("Population_time", side = 1, line = 3.5, cex = 3 )
#dev.off()

boxplot(het_data$O.HET~het_data$CO_TI)
boxplot(het_data$E.HET~het_data$CO_TI)

t.test(het_data$het_frac[het_data$CO_TI=="DK_C"], het_data$het_frac[het_data$CO_TI=="DK_H"], alternative = "two.sided")
mean(het_data$F[het_data$CO_TI=="DK_C"])

wilcox.test(het_data$het_frac[het_data$CO_TI=="GL_C"], het_data$het_frac[het_data$CO_TI=="GL_H"], exact = T)
wilcox.test(het_data$het_frac[het_data$CO_TI=="IS_C"], het_data$het_frac[het_data$CO_TI=="IS_H"], exact = T)
wilcox.test(het_data$het_frac[het_data$CO_TI=="NO_C"], het_data$het_frac[het_data$CO_TI=="NO_H"], exact = T)
wilcox.test(het_data$het_frac[het_data$CO_TI=="DK_C"], het_data$het_frac[het_data$CO_TI=="DK_H"], exact = T)

wilcox.test(het_data$het_frac[het_data$CO_TI=="NO_C" |het_data$CO_TI=="NO_H" |
                                het_data$CO_TI=="DK_C" |het_data$CO_TI=="NO_H" | 
                                het_data$CO_TI=="EE_C"| het_data$CO_TI=="TU"], 
            het_data$het_frac[het_data$CO_TI=="GL_C" | het_data$CO_TI=="GL_H" |
                                het_data$CO_TI=="IS_C" | het_data$CO_TI=="IS_H"
                              ], exact = T)
#t-test of Ts/Tv ratio
#modern 2.83
#historic 2.84 
#all 2.87
t.test(2.83, 2.84)

### MapDamage
mapd_5<-read.table("all_merged_mapDamage_5pCtoT_freq_indcol_92_altered.txt", header = F)
mapd_3<-read.table("all_merged_mapDamage_3pGtoA_freq_indcol_92_altered.txt", header = F)
colnames(mapd_3)<-c("pos","val", "ind", "pop_time", "time")
colnames(mapd_5)<-c("pos","val", "ind", "pop_time", "time")
mapdamage3_contemp_hist<-
  ggplot(data = mapd_3, aes(x=pos, y=val)) + geom_line(aes(colour=time)) +
  ylim(0.01,0.045) + 
  labs(title = "", 
     x ="Nucleotide position from the 3'end", 
     y = "Frequency of G to A substitutions",
     color="") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        title = element_text(size = 16), legend.text = element_text(size = 14)) +
  scale_color_manual(values = c("#aaaaaa", "#000000"))


mapdamage5_contemp_hist<-ggplot(data = mapd_5, aes(x=pos*-1, y=val)) + geom_line(aes(colour=time)) +
  scale_y_continuous(position = "right", limits = c(0.01, 0.045)) + 
  scale_x_continuous(labels=c("25","20","15","10", "5","0")) +
  labs(title = "", 
       x ="Nucleotide position from the 5'end", 
       y = "Frequency of C to T substitutions",
       color="") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        title = element_text(size = 16), legend.text = element_text(size = 14)) +
  scale_color_manual(values = c("#aaaaaa", "#000000"))


mapdamage_2figures_290321<-
  ggarrange(mapdamage3_contemp_hist, mapdamage5_contemp_hist,
            labels = c("A", "B"),
            ncol = 2, nrow = 1, common.legend = TRUE)

#ggsave(plot = mapdamage_2figures_290321, filename = "MapDamage_3and5_300321.png", width = 30, height = 20, units = "cm", device = "png")

ggplot(data = mapd_5, aes(x=pos, y=val)) + geom_line(aes(colour=pop_time))



map_109<-read.table("all_results_mac1_109ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05.map")

#### PCA 

pcadata<-read.table("smartpca/out_smartPCA_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26_project_dirfromplink_evec", header=F)
pcadata<-pcadata[,c(1:11)]
colnames(pcadata)<-c("IND","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
fam<-read.table("smartpca/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26_goodinds.fam")
pcadata_name<-as.data.frame(fam[,c(1)])
pcadata_name[,c(2:12)]<-pcadata[,c(1:11)]
colnames(pcadata_name)[1]<-c("POP")

pop<-as.factor(pcadata_name$POP)

myCol<-c(pop)
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}
myColAlpha<-add.alpha(myCol, alpha=0.5)
library(maptools)
#pdf("PCA_133ind_PCangsd_sexsplit_190220.pdf", width = 13, height = 8)
plot(pcadata_name$PC1, pcadata_name$PC2, col=pop, pch=20, cex=2.5, cex.lab=1.5) #, xlab=paste("PC1, variance", percent[9], "%", sep=" "), ylab=paste("PC2, variance", percent[10], "%", sep=" ")
legend("bottomright", legend=levels(pop), text.col=1:nlevels(pop), cex = .5)
#pointLabel(PC$PC1, PC$PC2, labels=sample, cex=0.5)
#dev.off()

### Projected individuals 

pcadata_p<-read.table("smartpca/out_smartPCA_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26_project_dirfromplink_evec", header=F)
pcadata_p<-pcadata_p[,c(1:11)]
colnames(pcadata_p)<-c("IND","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
fam<-read.table("smartpca/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26_goodinds.fam")
pcadata_p_name<-as.data.frame(fam[,c(1)])
pcadata_p_name[,c(2:12)]<-pcadata_p[,c(1:11)]
colnames(pcadata_p_name)[1]<-c("POP")

pop<-as.factor(pcadata_p_name$POP)

myCol<-c(pop)
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}
myColAlpha<-add.alpha(myCol, alpha=0.5)
library(maptools)
#pdf("PCA_133ind_PCangsd_sexsplit_190220.pdf", width = 13, height = 8)
plot(pcadata_p_name$PC1, pcadata_p_name$PC2, col=pop, pch=20, cex=2.5, cex.lab=1.5) #, xlab=paste("PC1, variance", percent[9], "%", sep=" "), ylab=paste("PC2, variance", percent[10], "%", sep=" ")
legend("bottomright", legend=levels(pop), text.col=1:nlevels(pop), cex = .5)

##plotting more PCs.
pairs(~ PC2 + PC2 + PC3 + PC4 + PC5, pch=19, data=pcadata_p, col=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"))
pairs(pcadata_name[3:7], col=myColAlpha)

eigenvals_project<-read.table("smartpca/out_smartPCA_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26_project_dirfromplink_evec", header = F)
eigenvals_project

#pca_92_130921 <- 
pcadata_p_name %>% 
  ggplot(aes(y=PC2, x=PC1, color=POP, size=POP)) + 
  scale_color_manual(values = c("#1F78B4", "#A6CEE3", "#6A3D9A", "#33A02C", "#B2DF8A", "#E31A1C", "#FB9A99", "#FF7F00", "#FDBF6F", "#CAB2D6"),  name="Pop_Time") + 
  scale_size_manual(values=c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5), element_blank(), name="") +
  guides(size = F) + 
  #ylim(0,500) +
  geom_point() +
  labs(x ="PC1 45.34%", y = "PC2 15.05%") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
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
    text = element_text(size = 14))

#ggsave(plot = pca_92_130921, filename = "smartPCA_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26_project_dirfromplink_100321_130921.png", width = 30, height = 20, units = "cm", device = "png")

pairs(~PC1 + PC2 + PC2 + PC3 + PC4 + PC5, pch=19, data=pcadata_p_name, col=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"))


library(RcppCNPy)
data <- as.matrix(read.table("pcangsd/PCangsd/test_from_plink_280121.cov")) # Reads in estimated covariance matrix
eig<-eigen(data, symm=T)
eig$val <- eig$val/sum(eig$val)
percent<-round(eig$val*100, digits=2)
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))
plot(PC$PC1, -1*PC$PC2, xlab=paste("PC1, variance", percent[1], "%", sep=" "), ylab=paste("PC2, variance", percent[2], "%", sep=" "))

e <- eigen(data)
plot(e$vectors[,1], e$vectors[,2]*-1,xlab="PC1",ylab="PC2", main="individual allele frequency")

#### ADMIXTURE 
fam<-read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26_goodinds.fam")
k1<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.1.Q")))
k2<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.2.Q")))
k3<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.3.Q")))
k4<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.4.Q")))
k5<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.5.Q")))
k6<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.6.Q")))
k7<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.7.Q")))
k8<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.8.Q")))
k9<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.9.Q")))
k10<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.10.Q")))
k11<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.11.Q")))
k12<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.12.Q")))
k13<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.13.Q")))
k14<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.14.Q")))
k15<-t(as.matrix(read.table("admixture/all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr26.15.Q")))

pop_order<-as.data.frame(c(het_data[,c(1:4)]))
pop_c_order<-pop_order[order(pop_order[,1]),]


# brewer.pal(12,"Paired")
# "#A6CEE3" - light blue
# "#1F78B4" - blue
# "#B2DF8A" - light green
# "#33A02C" - green
# "#FB9A99" - light red 
# "#E31A1C" - red
# "#FDBF6F" - light orange
# "#FF7F00" - orange 
# "#CAB2D6" - light purple
# "#6A3D9A" - purple
# "#FFFF99" - light yellow 
# "#B15928" - brown

pdf("admixture_92ind_6together_140921.pdf", width=13, height=10)
admix.par<-par(mfrow=c(3,2))

#tissemand2<-palette(c("#A6CEE3", "#1F78B4"))
#tissemand2<-palette(c("#A6CEE3", "#1F78B4"))

tissemand2<-palette(c("#FF7F00", "#33A02C"))
tissemand2<-palette(c("#FF7F00", "#33A02C"))

#k2
admixk2<-k2[,order(pop_order[,1])]
par(mar=c(1,2,1.5,0.5))
#png("Rplot_admixture_k2_109ind_miss75_projection_270121.png", height = 800, width = 1300)
barplot(admixk2, col=tissemand2, space=0, border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
title("A", adj=0, line = 0.2, cex.main=2)
#text(x=tapply(x1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
abline(v=12)
abline(v=20)
abline(v=45)
abline(v=47)
abline(v=59)
abline(v=72)
abline(v=83)
abline(v=88)
abline(v=91)
#dev.off()

# tissemand3<-palette(c("#1F78B4", "#A6CEE3", "#B2DF8A"))
# tissemand3<-palette(c("#1F78B4", "#A6CEE3", "#B2DF8A"))

tissemand3<-palette(c("#33A02C", "#FF7F00", "#1F78B4"))
tissemand3<-palette(c("#33A02C", "#FF7F00", "#1F78B4")) 

#k3
admixk3<-k3[,order(pop_order[,1])]
par(mar=c(1,0.5,1.5,1))
#png("Rplot_admixture_k3_109ind_miss75_projection_270121.png", height = 800, width = 1300)
barplot(admixk3,col=tissemand3, space=0, border=NA, yaxt="n") # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
title("B", adj=0, line = 0.2, cex.main=2)
#text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
abline(v=12)
abline(v=20)
abline(v=45)
abline(v=47)
abline(v=59)
abline(v=72)
abline(v=83)
abline(v=88)
abline(v=91)
#dev.off()

# "#A6CEE3" - light blue
# "#1F78B4" - blue
# "#B2DF8A" - light green
# "#33A02C" - green

# tissemand4<-palette(c("#A6CEE3", "#B2DF8A", "#1F78B4", "#33A02C"))
# tissemand4<-palette(c("#A6CEE3", "#B2DF8A", "#1F78B4", "#33A02C"))

tissemand4<-palette(c("#FF7F00", "#1F78B4", "#E31A1C", "#33A02C"))
tissemand4<-palette(c("#FF7F00", "#1F78B4", "#E31A1C", "#33A02C"))

#k4
admixk4<-k4[,order(pop_order[,1])]
par(mar=c(1,2,1.5,0.5))
#png("Rplot_admixture_k4_109ind_miss75_projection_270121.png", height = 800, width = 1300)
barplot(admixk4,col=tissemand4,space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
title("C", adj=0, line = 0.2, cex.main=2)
#text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
abline(v=12)
abline(v=20)
abline(v=45)
abline(v=47)
abline(v=59)
abline(v=72)
abline(v=83)
abline(v=88)
abline(v=91)
#dev.off()

# tissemand5<-palette(c("#A6CEE3", "#B2DF8A", "#1F78B4", "#33A02C", "#FB9A99"))
# tissemand5<-palette(c("#A6CEE3", "#B2DF8A", "#1F78B4", "#33A02C", "#FB9A99"))

tissemand5<-palette(c("#FF7F00", "#1F78B4", "#FB9A99", "#33A02C", "#E31A1C"))
tissemand5<-palette(c("#FF7F00", "#1F78B4", "#FB9A99", "#33A02C", "#E31A1C"))

#k5
admixk5<-k5[,order(pop_order[,1])]
par(mar=c(1,0.5,1.5,1))
#png("Rplot_admixture_k5_109ind_miss75_projection_270121.png", height = 800, width = 1300)
barplot(admixk5,col=tissemand5,space=0,border=NA, yaxt="n") # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
title("D", adj=0, line = 0.2, cex.main=2)
#text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
abline(v=12)
abline(v=20)
abline(v=45)
abline(v=47)
abline(v=59)
abline(v=72)
abline(v=83)
abline(v=88)
abline(v=91)
#dev.off()

# "#A6CEE3" - light blue
# "#1F78B4" - blue
# "#B2DF8A" - light green
# "#33A02C" - green
# "#FB9A99" - light red 
# "#E31A1C" - red

# tissemand6<-palette(c("#33A02C", "#A6CEE3", "#FB9A99", "#E31A1C", "#B2DF8A", "#1F78B4"))
# tissemand6<-palette(c("#33A02C", "#A6CEE3", "#FB9A99", "#E31A1C", "#B2DF8A", "#1F78B4"))

tissemand6<-palette(c("#33A02C", "#FF7F00", "#E31A1C",  "#6A3D9A","#1F78B4", "#FB9A99"))
tissemand6<-palette(c("#33A02C", "#FF7F00", "#E31A1C",  "#6A3D9A","#1F78B4", "#FB9A99"))

#k6
admixk6<-k6[,order(pop_order[,1])]
par(mar=c(4.5,2,1.5,0.5))
#png("Rplot_admixture_k6_109ind_miss75_projection_270121.png", height = 800, width = 1300)
barplot(admixk6,col=tissemand6,space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
title("E", adj=0, line = 0.2, cex.main=2)
text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=2, offset=0, pos=2)
abline(v=12)
abline(v=20)
abline(v=45)
abline(v=47)
abline(v=59)
abline(v=72)
abline(v=83)
abline(v=88)
abline(v=91)
#dev.off()

# "#A6CEE3" - light blue
# "#1F78B4" - blue
# "#B2DF8A" - light green
# "#33A02C" - green
# "#FB9A99" - light red 
# "#E31A1C" - red
# "#FDBF6F" - light orange

# tissemand7<-palette(c("#33A02C", "#B2DF8A",  "#E31A1C", "#FDBF6F", "#FB9A99",  "#A6CEE3","#1F78B4" ))
# tissemand7<-palette(c("#33A02C", "#B2DF8A",  "#E31A1C", "#FDBF6F", "#FB9A99", "#A6CEE3", "#1F78B4" ))

tissemand7<-palette(c("#33A02C", "#1F78B4", "#6A3D9A", "#FDBF6F", "#E31A1C", "#FF7F00", "#FB9A99"))
tissemand7<-palette(c("#33A02C", "#1F78B4", "#6A3D9A", "#FDBF6F", "#E31A1C", "#FF7F00", "#FB9A99"))

#k7
admixk7<-k7[,order(pop_order[,1])]
par(mar=c(4.5,0.5,1.5,1))
#png("Rplot_admixture_k7_109ind_miss75_projection_270121.png", height = 800, width = 1300)
barplot(admixk7,col=tissemand7,space=0,border=NA, yaxt="n") # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
title("F", adj=0, line = 0.2, cex.main=2)
text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=2, offset=0, pos=2)
abline(v=12)
abline(v=20)
abline(v=45)
abline(v=47)
abline(v=59)
abline(v=72)
abline(v=83)
abline(v=88)
abline(v=91)
#dev.off()

dev.off()

#k8
admixk8<-k8[,order(pop_order[,1])]
par(mar=c(3.6,2,1,1))
#png("Rplot_admixture_k8_92ind_miss75_projection_110315.png", height = 800, width = 1300)
barplot(admixk8,col=brewer.pal(8,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
abline(v=12)
abline(v=20)
abline(v=45)
abline(v=47)
abline(v=59)
abline(v=72)
abline(v=83)
abline(v=88)
abline(v=91)
#dev.off()


#k9
admixk9<-k9[,order(pop_order[,1])]
par(mar=c(3.6,2,1,1))
#png("Rplot_admixture_k9_92ind_miss75_projection_110315.png", height = 800, width = 1300)
barplot(admixk9,col=brewer.pal(9,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
abline(v=12)
abline(v=20)
abline(v=45)
abline(v=47)
abline(v=59)
abline(v=72)
abline(v=83)
abline(v=88)
abline(v=91)
#dev.off()

#k10
admixk10<-k10[,order(pop_order[,1])]
par(mar=c(3.6,2,1,1))
#png("Rplot_admixture_k10_92ind_miss75_projection_110315.png", height = 800, width = 1300)
barplot(admixk10,col=brewer.pal(10,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
abline(v=12)
abline(v=20)
abline(v=45)
abline(v=47)
abline(v=59)
abline(v=72)
abline(v=83)
abline(v=88)
abline(v=91)
#dev.off()

#k11
admixk11<-k11[,order(pop_order[,1])]
par(mar=c(3.6,2,1,1))
#png("Rplot_admixture_k11_92ind_miss75_projection_110315.png", height = 800, width = 1300)
barplot(admixk11,col=brewer.pal(11,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
abline(v=12)
abline(v=20)
abline(v=45)
abline(v=47)
abline(v=59)
abline(v=72)
abline(v=83)
abline(v=88)
abline(v=91)
#dev.off()

#k12
admixk12<-k12[,order(pop_order[,1])]
par(mar=c(3.6,2,1,1))
#png("Rplot_admixture_k12_92ind_miss75_projection_110315.png", height = 800, width = 1300)
barplot(admixk12,col=brewer.pal(12,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
abline(v=12)
abline(v=20)
abline(v=45)
abline(v=47)
abline(v=59)
abline(v=72)
abline(v=83)
abline(v=88)
abline(v=91)
#dev.off()

#k13
admixk13<-k13[,order(pop_order[,1])]
par(mar=c(3.6,2,1,1))
#png("Rplot_admixture_k13_92ind_miss75_projection_110315.png", height = 800, width = 1300)
barplot(admixk13,col=viridis(13),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
abline(v=12)
abline(v=20)
abline(v=45)
abline(v=47)
abline(v=59)
abline(v=72)
abline(v=83)
abline(v=88)
abline(v=91)
#dev.off()

#k14
admixk14<-k14[,order(pop_order[,1])]
par(mar=c(3.6,2,1,1))
#png("Rplot_admixture_k14_92ind_miss75_projection_110315.png", height = 800, width = 1300)
barplot(admixk14,col=viridis(14),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
abline(v=12)
abline(v=20)
abline(v=45)
abline(v=47)
abline(v=59)
abline(v=72)
abline(v=83)
abline(v=88)
abline(v=91)
#dev.off()

#k15
admixk15<-k15[,order(pop_order[,1])]
par(mar=c(3.6,2,1,1))
#png("Rplot_admixture_k15_92ind_miss75_projection_110315.png", height = 800, width = 1300)
barplot(admixk15,col=viridis(15),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
abline(v=12)
abline(v=20)
abline(v=45)
abline(v=47)
abline(v=59)
abline(v=72)
abline(v=83)
abline(v=88)
abline(v=91)
#dev.off()

# barplot(k8,col=brewer.pal(8,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop),pop[,8],mean),-0.05,unique(pop[,6]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=42)
# abline(v=54)
# abline(v=59)
# abline(v=70)
# abline(v=73)
# abline(v=85)
# abline(v=93)
# abline(v=95)
# abline(v=108)

# ############################
# ############################
# ############################
# ### test only 15 Icelandic
# 
# fam<-read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.fam")
# k1<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.1.Q")))
# k2<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.2.Q")))
# k3<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.3.Q")))
# k4<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.4.Q")))
# k5<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.5.Q")))
# k6<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.6.Q")))
# k7<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.7.Q")))
# k8<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.8.Q")))
# k9<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.9.Q")))
# k10<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.10.Q")))
# k11<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.11.Q")))
# k12<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.12.Q")))
# k13<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.13.Q")))
# k14<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.14.Q")))
# k15<-t(as.matrix(read.table("admixture/all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newchr.15.Q")))
# 
# het_data_82<-read.csv("all_results_mac1_82ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.het.csv", sep = "," , header = T)
# 
# pop_order<-as.data.frame(c(het_data_82[,c(1:4)]))
# pop_c_order<-pop_order[order(pop_order[,1]),]
# 
# #k2
# admixk2<-k2[,order(pop_order[,1])]
# #par(mar=c(10,5,5,1))
# #png("Rplot_admixture_k2_109ind_miss75_projection_270121.png", height = 800, width = 1300)
# barplot(admixk2,col=brewer.pal(8,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=12)
# abline(v=20)
# abline(v=45)
# abline(v=47)
# abline(v=59)
# abline(v=72)
# abline(v=83)
# abline(v=88)
# abline(v=91)
# #dev.off()
# 
# #k3
# admixk3<-k3[,order(pop_order[,1])]
# par(mar=c(3.6,2,1,1))
# #png("Rplot_admixture_k3_109ind_miss75_projection_270121.png", height = 800, width = 1300)
# barplot(admixk3,col=brewer.pal(8,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=12)
# abline(v=20)
# abline(v=45)
# abline(v=47)
# abline(v=59)
# abline(v=72)
# abline(v=83)
# abline(v=88)
# abline(v=91)
# #dev.off()
# 
# #k4
# admixk4<-k4[,order(pop_order[,1])]
# par(mar=c(3.6,2,1,1))
# #png("Rplot_admixture_k4_109ind_miss75_projection_270121.png", height = 800, width = 1300)
# barplot(admixk4,col=brewer.pal(8,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=12)
# abline(v=20)
# abline(v=45)
# abline(v=47)
# abline(v=59)
# abline(v=72)
# abline(v=83)
# abline(v=88)
# abline(v=91)
# #dev.off()
# 
# #k5
# admixk5<-k5[,order(pop_order[,1])]
# par(mar=c(3.6,2,1,1))
# #png("Rplot_admixture_k5_109ind_miss75_projection_270121.png", height = 800, width = 1300)
# barplot(admixk5,col=brewer.pal(8,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=12)
# abline(v=20)
# abline(v=45)
# abline(v=47)
# abline(v=59)
# abline(v=72)
# abline(v=83)
# abline(v=88)
# abline(v=91)
# #dev.off()
# 
# #k6
# admixk6<-k6[,order(pop_order[,1])]
# par(mar=c(3.6,2,1,1))
# #png("Rplot_admixture_k6_109ind_miss75_projection_270121.png", height = 800, width = 1300)
# barplot(admixk6,col=brewer.pal(8,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=12)
# abline(v=20)
# abline(v=45)
# abline(v=47)
# abline(v=59)
# abline(v=72)
# abline(v=83)
# abline(v=88)
# abline(v=91)
# #dev.off()
# 
# #k7
# admixk7<-k7[,order(pop_order[,1])]
# par(mar=c(3.6,2,1,1))
# #png("Rplot_admixture_k7_109ind_miss75_projection_270121.png", height = 800, width = 1300)
# barplot(admixk7,col=brewer.pal(8,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=12)
# abline(v=20)
# abline(v=45)
# abline(v=47)
# abline(v=59)
# abline(v=72)
# abline(v=83)
# abline(v=88)
# abline(v=91)
# #dev.off()
# 
# #k8
# admixk8<-k8[,order(pop_order[,1])]
# par(mar=c(3.6,2,1,1))
# #png("Rplot_admixture_k8_109ind_miss75_projection_270121.png", height = 800, width = 1300)
# barplot(admixk8,col=brewer.pal(8,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=12)
# abline(v=20)
# abline(v=45)
# abline(v=47)
# abline(v=59)
# abline(v=72)
# abline(v=83)
# abline(v=88)
# abline(v=91)
# #dev.off()
# 
# #k9
# admixk9<-k9[,order(pop_order[,1])]
# par(mar=c(3.6,2,1,1))
# #png("Rplot_admixture_k9_109ind_miss75_projection_270121.png", height = 800, width = 1300)
# barplot(admixk9,col=brewer.pal(9,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=12)
# abline(v=20)
# abline(v=45)
# abline(v=47)
# abline(v=59)
# abline(v=72)
# abline(v=83)
# abline(v=88)
# abline(v=91)
# #dev.off()
# 
# #k10
# admixk10<-k10[,order(pop_order[,1])]
# par(mar=c(3.6,2,1,1))
# #png("Rplot_admixture_k10_109ind_miss75_projection_270121.png", height = 800, width = 1300)
# barplot(admixk10,col=brewer.pal(10,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=12)
# abline(v=20)
# abline(v=45)
# abline(v=47)
# abline(v=59)
# abline(v=72)
# abline(v=83)
# abline(v=88)
# abline(v=91)
# #dev.off()
# 
# #k11
# admixk11<-k11[,order(pop_order[,1])]
# par(mar=c(3.6,2,1,1))
# #png("Rplot_admixture_k11_109ind_miss75_projection_270121.png", height = 800, width = 1300)
# barplot(admixk11,col=brewer.pal(11,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=12)
# abline(v=20)
# abline(v=45)
# abline(v=47)
# abline(v=59)
# abline(v=72)
# abline(v=83)
# abline(v=88)
# abline(v=91)
# #dev.off()
# 
# #k12
# admixk12<-k12[,order(pop_order[,1])]
# par(mar=c(3.6,2,1,1))
# #png("Rplot_admixture_k12_109ind_miss75_projection_270121.png", height = 800, width = 1300)
# barplot(admixk12,col=brewer.pal(12,"Paired"),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=12)
# abline(v=20)
# abline(v=45)
# abline(v=47)
# abline(v=59)
# abline(v=72)
# abline(v=83)
# abline(v=88)
# abline(v=91)
# #dev.off()
# 
# #k13
# admixk13<-k13[,order(pop_order[,1])]
# par(mar=c(3.6,2,1,1))
# #png("Rplot_admixture_k13_109ind_miss75_projection_270121.png", height = 800, width = 1300)
# barplot(admixk13,col=viridis(13),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=12)
# abline(v=20)
# abline(v=45)
# abline(v=47)
# abline(v=59)
# abline(v=72)
# abline(v=83)
# abline(v=88)
# abline(v=91)
# #dev.off()
# 
# #k14
# admixk14<-k14[,order(pop_order[,1])]
# par(mar=c(3.6,2,1,1))
# #png("Rplot_admixture_k14_109ind_miss75_projection_270121.png", height = 800, width = 1300)
# barplot(admixk14,col=viridis(14),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=12)
# abline(v=20)
# abline(v=45)
# abline(v=47)
# abline(v=59)
# abline(v=72)
# abline(v=83)
# abline(v=88)
# abline(v=91)
# #dev.off()
# 
# #k10
# admixk15<-k15[,order(pop_order[,1])]
# par(mar=c(3.6,2,1,1))
# #png("Rplot_admixture_k10_109ind_miss75_projection_270121.png", height = 800, width = 1300)
# barplot(admixk15,col=viridis(15),space=0,border=NA) # ,space=0,border=NA,ylab="Admixture", cex.lab=1.8 ,main="NGSadmix K3",cex.main=2)
# text(tapply(1:nrow(pop_c_order),pop_c_order[,3],mean),-0.05,unique(pop_c_order[,3]),xpd=T,srt=70,cex=1.3, offset=0, pos=2)
# abline(v=12)
# abline(v=20)
# abline(v=45)
# abline(v=47)
# abline(v=59)
# abline(v=72)
# abline(v=83)
# abline(v=88)
# abline(v=91)
# #dev.off()


#####################################
################ nucleotide diversity

nuc_div_all_50k<-read.table("Nucleotide diversity/nucleotide_diversity_win50k_all.windowed.pi", header = T)
nuc_div_GL_modern_50k<-read.table("Nucleotide diversity/nucleotide_diversity_win50k_modern_GL.windowed.pi", header = T)
nuc_div_GL_hist_50k<-read.table("Nucleotide diversity/nucleotide_diversity_win50k_historic_GL.windowed.pi", header = T)
nuc_div_IS_modern_50k<-read.table("Nucleotide diversity/nucleotide_diversity_win50k_modern_IS.windowed.pi", header = T)
nuc_div_IS_hist_50k<-read.table("Nucleotide diversity/nucleotide_diversity_win50k_historic_IS.windowed.pi", header = T)
nuc_div_NO_modern_50k<-read.table("Nucleotide diversity/nucleotide_diversity_win50k_modern_NO.windowed.pi", header = T)
nuc_div_NO_hist_50k<-read.table("Nucleotide diversity/nucleotide_diversity_win50k_historic_NO.windowed.pi", header = T)
nuc_div_DK_modern_50k<-read.table("Nucleotide diversity/nucleotide_diversity_win50k_modern_DK.windowed.pi", header = T)
nuc_div_DK_hist_50k<-read.table("Nucleotide diversity/nucleotide_diversity_win50k_historic_DK.windowed.pi", header = T)
nuc_div_EE_modern_50k<-read.table("Nucleotide diversity/nucleotide_diversity_win50k_modern_EE.windowed.pi", header = T)
nuc_div_TU_hist_50k<-read.table("Nucleotide diversity/nucleotide_diversity_win50k_historic_TU.windowed.pi", header = T)

sum(nuc_div_all_50k$N_VARIANTS)
sum(nuc_div_GL_modern_50k$N_VARIANTS)
sum(nuc_div_GL_hist_50k$N_VARIANTS)
sum(nuc_div_IS_modern_50k$N_VARIANTS)
sum(nuc_div_IS_hist_50k$N_VARIANTS)
sum(nuc_div_NO_modern_50k$N_VARIANTS)
sum(nuc_div_NO_hist_50k$N_VARIANTS)
sum(nuc_div_DK_modern_50k$N_VARIANTS)
sum(nuc_div_DK_hist_50k$N_VARIANTS)
sum(nuc_div_EE_modern_50k$N_VARIANTS)
sum(nuc_div_TU_hist_50k$N_VARIANTS)

mean(nuc_div_all_50k$PI)*10e3
sd(nuc_div_all_50k$PI)*10e3
mean(nuc_div_GL_modern_50k$PI)*10e3
sd(nuc_div_GL_modern_50k$PI)*10e3
mean(nuc_div_GL_hist_50k$PI)*10e3
sd(nuc_div_GL_hist_50k$PI)*10e3
mean(nuc_div_IS_modern_50k$PI)*10e3
sd(nuc_div_IS_modern_50k$PI)*10e3
mean(nuc_div_IS_hist_50k$PI)*10e3
sd(nuc_div_IS_hist_50k$PI)*10e3
mean(nuc_div_NO_modern_50k$PI)*10e3
sd(nuc_div_NO_modern_50k$PI)*10e3
mean(nuc_div_NO_hist_50k$PI)*10e3
sd(nuc_div_NO_hist_50k$PI)*10e3
mean(nuc_div_DK_modern_50k$PI)*10e3
sd(nuc_div_DK_modern_50k$PI)*10e3
mean(nuc_div_DK_hist_50k$PI)*10e3
sd(nuc_div_DK_hist_50k$PI)*10e3
mean(nuc_div_EE_modern_50k$PI)*10e3
sd(nuc_div_EE_modern_50k$PI)*10e3
mean(nuc_div_TU_hist_50k$PI)*10e3
sd(nuc_div_TU_hist_50k$PI)*10e3

mean((nuc_div_all_50k$PI/50000)*nuc_div_all_50k$N_VARIANTS)
sd(nuc_div_all_50k$PI)*10e3
mean((nuc_div_GL_modern_50k$PI/50000)*nuc_div_GL_modern_50k$N_VARIANTS)
sd(nuc_div_GL_modern_50k$PI)*10e3
mean((nuc_div_GL_hist_50k$PI/50000)*nuc_div_GL_hist_50k$N_VARIANTS)
sd(nuc_div_GL_hist_50k$PI)*10e3
mean((nuc_div_IS_modern_50k$PI/50000)*nuc_div_IS_modern_50k$N_VARIANTS)
sd(nuc_div_IS_modern_50k$PI)*10e3
mean((nuc_div_IS_hist_50k$PI/50000)*nuc_div_IS_hist_50k$N_VARIANTS)
sd(nuc_div_IS_hist_50k$PI)*10e3
mean((nuc_div_NO_modern_50k$PI/50000)*nuc_div_NO_modern_50k$N_VARIANTS)
sd(nuc_div_NO_modern_50k$PI)*10e3
mean((nuc_div_NO_hist_50k$PI/50000)*nuc_div_NO_hist_50k$N_VARIANTS)
sd(nuc_div_NO_hist_50k$PI)*10e3
mean((nuc_div_DK_modern_50k$PI/50000)*nuc_div_DK_modern_50k$N_VARIANTS)
sd(nuc_div_DK_modern_50k$PI)*10e3
mean((nuc_div_DK_hist_50k$PI/50000)*nuc_div_DK_hist_50k$N_VARIANTS)
sd(nuc_div_DK_hist_50k$PI)*10e3
mean((nuc_div_EE_modern_50k$PI/50000)*nuc_div_EE_modern_50k$N_VARIANTS)
sd(nuc_div_EE_modern_50k$PI)*10e3
mean((nuc_div_TU_hist_50k$PI/50000)*nuc_div_TU_hist_50k$N_VARIANTS)
sd(nuc_div_TU_hist_50k$PI)*10e3


nuc_div_all_whole<-read.table("Nucleotide diversity/nucleotide_diversity_wholegenome_all.windowed.pi", header = T)
nuc_div_GL_modern_whole<-read.table("Nucleotide diversity/nucleotide_diversity_wholegenome_modern_GL.windowed.pi", header = T)
nuc_div_GL_hist_whole<-read.table("Nucleotide diversity/nucleotide_diversity_wholegenome_historic_GL.windowed.pi", header = T)
nuc_div_IS_modern_whole<-read.table("Nucleotide diversity/nucleotide_diversity_wholegenome_modern_IS.windowed.pi", header = T)
nuc_div_IS_hist_whole<-read.table("Nucleotide diversity/nucleotide_diversity_wholegenome_historic_IS.windowed.pi", header = T)
nuc_div_NO_modern_whole<-read.table("Nucleotide diversity/nucleotide_diversity_wholegenome_modern_NO.windowed.pi", header = T)
nuc_div_NO_hist_whole<-read.table("Nucleotide diversity/nucleotide_diversity_wholegenome_historic_NO.windowed.pi", header = T)
nuc_div_DK_modern_whole<-read.table("Nucleotide diversity/nucleotide_diversity_wholegenome_modern_DK.windowed.pi", header = T)
nuc_div_DK_hist_whole<-read.table("Nucleotide diversity/nucleotide_diversity_wholegenome_historic_DK.windowed.pi", header = T)
nuc_div_EE_modern_whole<-read.table("Nucleotide diversity/nucleotide_diversity_wholegenome_modern_EE.windowed.pi", header = T)
nuc_div_TU_hist_whole<-read.table("Nucleotide diversity/nucleotide_diversity_wholegenome_historic_TU.windowed.pi", header = T)

mean(nuc_div_all_whole$PI)
sd(nuc_div_all_whole$PI)
mean(nuc_div_GL_modern_whole$PI)
sd(nuc_div_GL_modern_whole$PI)
mean(nuc_div_GL_hist_whole$PI)
sd(nuc_div_GL_hist_whole$PI)
mean(nuc_div_IS_modern_whole$PI)
sd(nuc_div_IS_modern_whole$PI)
mean(nuc_div_IS_hist_whole$PI)
sd(nuc_div_IS_hist_whole$PI)
mean(nuc_div_NO_modern_whole$PI)
sd(nuc_div_NO_modern_whole$PI)
mean(nuc_div_NO_hist_whole$PI)
sd(nuc_div_NO_hist_whole$PI)
mean(nuc_div_DK_modern_whole$PI)
sd(nuc_div_DK_modern_whole$PI)
mean(nuc_div_DK_hist_whole$PI)
sd(nuc_div_DK_hist_whole$PI)
mean(nuc_div_EE_modern_whole$PI)
sd(nuc_div_EE_modern_whole$PI)
mean(nuc_div_TU_hist_whole$PI)
sd(nuc_div_TU_hist_whole$PI)

nuc_div_all_sites<-read.table("Nucleotide diversity/nucleotide_diversity_sites_all.sites.pi", header = T)
nuc_div_GL_modern_sites<-read.table("Nucleotide diversity/nucleotide_diversity_sites_modern_GL.sites.pi", header = T)
nuc_div_GL_hist_sites<-read.table("Nucleotide diversity/nucleotide_diversity_sites_historic_GL.sites.pi", header = T)
nuc_div_IS_modern_sites<-read.table("Nucleotide diversity/nucleotide_diversity_sites_modern_IS.sites.pi", header = T)
nuc_div_IS_hist_sites<-read.table("Nucleotide diversity/nucleotide_diversity_sites_historic_IS.sites.pi", header = T)
nuc_div_NO_modern_sites<-read.table("Nucleotide diversity/nucleotide_diversity_sites_modern_NO.sites.pi", header = T)
nuc_div_NO_hist_sites<-read.table("Nucleotide diversity/nucleotide_diversity_sites_historic_NO.sites.pi", header = T)
nuc_div_DK_modern_sites<-read.table("Nucleotide diversity/nucleotide_diversity_sites_modern_DK.sites.pi", header = T)
nuc_div_DK_hist_sites<-read.table("Nucleotide diversity/nucleotide_diversity_sites_historic_DK.sites.pi", header = T)
nuc_div_EE_modern_sites<-read.table("Nucleotide diversity/nucleotide_diversity_sites_modern_EE.sites.pi", header = T)
nuc_div_TU_hist_sites<-read.table("Nucleotide diversity/nucleotide_diversity_sites_historic_TU.sites.pi", header = T)

mean(nuc_div_all_sites$PI, na.rm = T)
sd(nuc_div_all_sites$PI, na.rm = T)
length(na.omit(nuc_div_all_sites$PI))
mean(nuc_div_GL_modern_sites$PI, na.rm = T)
sd(nuc_div_GL_modern_sites$PI, na.rm = T)
length(na.omit(nuc_div_GL_modern_sites$PI))
mean(nuc_div_GL_hist_sites$PI, na.rm = T)
sd(nuc_div_GL_hist_sites$PI, na.rm = T)
length(na.omit(nuc_div_GL_hist_sites$PI))
mean(nuc_div_IS_modern_sites$PI, na.rm = T)
sd(nuc_div_IS_modern_sites$PI, na.rm = T)
length(na.omit(nuc_div_IS_modern_sites$PI))
mean(nuc_div_IS_hist_sites$PI, na.rm = T)
sd(nuc_div_IS_hist_sites$PI, na.rm = T)
length(na.omit(nuc_div_IS_hist_sites$PI))
mean(nuc_div_NO_modern_sites$PI, na.rm = T)
sd(nuc_div_NO_modern_sites$PI, na.rm = T)
length(na.omit(nuc_div_NO_modern_sites$PI))
mean(nuc_div_NO_hist_sites$PI, na.rm = T)
sd(nuc_div_NO_hist_sites$PI, na.rm = T)
length(na.omit(nuc_div_NO_hist_sites$PI))
mean(nuc_div_DK_modern_sites$PI, na.rm = T)
sd(nuc_div_DK_modern_sites$PI, na.rm = T)
length(na.omit(nuc_div_DK_modern_sites$PI))
mean(nuc_div_DK_hist_sites$PI, na.rm = T)
sd(nuc_div_DK_hist_sites$PI, na.rm = T)
length(na.omit(nuc_div_DK_hist_sites$PI))
mean(nuc_div_EE_modern_sites$PI, na.rm = T)
sd(nuc_div_EE_modern_sites$PI, na.rm = T)
length(na.omit(nuc_div_EE_modern_sites$PI))
mean(nuc_div_TU_hist_sites$PI, na.rm = T)
sd(nuc_div_TU_hist_sites$PI, na.rm = T)
length(na.omit(nuc_div_TU_hist_sites$PI))

####################################################
############################ tajimasD 

tajD_all<-read.table("tajimasD/tajimasD_win50k_120321_all.Tajima.D", header = T)
tajD_GL_modern<-read.table("tajimasD/tajimasD_win50k_modern_120321_GL.Tajima.D", header = T)
tajD_GL_hist<-read.table("tajimasD/tajimasD_win50k_historic_120321_GL.Tajima.D", header = T)
tajD_IS_modern<-read.table("tajimasD/tajimasD_win50k_modern_120321_IS.Tajima.D", header = T)
tajD_IS_hist<-read.table("tajimasD/tajimasD_win50k_historic_120321_IS.Tajima.D", header = T)
tajD_NO_modern<-read.table("tajimasD/tajimasD_win50k_modern_120321_NO.Tajima.D", header = T)
tajD_NO_hist<-read.table("tajimasD/tajimasD_win50k_historic_120321_NO.Tajima.D", header = T)
tajD_DK_modern<-read.table("tajimasD/tajimasD_win50k_modern_120321_DK.Tajima.D", header = T)
tajD_DK_hist<-read.table("tajimasD/tajimasD_win50k_historic_120321_DK.Tajima.D", header = T)
tajD_EE_modern<-read.table("tajimasD/tajimasD_win50k_modern_120321_EE.Tajima.D", header = T)
tajD_TU_hist<-read.table("tajimasD/tajimasD_win50k_historic_120321_TU.Tajima.D", header = T)



sum(na.omit(tajD_GL_modern$N_SNPS))
sum(na.omit(tajD_GL_hist$N_SNPS))
sum(na.omit(tajD_IS_modern$N_SNPS))
sum(na.omit(tajD_IS_hist$N_SNPS))
sum(na.omit(tajD_NO_modern$N_SNPS))
sum(na.omit(tajD_NO_hist$N_SNPS))
sum(na.omit(tajD_DK_modern$N_SNPS))
sum(na.omit(tajD_DK_hist$N_SNPS))
sum(na.omit(tajD_EE_modern$N_SNPS))
sum(na.omit(tajD_TU_hist$N_SNPS))
sum(na.omit(tajD_all$N_SNPS))

mean(na.omit(tajD_all$TajimaD))
sd(na.omit(tajD_all$TajimaD))
mean(na.omit(tajD_GL_modern$TajimaD))
sd(na.omit(tajD_GL_modern$TajimaD))
mean(na.omit(tajD_GL_hist$TajimaD))
sd(na.omit(tajD_GL_hist$TajimaD))
mean(na.omit(tajD_IS_modern$TajimaD))
sd(na.omit(tajD_IS_modern$TajimaD))
mean(na.omit(tajD_IS_hist$TajimaD))
sd(na.omit(tajD_IS_hist$TajimaD))
mean(na.omit(tajD_NO_modern$TajimaD))
sd(na.omit(tajD_NO_modern$TajimaD))
mean(na.omit(tajD_NO_hist$TajimaD))
sd(na.omit(tajD_NO_hist$TajimaD))
mean(na.omit(tajD_DK_modern$TajimaD))
sd(na.omit(tajD_DK_modern$TajimaD))
mean(na.omit(tajD_DK_hist$TajimaD))
sd(na.omit(tajD_DK_hist$TajimaD))
mean(na.omit(tajD_EE_modern$TajimaD))
sd(na.omit(tajD_EE_modern$TajimaD))
mean(na.omit(tajD_TU_hist$TajimaD))
sd(na.omit(tajD_TU_hist$TajimaD))

##### Gst and heterozygosity with vcfR 
vcf_92_0.05<-read.vcfR("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05.recode.vcf")
inds_vcf_92_0.05<-colnames(vcf_92_0.05@gt)[2:93]
pops_vcf_92_0.05<-colnames(vcf_92_0.05@gt)[2:93]
pops_vcf_92_0.05[1:25]='IS_C'
pops_vcf_92_0.05[26:37]='NO_C'
pops_vcf_92_0.05[38:42]='DK_H'
pops_vcf_92_0.05[43:53]='DK_C'
pops_vcf_92_0.05[54:56]='EE_C'
pops_vcf_92_0.05[57:68]='GL_C'
pops_vcf_92_0.05[69:76]='GL_H'
pops_vcf_92_0.05[77:78]='IS_H'
pops_vcf_92_0.05[79:91]='NO_H'
pops_vcf_92_0.05[92]='TU_H'
pops_vcf_92_0.05_fac<-as.factor(pops_vcf_92_0.05)
genetic_diff_92<-genetic_diff(vcf_92_0.05, pops=pops_vcf_92_0.05_fac, method = "nei")

popsum_GL_C<-gt.to.popsum(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05_fac=="GL_C"]+1)])
popsum_GL_H<-gt.to.popsum(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05_fac=="GL_H"]+1)])
popsum_IS_C<-gt.to.popsum(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05_fac=="IS_C"]+1)])
popsum_IS_H<-gt.to.popsum(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05_fac=="IS_H"]+1)])
popsum_NO_C<-gt.to.popsum(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05_fac=="NO_C"]+1)])
popsum_NO_H<-gt.to.popsum(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05_fac=="NO_H"]+1)])
popsum_DK_C<-gt.to.popsum(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05_fac=="DK_C"]+1)])
popsum_DK_H<-gt.to.popsum(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05_fac=="DK_H"]+1)])
popsum_EE_C<-gt.to.popsum(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05_fac=="EE_C"]+1)])
popsum_TU_H<-gt.to.popsum(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05_fac=="TU_H"]+1)])


mean(na.omit(popsum_GL_C$He/length(na.omit(popsum_GL_C$He))))
mean(na.omit(popsum_GL_H$He/length(na.omit(popsum_GL_H$He))))
mean(na.omit(popsum_IS_C$He/length(na.omit(popsum_IS_C$He))))
mean(na.omit(popsum_IS_H$He/length(na.omit(popsum_IS_H$He))))
mean(na.omit(popsum_NO_C$He/length(na.omit(popsum_NO_C$He))))
mean(na.omit(popsum_NO_H$He/length(na.omit(popsum_NO_H$He))))
mean(na.omit(popsum_DK_C$He/length(na.omit(popsum_DK_C$He))))
mean(na.omit(popsum_DK_H$He/length(na.omit(popsum_DK_H$He))))
mean(na.omit(popsum_EE_C$He/length(na.omit(popsum_EE_C$He))))
mean(na.omit(popsum_TU_H$He/length(na.omit(popsum_TU_H$He))))

sd(na.omit(popsum_GL_C$He/length(na.omit(popsum_GL_C$He))))
sd(na.omit(popsum_GL_H$He/length(na.omit(popsum_GL_H$He))))
sd(na.omit(popsum_IS_C$He/length(na.omit(popsum_IS_C$He))))
sd(na.omit(popsum_IS_H$He/length(na.omit(popsum_IS_H$He))))
sd(na.omit(popsum_NO_C$He/length(na.omit(popsum_NO_C$He))))
sd(na.omit(popsum_NO_H$He/length(na.omit(popsum_NO_H$He))))
sd(na.omit(popsum_DK_C$He/length(na.omit(popsum_DK_C$He))))
sd(na.omit(popsum_DK_H$He/length(na.omit(popsum_DK_H$He))))
sd(na.omit(popsum_EE_C$He/length(na.omit(popsum_EE_C$He))))
sd(na.omit(popsum_TU_H$He/length(na.omit(popsum_TU_H$He))))

## Gst perhaps with p-val new
x <- vcfR2genlight(vcf_92_0.05)
dist(x)-> x.dist
xdist.mat= as.matrix(x.dist)
pop_dist<-as.vector(read.table("Populationlist_perind.txt"))

dist.perm.mat=function(mat,pop,i,j,perm)
{
  # calulcates net distance between two groups
  # mat is a distance matrix
  # pop is a vector with two pop indexes
  # perm is a number of permutations
  # byr til töflu af 2 synum ij
  mat=mat[(1:length(pop))[pop==i|pop==j],(1:length(pop))[pop==i|pop==j]]
  pop=pop[pop==i|pop==j]
  pop.i=(1:length(pop))[pop==i]
  pop.j=(1:length(pop))[pop==j]
  mat.ii=mat[pop.i,pop.i]
  mat.jj=mat[pop.j,pop.j]
  mat.ij=mat[pop.i,pop.j]
  m.ii=mean(mat.ii[lower.tri(mat.ii)])
  m.jj=mean(mat.jj[lower.tri(mat.jj)])
  m.ij=mean(unlist(mat.ij))
  obs.diff=abs(m.ij-(m.ii+m.jj)/2)
  perm.diff=1:perm
  for(k in 1:perm)
  {
    pop.index=sample(1:nrow(mat))
    mat=mat[pop.index,pop.index]
    mat.ii=mat[pop.i,pop.i]
    mat.jj=mat[pop.j,pop.j]
    mat.ij=mat[pop.i,pop.j]
    m.ii=mean(mat.ii[lower.tri(mat.ii)])
    m.jj=mean(mat.jj[lower.tri(mat.jj)])
    m.ij=mean(unlist(mat.ij))
    perm.diff[k]=abs(m.ij-(m.ii+m.jj)/2)
  }
  p.obs=length(perm.diff[perm.diff>=obs.diff])/perm
  return(c("D.obs"=obs.diff,"P"=p.obs))
}

dist.perm.mat(xdist.mat,pop_dist,"IS_C","IS_H",100)
dist.perm.mat(xdist.mat,pop_dist,"IS_C","GL_C",100)
dist.perm.mat(xdist.mat,pop_dist,"IS_C","NO_C",100)


for(i in 1:10){for(j in (i+1):10){
  tmp=dist.perm.mat(xdist.mat,pop,i,j,100);Dxy.mat[i,j]=tmp[1];Dxy.mat[j,i]=
    tmp[2]}}

library(vegan)
adonis2(xdist.mat~pop_dist) 

cont_dist_020821<-read.table("dist_contemporary_snæbjörn_020821.txt")
pairghist.cont_dist_020821<-cmdscale(as.dist(t(cont_dist_020821)))
plot(pairghist.cont_dist_020821)

Dxy.mat<-read.table("Dist_above_to_compare_with_Gst_050821.txt", header = T)
gst.mat<-read.table("Gst_below_to_compare_with_dist_050821.txt", header = T)
t(Dxy.mat) -> Dxy.mat.t # to have the distances in the lower triangle
Dxy.mat.t[lower.tri(Dxy.mat.t)]
gst.mat[lower.tri(gst.mat)]
plot(Dxy.mat.t[lower.tri(Dxy.mat.t)],gst.mat[lower.tri(gst.mat)])

dist_mat<-read.table("Dist_above_to_compare_with_Gst_050821_diss.txt",fill = T, header = T)
gst_mat<-read.table("Gst_below_to_compare_with_dist_050821_diss.txt", fill = T, header = T)
as.dist(dist_mat)

# GL_C others
# genedif_GL_C_GL_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="GL_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_C"|
#                                                              pops_vcf_92_0.05=="GL_H"]),
#                                 method = "nei")
# 
# genedif_GL_C_IS_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="IS_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_C"|
#                                                              pops_vcf_92_0.05=="IS_C"]),
#                                 method = "nei")
# 
# genedif_GL_C_IS_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="IS_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_C"|
#                                                              pops_vcf_92_0.05=="IS_H"]),
#                                 method = "nei")
# 
# genedif_GL_C_NO_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="NO_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_C"|
#                                                              pops_vcf_92_0.05=="NO_C"]),
#                                 method = "nei")
# 
# genedif_GL_C_NO_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="NO_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_C"|
#                                                              pops_vcf_92_0.05=="NO_H"]),
#                                 method = "nei")
# 
# genedif_GL_C_DK_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="DK_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_C"|
#                                                              pops_vcf_92_0.05=="DK_C"]),
#                                 method = "nei")
# 
# genedif_GL_C_DK_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="DK_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_C"|
#                                                              pops_vcf_92_0.05=="DK_H"]),
#                                 method = "nei")
# 
# genedif_GL_C_EE_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="EE_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_C"|
#                                                              pops_vcf_92_0.05=="EE_C"]),
#                                 method = "nei")
# 
# genedif_GL_C_TU_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="TU_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_C"|
#                                                              pops_vcf_92_0.05=="TU_H"]),
#                                 method = "nei")
# 
# #GL_H others
# genedif_GL_H_IS_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="IS_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_H"|
#                                                              pops_vcf_92_0.05=="IS_C"]),
#                                 method = "nei")
# 
# genedif_GL_H_IS_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="IS_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_H"|
#                                                              pops_vcf_92_0.05=="IS_H"]),
#                                 method = "nei")
# 
# genedif_GL_H_NO_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="NO_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_H"|
#                                                              pops_vcf_92_0.05=="NO_C"]),
#                                 method = "nei")
# 
# genedif_GL_H_NO_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="NO_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_H"|
#                                                              pops_vcf_92_0.05=="NO_H"]),
#                                 method = "nei")
# 
# genedif_GL_H_DK_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="DK_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_H"|
#                                                              pops_vcf_92_0.05=="DK_C"]),
#                                 method = "nei")
# 
# genedif_GL_H_DK_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="DK_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_H"|
#                                                              pops_vcf_92_0.05=="DK_H"]),
#                                 method = "nei")
# 
# genedif_GL_H_EE_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="EE_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_H"|
#                                                              pops_vcf_92_0.05=="EE_C"]),
#                                 method = "nei")
# 
# genedif_GL_H_TU_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="GL_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="TU_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="GL_H"|
#                                                              pops_vcf_92_0.05=="TU_H"]),
#                                 method = "nei")
# #
# # IS_C others
# genedif_IS_C_IS_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="IS_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="IS_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="IS_C"|
#                                                              pops_vcf_92_0.05=="IS_H"]),
#                                 method = "nei")
# 
# genedif_IS_C_NO_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="IS_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="NO_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="IS_C"|
#                                                              pops_vcf_92_0.05=="NO_C"]),
#                                 method = "nei")
# 
# genedif_IS_C_NO_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="IS_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="IS_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="IS_C"|
#                                                              pops_vcf_92_0.05=="IS_H"]),
#                                 method = "nei")
# 
# genedif_IS_C_DK_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="IS_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="DK_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="IS_C"|
#                                                              pops_vcf_92_0.05=="DK_C"]),
#                                 method = "nei")
# 
# genedif_IS_C_DK_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="IS_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="DK_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="IS_C"|
#                                                              pops_vcf_92_0.05=="DK_H"]),
#                                 method = "nei")
# 
# genedif_IS_C_EE_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="IS_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="EE_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="IS_C"|
#                                                              pops_vcf_92_0.05=="EE_C"]),
#                                 method = "nei")
# 
# genedif_IS_C_TU_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="IS_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="TU_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="IS_C"|
#                                                              pops_vcf_92_0.05=="TU_H"]),
#                                 method = "nei")
# 
# # IS_H and others
# genedif_IS_H_NO_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="IS_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="NO_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="IS_H"|
#                                                              pops_vcf_92_0.05=="NO_C"]),
#                                 method = "nei")
# 
# genedif_IS_H_NO_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="IS_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="NO_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="IS_H"|
#                                                              pops_vcf_92_0.05=="NO_H"]),
#                                 method = "nei")
# 
# genedif_IS_H_DK_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="IS_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="DK_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="IS_H"|
#                                                              pops_vcf_92_0.05=="DK_C"]),
#                                 method = "nei")
# 
# genedif_IS_H_DK_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="IS_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="DK_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="IS_H"|
#                                                              pops_vcf_92_0.05=="DK_H"]),
#                                 method = "nei")
# 
# genedif_IS_H_EE_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="IS_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="EE_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="IS_H"|
#                                                              pops_vcf_92_0.05=="EE_C"]),
#                                 method = "nei")
# 
# genedif_IS_H_TU_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="IS_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="TU_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="IS_H"|
#                                                              pops_vcf_92_0.05=="TU_H"]),
#                                 method = "nei")
# 
# # #NO_C and others
# genedif_NO_C_NO_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="NO_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="NO_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="NO_C"|
#                                                              pops_vcf_92_0.05=="NO_H"]),
#                                 method = "nei")
# 
# genedif_NO_C_DK_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="NO_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="DK_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="NO_C"|
#                                                              pops_vcf_92_0.05=="DK_C"]),
#                                 method = "nei")
# 
# genedif_NO_C_DK_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="NO_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="DK_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="NO_C"|
#                                                              pops_vcf_92_0.05=="DK_H"]),
#                                 method = "nei")
# 
# genedif_NO_C_EE_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="NO_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="EE_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="NO_C"|
#                                                              pops_vcf_92_0.05=="EE_C"]),
#                                 method = "nei")
# 
# genedif_NO_C_TU_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="NO_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="TU_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="NO_C"|
#                                                              pops_vcf_92_0.05=="TU_H"]),
#                                 method = "nei")
# 
# #NO_H and other
# genedif_NO_H_DK_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="NO_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="DK_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="NO_H"|
#                                                              pops_vcf_92_0.05=="DK_C"]),
#                                 method = "nei")
# 
# genedif_NO_H_DK_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="NO_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="DK_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="NO_H"|
#                                                              pops_vcf_92_0.05=="DK_H"]),
#                                 method = "nei")
# 
# genedif_NO_H_EE_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="NO_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="EE_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="NO_H"|
#                                                              pops_vcf_92_0.05=="EE_C"]),
#                                 method = "nei")
# 
# genedif_NO_H_TU_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="NO_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="TU_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="NO_H"|
#                                                              pops_vcf_92_0.05=="TU_H"]),
#                                 method = "nei")
# 
# ##DK_C and others
# genedif_DK_C_DK_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="DK_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="DK_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="DK_C"|
#                                                              pops_vcf_92_0.05=="DK_H"]),
#                                 method = "nei")
# 
# genedif_DK_C_EE_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="DK_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="EE_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="DK_C"|
#                                                              pops_vcf_92_0.05=="EE_C"]),
#                                 method = "nei")
# 
# genedif_DK_C_TU_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="DK_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="TU_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="DK_C"|
#                                                              pops_vcf_92_0.05=="TU_H"]),
#                                 method = "nei")
# 
# #DK_H and other
# genedif_DK_H_EE_C<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="DK_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="EE_C"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="DK_H"|
#                                                              pops_vcf_92_0.05=="EE_C"]),
#                                 method = "nei")
# 
# genedif_DK_H_TU_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="DK_H"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="TU_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="DK_H"|
#                                                              pops_vcf_92_0.05=="TU_H"]),
#                                 method = "nei")
# 
# #EE_C and TU_H
# genedif_EE_C_TU_H<-genetic_diff(vcf_92_0.05[,c(1,(1:92)[pops_vcf_92_0.05=="EE_C"]+1,
#                                                (1:92)[pops_vcf_92_0.05=="TU_H"]+1)],
#                                 as.factor(pops_vcf_92_0.05[pops_vcf_92_0.05=="EE_C"|
#                                                              pops_vcf_92_0.05=="TU_H"]),
#                                 method = "nei")
# #GL_C
# write.table(genedif_GL_C_GL_H, "Gst_vcfR/Gst_genedif_GL_C_GL_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_C_IS_C, "Gst_vcfR/Gst_genedif_GL_C_IS_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_C_IS_H, "Gst_vcfR/Gst_genedif_GL_C_IS_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_C_NO_C, "Gst_vcfR/Gst_genedif_GL_C_NO_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_C_NO_H, "Gst_vcfR/Gst_genedif_GL_C_NO_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_C_DK_C, "Gst_vcfR/Gst_genedif_GL_C_DK_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_C_DK_H, "Gst_vcfR/Gst_genedif_GL_C_DK_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_C_EE_C, "Gst_vcfR/Gst_genedif_GL_C_EE_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_C_TU_H, "Gst_vcfR/Gst_genedif_GL_C_TU_H.txt", quote = F, row.names = F, sep = "\t")
# 
# #GL_H
# write.table(genedif_GL_H_IS_C, "Gst_vcfR/Gst_genedif_GL_H_IS_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_H_IS_H, "Gst_vcfR/Gst_genedif_GL_H_IS_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_H_NO_C, "Gst_vcfR/Gst_genedif_GL_H_NO_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_H_NO_H, "Gst_vcfR/Gst_genedif_GL_H_NO_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_H_DK_C, "Gst_vcfR/Gst_genedif_GL_H_DK_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_H_DK_H, "Gst_vcfR/Gst_genedif_GL_H_DK_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_H_EE_C, "Gst_vcfR/Gst_genedif_GL_H_EE_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_GL_H_TU_H, "Gst_vcfR/Gst_genedif_GL_H_TU_H.txt", quote = F, row.names = F, sep = "\t")
# 
# # IS_C others
# write.table(genedif_IS_C_IS_H, "Gst_vcfR/Gst_genedif_IS_C_IS_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_IS_C_NO_C, "Gst_vcfR/Gst_genedif_IS_C_NO_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_IS_C_NO_H, "Gst_vcfR/Gst_genedif_IS_C_NO_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_IS_C_DK_C, "Gst_vcfR/Gst_genedif_IS_C_DK_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_IS_C_DK_H, "Gst_vcfR/Gst_genedif_IS_C_DK_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_IS_C_EE_C, "Gst_vcfR/Gst_genedif_IS_C_EE_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_IS_C_TU_H, "Gst_vcfR/Gst_genedif_IS_C_TU_H.txt", quote = F, row.names = F, sep = "\t")
# 
# # IS_H, "Gst_vcfR/Gst_ and others
# write.table(genedif_IS_H_NO_C, "Gst_vcfR/Gst_genedif_IS_H_NO_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_IS_H_NO_H, "Gst_vcfR/Gst_genedif_IS_H_NO_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_IS_H_DK_C, "Gst_vcfR/Gst_genedif_IS_H_DK_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_IS_H_DK_H, "Gst_vcfR/Gst_genedif_IS_H_DK_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_IS_H_EE_C, "Gst_vcfR/Gst_genedif_IS_H_EE_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_IS_H_TU_H, "Gst_vcfR/Gst_genedif_IS_H_TU_H.txt", quote = F, row.names = F, sep = "\t")
# 
# # #NO_C and others
# write.table(genedif_NO_C_NO_H, "Gst_vcfR/Gst_genedif_NO_C_NO_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_NO_C_DK_C, "Gst_vcfR/Gst_genedif_NO_C_DK_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_NO_C_DK_H, "Gst_vcfR/Gst_genedif_NO_C_DK_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_NO_C_EE_C, "Gst_vcfR/Gst_genedif_NO_C_EE_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_NO_C_TU_H, "Gst_vcfR/Gst_genedif_NO_C_TU_H.txt", quote = F, row.names = F, sep = "\t")
# 
# #NO_H, "Gst_vcfR/Gst_ and other
# write.table(genedif_NO_H_DK_C, "Gst_vcfR/Gst_genedif_NO_H_DK_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_NO_H_DK_H, "Gst_vcfR/Gst_genedif_NO_H_DK_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_NO_H_EE_C, "Gst_vcfR/Gst_genedif_NO_H_EE_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_NO_H_TU_H, "Gst_vcfR/Gst_genedif_NO_H_TU_H.txt", quote = F, row.names = F, sep = "\t")
# 
# ##DK_C and others
# write.table(genedif_DK_C_DK_H, "Gst_vcfR/Gst_genedif_DK_C_DK_H.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_DK_C_EE_C, "Gst_vcfR/Gst_genedif_DK_C_EE_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_DK_C_TU_H, "Gst_vcfR/Gst_genedif_DK_C_TU_H.txt", quote = F, row.names = F, sep = "\t")
# 
# #DK_H, "Gst_vcfR/Gst_ and other
# write.table(genedif_DK_H_EE_C, "Gst_vcfR/Gst_genedif_DK_H_EE_C.txt", quote = F, row.names = F, sep = "\t")
# write.table(genedif_DK_H_TU_H, "Gst_vcfR/Gst_genedif_DK_H_TU_H.txt", quote = F, row.names = F, sep = "\t")
# 
# #EE_C and TU_H, "Gst_vcfR/Gst_
# write.table(genedif_EE_C_TU_H, "Gst_vcfR/Gst_genedif_EE_C_TU_H.txt", quote = F, row.names = F, sep = "\t")

# # #GL_C and others
# mean(genedif_GL_C_GL_H$Gst, na.rm = T)
# quantile(genedif_GL_C_GL_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_C_GL_H$Gst, na.rm = T)
# mean(genedif_GL_C_IS_C$Gst, na.rm = T)
# quantile(genedif_GL_C_IS_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_C_IS_C$Gst, na.rm = T)
# mean(genedif_GL_C_IS_H$Gst, na.rm = T)
# quantile(genedif_GL_C_IS_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_C_IS_H$Gst, na.rm = T)
# mean(genedif_GL_C_NO_C$Gst, na.rm = T)
# quantile(genedif_GL_C_NO_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_C_NO_C$Gst, na.rm = T)
# mean(genedif_GL_C_NO_H$Gst, na.rm = T)
# quantile(genedif_GL_C_NO_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_C_NO_H$Gst, na.rm = T)
# mean(genedif_GL_C_DK_C$Gst, na.rm = T)
# quantile(genedif_GL_C_DK_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_C_DK_C$Gst, na.rm = T)
# mean(genedif_GL_C_DK_H$Gst, na.rm = T)
# quantile(genedif_GL_C_DK_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_C_DK_H$Gst, na.rm = T)
mean(genedif_GL_C_EE_C$Gst, na.rm = T)
quantile(genedif_GL_C_EE_C$Gst, c(0.025, 0.975),na.rm=T)
sd(genedif_GL_C_EE_C$Gst, na.rm = T)
# mean(genedif_GL_C_TU_H$Gst, na.rm = T)
# quantile(genedif_GL_C_TU_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_C_TU_H$Gst, na.rm = T)
# 
# #GL_H and others
# mean(genedif_GL_H_IS_C$Gst, na.rm = T)
# quantile(genedif_GL_H_IS_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_H_IS_C$Gst, na.rm = T)
# mean(genedif_GL_H_IS_H$Gst, na.rm = T)
# quantile(genedif_GL_H_IS_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_H_IS_H$Gst, na.rm = T)
# mean(genedif_GL_H_NO_C$Gst, na.rm = T)
# quantile(genedif_GL_H_NO_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_H_NO_C$Gst, na.rm = T)
# mean(genedif_GL_H_NO_H$Gst, na.rm = T)
# quantile(genedif_GL_H_NO_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_H_NO_H$Gst, na.rm = T)
# mean(genedif_GL_H_DK_C$Gst, na.rm = T)
# quantile(genedif_GL_H_DK_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_H_DK_C$Gst, na.rm = T)
# mean(genedif_GL_H_DK_H$Gst, na.rm = T)
# quantile(genedif_GL_H_DK_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_H_DK_H$Gst, na.rm = T)
# mean(genedif_GL_H_EE_C$Gst, na.rm = T)
# quantile(genedif_GL_H_EE_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_H_EE_C$Gst, na.rm = T)
# mean(genedif_GL_H_TU_H$Gst, na.rm = T)
# quantile(genedif_GL_H_TU_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_GL_H_TU_H$Gst, na.rm = T)
# 
# #IS_C and others
# mean(genedif_IS_C_IS_H$Gst, na.rm = T)
# quantile(genedif_IS_C_IS_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_IS_C_IS_H$Gst, na.rm = T)
# mean(genedif_IS_C_NO_C$Gst, na.rm = T)
# quantile(genedif_IS_C_NO_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_IS_C_NO_C$Gst, na.rm = T)
# mean(genedif_IS_C_NO_H$Gst, na.rm = T)
# quantile(genedif_IS_C_NO_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_IS_C_NO_H$Gst, na.rm = T)
# mean(genedif_IS_C_DK_C$Gst, na.rm = T)
# quantile(genedif_IS_C_DK_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_IS_C_DK_C$Gst, na.rm = T)
# mean(genedif_IS_C_DK_H$Gst, na.rm = T)
# quantile(genedif_IS_C_DK_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_IS_C_DK_H$Gst, na.rm = T)
# mean(genedif_IS_C_EE_C$Gst, na.rm = T)
# quantile(genedif_IS_C_EE_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_IS_C_EE_C$Gst, na.rm = T)
# mean(genedif_IS_C_TU_H$Gst, na.rm = T)
# quantile(genedif_IS_C_TU_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_IS_C_TU_H$Gst, na.rm = T)
# 
# #IS_H and others
# mean(genedif_IS_H_NO_C$Gst, na.rm = T)
# quantile(genedif_IS_H_NO_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_IS_H_NO_C$Gst, na.rm = T)
# mean(genedif_IS_H_NO_H$Gst, na.rm = T)
# quantile(genedif_IS_H_NO_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_IS_H_NO_H$Gst, na.rm = T)
# mean(genedif_IS_H_DK_C$Gst, na.rm = T)
# quantile(genedif_IS_H_DK_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_IS_H_DK_C$Gst, na.rm = T)
# mean(genedif_IS_H_DK_H$Gst, na.rm = T)
# quantile(genedif_IS_H_DK_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_IS_H_DK_H$Gst, na.rm = T)
# mean(genedif_IS_H_EE_C$Gst, na.rm = T)
# quantile(genedif_IS_H_EE_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_IS_H_EE_C$Gst, na.rm = T)
# mean(genedif_IS_H_TU_H$Gst, na.rm = T)
# quantile(genedif_IS_H_TU_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_IS_H_TU_H$Gst, na.rm = T)
# 
# #NO_C and others
# mean(genedif_NO_C_NO_H$Gst, na.rm = T)
# quantile(genedif_NO_C_NO_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_NO_C_NO_H$Gst, na.rm = T)
# mean(genedif_NO_C_DK_C$Gst, na.rm = T)
# quantile(genedif_NO_C_DK_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_NO_C_DK_C$Gst, na.rm = T)
# mean(genedif_NO_C_DK_H$Gst, na.rm = T)
# quantile(genedif_NO_C_DK_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_NO_C_DK_H$Gst, na.rm = T)
# mean(genedif_NO_C_EE_C$Gst, na.rm = T)
# quantile(genedif_NO_C_EE_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_NO_C_EE_C$Gst, na.rm = T)
# mean(genedif_NO_C_TU_H$Gst, na.rm = T)
# quantile(genedif_NO_C_TU_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_NO_C_TU_H$Gst, na.rm = T)
# 
# #NO_H and others
# mean(genedif_NO_H_DK_C$Gst, na.rm = T)
# quantile(genedif_NO_H_DK_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_NO_H_DK_C$Gst, na.rm = T)
# mean(genedif_NO_H_DK_H$Gst, na.rm = T)
# quantile(genedif_NO_H_DK_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_NO_H_DK_H$Gst, na.rm = T)
# mean(genedif_NO_H_EE_C$Gst, na.rm = T)
# quantile(genedif_NO_H_EE_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_NO_H_EE_C$Gst, na.rm = T)
# mean(genedif_NO_H_TU_H$Gst, na.rm = T)
# quantile(genedif_NO_H_TU_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_NO_H_TU_H$Gst, na.rm = T)
# 
# #DK_C
# mean(genedif_DK_C_DK_H$Gst, na.rm = T)
# quantile(genedif_DK_C_DK_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_DK_C_DK_H$Gst, na.rm = T)
# mean(genedif_DK_C_EE_C$Gst, na.rm = T)
# quantile(genedif_DK_C_EE_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_DK_C_EE_C$Gst, na.rm = T)
# mean(genedif_DK_C_TU_H$Gst, na.rm = T)
# quantile(genedif_DK_C_TU_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_DK_C_TU_H$Gst, na.rm = T)
# 
# #DK_H
# mean(genedif_DK_H_EE_C$Gst, na.rm = T)
# quantile(genedif_DK_H_EE_C$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_DK_H_EE_C$Gst, na.rm = T)
# mean(genedif_DK_H_TU_H$Gst, na.rm = T)
# quantile(genedif_DK_H_TU_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_DK_H_TU_H$Gst, na.rm = T)
# 
# #EE_C and TU_H
# mean(genedif_EE_C_TU_H$Gst, na.rm = T)
# quantile(genedif_EE_C_TU_H$Gst, c(0.025, 0.975),na.rm=T)
# sd(genedif_EE_C_TU_H$Gst, na.rm = T)
# 
# 
# ## GL_C
# sd(genedif_GL_C_GL_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_C_GL_H$Gst)))
# sd(genedif_GL_C_IS_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_C_IS_C$Gst)))
# sd(genedif_GL_C_IS_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_C_IS_H$Gst)))
# sd(genedif_GL_C_NO_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_C_NO_C$Gst)))
# sd(genedif_GL_C_NO_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_C_NO_H$Gst)))
# sd(genedif_GL_C_DK_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_C_DK_C$Gst)))
# sd(genedif_GL_C_DK_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_C_DK_H$Gst)))
# sd(genedif_GL_C_EE_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_C_EE_C$Gst)))
# sd(genedif_GL_C_TU_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_C_TU_H$Gst)))
# 
# ## GL_H
# sd(genedif_GL_H_IS_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_H_IS_C$Gst)))
# sd(genedif_GL_H_IS_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_H_IS_H$Gst)))
# sd(genedif_GL_H_NO_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_H_NO_C$Gst)))
# sd(genedif_GL_H_NO_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_H_NO_H$Gst)))
# sd(genedif_GL_H_DK_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_H_DK_C$Gst)))
# sd(genedif_GL_H_DK_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_H_DK_H$Gst)))
# sd(genedif_GL_H_EE_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_H_EE_C$Gst)))
# sd(genedif_GL_H_TU_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_GL_H_TU_H$Gst)))
# 
# ## IS_C
# sd(genedif_IS_C_IS_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_IS_C_IS_H$Gst)))
# sd(genedif_IS_C_NO_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_IS_C_NO_C$Gst)))
# sd(genedif_IS_C_NO_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_IS_C_NO_H$Gst)))
# sd(genedif_IS_C_DK_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_IS_C_DK_C$Gst)))
# sd(genedif_IS_C_DK_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_IS_C_DK_H$Gst)))
# sd(genedif_IS_C_EE_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_IS_C_EE_C$Gst)))
# sd(genedif_IS_C_TU_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_IS_C_TU_H$Gst)))
# 
# ## IS_H
# sd(genedif_IS_H_NO_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_IS_H_NO_C$Gst)))
# sd(genedif_IS_H_NO_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_IS_H_NO_H$Gst)))
# sd(genedif_IS_H_DK_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_IS_H_DK_C$Gst)))
# sd(genedif_IS_H_DK_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_IS_H_DK_H$Gst)))
# sd(genedif_IS_H_EE_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_IS_H_EE_C$Gst)))
# sd(genedif_IS_H_TU_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_IS_H_TU_H$Gst)))
# 
# ## NO_C
# sd(genedif_NO_C_NO_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_NO_C_NO_H$Gst)))
# sd(genedif_NO_C_DK_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_NO_C_DK_C$Gst)))
# sd(genedif_NO_C_DK_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_NO_C_DK_H$Gst)))
# sd(genedif_NO_C_EE_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_NO_C_EE_C$Gst)))
# sd(genedif_NO_C_TU_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_NO_C_TU_H$Gst)))
# 
# ## NO_H
# sd(genedif_NO_H_DK_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_NO_H_DK_C$Gst)))
# sd(genedif_NO_H_DK_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_NO_H_DK_H$Gst)))
# sd(genedif_NO_H_EE_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_NO_H_EE_C$Gst)))
# sd(genedif_NO_H_TU_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_NO_H_TU_H$Gst)))
# 
# ## DK_C
# sd(genedif_DK_C_DK_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_DK_C_DK_H$Gst)))
# sd(genedif_DK_C_EE_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_DK_C_EE_C$Gst)))
# sd(genedif_DK_C_TU_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_DK_C_TU_H$Gst)))
# 
# ## DK_H
# sd(genedif_DK_H_EE_C$Gst, na.rm = T)/sqrt(na.omit(length(genedif_DK_H_EE_C$Gst)))
# sd(genedif_DK_H_TU_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_DK_H_TU_H$Gst)))
# 
# ## EE_C
# sd(genedif_EE_C_TU_H$Gst, na.rm = T)/sqrt(na.omit(length(genedif_EE_C_TU_H$Gst)))
# 
# ## GL_C
# length(na.omit(genedif_GL_C_GL_H$Gst[genedif_GL_C_GL_H$Gst>c(quantile(na.omit(genedif_GL_C_GL_H$Gst),0.975)) | 
#                                        genedif_GL_C_GL_H$Gst<c(quantile(na.omit(genedif_GL_C_GL_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_C_GL_H$Gst))
# 
# length(na.omit(genedif_GL_C_IS_C$Gst[genedif_GL_C_IS_C$Gst>c(quantile(na.omit(genedif_GL_C_IS_C$Gst),0.975)) | 
#                                        genedif_GL_C_IS_C$Gst<c(quantile(na.omit(genedif_GL_C_IS_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_C_IS_C$Gst))
# 
# length(na.omit(genedif_GL_C_IS_H$Gst[genedif_GL_C_IS_H$Gst>c(quantile(na.omit(genedif_GL_C_IS_H$Gst),0.975)) | 
#                                        genedif_GL_C_IS_H$Gst<c(quantile(na.omit(genedif_GL_C_IS_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_C_IS_H$Gst))
# 
# length(na.omit(genedif_GL_C_NO_C$Gst[genedif_GL_C_NO_C$Gst>c(quantile(na.omit(genedif_GL_C_NO_C$Gst),0.975)) | 
#                                        genedif_GL_C_NO_C$Gst<c(quantile(na.omit(genedif_GL_C_NO_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_C_NO_C$Gst))
# 
# length(na.omit(genedif_GL_C_NO_H$Gst[genedif_GL_C_NO_H$Gst>c(quantile(na.omit(genedif_GL_C_NO_H$Gst),0.975)) | 
#                                        genedif_GL_C_NO_H$Gst<c(quantile(na.omit(genedif_GL_C_NO_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_C_NO_H$Gst))
# 
# length(na.omit(genedif_GL_C_DK_C$Gst[genedif_GL_C_DK_C$Gst>c(quantile(na.omit(genedif_GL_C_DK_C$Gst),0.975)) | 
#                                        genedif_GL_C_DK_C$Gst<c(quantile(na.omit(genedif_GL_C_DK_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_C_DK_C$Gst))
# 
# length(na.omit(genedif_GL_C_DK_H$Gst[genedif_GL_C_DK_H$Gst>c(quantile(na.omit(genedif_GL_C_DK_H$Gst),0.975)) | 
#                                        genedif_GL_C_DK_H$Gst<c(quantile(na.omit(genedif_GL_C_DK_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_C_DK_H$Gst))
# 
# length(na.omit(genedif_GL_C_EE_C$Gst[genedif_GL_C_EE_C$Gst>c(quantile(na.omit(genedif_GL_C_EE_C$Gst),0.975)) | 
#                                        genedif_GL_C_EE_C$Gst<c(quantile(na.omit(genedif_GL_C_EE_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_C_EE_C$Gst))
# 
# length(na.omit(genedif_GL_C_TU_H$Gst[genedif_GL_C_TU_H$Gst>c(quantile(na.omit(genedif_GL_C_TU_H$Gst),0.975)) | 
#                                        genedif_GL_C_TU_H$Gst<c(quantile(na.omit(genedif_GL_C_TU_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_C_TU_H$Gst))
# 
# ## GL_H
# length(na.omit(genedif_GL_H_IS_C$Gst[genedif_GL_H_IS_C$Gst>c(quantile(na.omit(genedif_GL_H_IS_C$Gst),0.975)) | 
#                                        genedif_GL_H_IS_C$Gst<c(quantile(na.omit(genedif_GL_H_IS_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_H_IS_C$Gst))
# 
# length(na.omit(genedif_GL_H_IS_H$Gst[genedif_GL_H_IS_H$Gst>c(quantile(na.omit(genedif_GL_H_IS_H$Gst),0.975)) | 
#                                        genedif_GL_H_IS_H$Gst<c(quantile(na.omit(genedif_GL_H_IS_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_H_IS_H$Gst))
# 
# length(na.omit(genedif_GL_H_NO_C$Gst[genedif_GL_H_NO_C$Gst>c(quantile(na.omit(genedif_GL_H_NO_C$Gst),0.975)) | 
#                                        genedif_GL_H_NO_C$Gst<c(quantile(na.omit(genedif_GL_H_NO_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_H_NO_C$Gst))
# 
# length(na.omit(genedif_GL_H_NO_H$Gst[genedif_GL_H_NO_H$Gst>c(quantile(na.omit(genedif_GL_H_NO_H$Gst),0.975)) | 
#                                        genedif_GL_H_NO_H$Gst<c(quantile(na.omit(genedif_GL_H_NO_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_H_NO_H$Gst))
# 
# length(na.omit(genedif_GL_H_DK_C$Gst[genedif_GL_H_DK_C$Gst>c(quantile(na.omit(genedif_GL_H_DK_C$Gst),0.975)) | 
#                                        genedif_GL_H_DK_C$Gst<c(quantile(na.omit(genedif_GL_H_DK_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_H_DK_C$Gst))
# 
# length(na.omit(genedif_GL_H_DK_H$Gst[genedif_GL_H_DK_H$Gst>c(quantile(na.omit(genedif_GL_H_DK_H$Gst),0.975)) | 
#                                        genedif_GL_H_DK_H$Gst<c(quantile(na.omit(genedif_GL_H_DK_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_H_DK_H$Gst))
# 
# length(na.omit(genedif_GL_H_EE_C$Gst[genedif_GL_H_EE_C$Gst>c(quantile(na.omit(genedif_GL_H_EE_C$Gst),0.975)) | 
#                                        genedif_GL_H_EE_C$Gst<c(quantile(na.omit(genedif_GL_H_EE_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_H_EE_C$Gst))
# 
# length(na.omit(genedif_GL_H_TU_H$Gst[genedif_GL_H_TU_H$Gst>c(quantile(na.omit(genedif_GL_H_TU_H$Gst),0.975)) | 
#                                        genedif_GL_H_TU_H$Gst<c(quantile(na.omit(genedif_GL_H_TU_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_GL_H_TU_H$Gst))
# 
# ## IS_C
# length(na.omit(genedif_IS_C_IS_H$Gst[genedif_IS_C_IS_H$Gst>c(quantile(na.omit(genedif_IS_C_IS_H$Gst),0.975)) | 
#                                        genedif_IS_C_IS_H$Gst<c(quantile(na.omit(genedif_IS_C_IS_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_IS_C_IS_H$Gst))
# 
# length(na.omit(genedif_IS_C_NO_C$Gst[genedif_IS_C_NO_C$Gst>c(quantile(na.omit(genedif_IS_C_NO_C$Gst),0.975)) | 
#                                        genedif_IS_C_NO_C$Gst<c(quantile(na.omit(genedif_IS_C_NO_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_IS_C_NO_C$Gst))
# 
# length(na.omit(genedif_IS_C_NO_H$Gst[genedif_IS_C_NO_H$Gst>c(quantile(na.omit(genedif_IS_C_NO_H$Gst),0.975)) | 
#                                        genedif_IS_C_NO_H$Gst<c(quantile(na.omit(genedif_IS_C_NO_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_IS_C_NO_H$Gst))
# 
# length(na.omit(genedif_IS_C_DK_C$Gst[genedif_IS_C_DK_C$Gst>c(quantile(na.omit(genedif_IS_C_DK_C$Gst),0.975)) | 
#                                        genedif_IS_C_DK_C$Gst<c(quantile(na.omit(genedif_IS_C_DK_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_IS_C_DK_C$Gst))
# 
# length(na.omit(genedif_IS_C_DK_H$Gst[genedif_IS_C_DK_H$Gst>c(quantile(na.omit(genedif_IS_C_DK_H$Gst),0.975)) | 
#                                        genedif_IS_C_DK_H$Gst<c(quantile(na.omit(genedif_IS_C_DK_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_IS_C_DK_H$Gst))
# 
# length(na.omit(genedif_IS_C_EE_C$Gst[genedif_IS_C_EE_C$Gst>c(quantile(na.omit(genedif_IS_C_EE_C$Gst),0.975)) | 
#                                        genedif_IS_C_EE_C$Gst<c(quantile(na.omit(genedif_IS_C_EE_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_IS_C_EE_C$Gst))
# 
# length(na.omit(genedif_IS_C_TU_H$Gst[genedif_IS_C_TU_H$Gst>c(quantile(na.omit(genedif_IS_C_TU_H$Gst),0.975)) | 
#                                        genedif_IS_C_TU_H$Gst<c(quantile(na.omit(genedif_IS_C_TU_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_IS_C_TU_H$Gst))
# 
# ## IS_H
# length(na.omit(genedif_IS_H_NO_C$Gst[genedif_IS_H_NO_C$Gst>c(quantile(na.omit(genedif_IS_H_NO_C$Gst),0.975)) | 
#                                        genedif_IS_H_NO_C$Gst<c(quantile(na.omit(genedif_IS_H_NO_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_IS_H_NO_C$Gst))
# 
# length(na.omit(genedif_IS_H_NO_H$Gst[genedif_IS_H_NO_H$Gst>c(quantile(na.omit(genedif_IS_H_NO_H$Gst),0.975)) | 
#                                        genedif_IS_H_NO_H$Gst<c(quantile(na.omit(genedif_IS_H_NO_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_IS_H_NO_H$Gst))
# 
# length(na.omit(genedif_IS_H_DK_C$Gst[genedif_IS_H_DK_C$Gst>c(quantile(na.omit(genedif_IS_H_DK_C$Gst),0.975)) | 
#                                        genedif_IS_H_DK_C$Gst<c(quantile(na.omit(genedif_IS_H_DK_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_IS_H_DK_C$Gst))
# 
# length(na.omit(genedif_IS_H_DK_H$Gst[genedif_IS_H_DK_H$Gst>c(quantile(na.omit(genedif_IS_H_DK_H$Gst),0.975)) | 
#                                        genedif_IS_H_DK_H$Gst<c(quantile(na.omit(genedif_IS_H_DK_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_IS_H_DK_H$Gst))
# 
# length(na.omit(genedif_IS_H_EE_C$Gst[genedif_IS_H_EE_C$Gst>c(quantile(na.omit(genedif_IS_H_EE_C$Gst),0.975)) | 
#                                        genedif_IS_H_EE_C$Gst<c(quantile(na.omit(genedif_IS_H_EE_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_IS_H_EE_C$Gst))
# 
# length(na.omit(genedif_IS_H_TU_H$Gst[genedif_IS_H_TU_H$Gst>c(quantile(na.omit(genedif_IS_H_TU_H$Gst),0.975)) | 
#                                        genedif_IS_H_TU_H$Gst<c(quantile(na.omit(genedif_IS_H_TU_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_IS_H_TU_H$Gst))
# 
# ## NO_C
# length(na.omit(genedif_NO_C_NO_H$Gst[genedif_NO_C_NO_H$Gst>c(quantile(na.omit(genedif_NO_C_NO_H$Gst),0.975)) | 
#                                        genedif_NO_C_NO_H$Gst<c(quantile(na.omit(genedif_NO_C_NO_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_NO_C_NO_H$Gst))
# 
# length(na.omit(genedif_NO_C_DK_C$Gst[genedif_NO_C_DK_C$Gst>c(quantile(na.omit(genedif_NO_C_DK_C$Gst),0.975)) | 
#                                        genedif_NO_C_DK_C$Gst<c(quantile(na.omit(genedif_NO_C_DK_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_NO_C_DK_C$Gst))
# 
# length(na.omit(genedif_NO_C_DK_H$Gst[genedif_NO_C_DK_H$Gst>c(quantile(na.omit(genedif_NO_C_DK_H$Gst),0.975)) | 
#                                        genedif_NO_C_DK_H$Gst<c(quantile(na.omit(genedif_NO_C_DK_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_NO_C_DK_H$Gst))
# 
# length(na.omit(genedif_NO_C_EE_C$Gst[genedif_NO_C_EE_C$Gst>c(quantile(na.omit(genedif_NO_C_EE_C$Gst),0.975)) | 
#                                        genedif_NO_C_EE_C$Gst<c(quantile(na.omit(genedif_NO_C_EE_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_NO_C_EE_C$Gst))
# 
# length(na.omit(genedif_NO_C_TU_H$Gst[genedif_NO_C_TU_H$Gst>c(quantile(na.omit(genedif_NO_C_TU_H$Gst),0.975)) | 
#                                        genedif_NO_C_TU_H$Gst<c(quantile(na.omit(genedif_NO_C_TU_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_NO_C_TU_H$Gst))
# 
# ## NO_H
# length(na.omit(genedif_NO_H_DK_C$Gst[genedif_NO_H_DK_C$Gst>c(quantile(na.omit(genedif_NO_H_DK_C$Gst),0.975)) | 
#                                        genedif_NO_H_DK_C$Gst<c(quantile(na.omit(genedif_NO_H_DK_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_NO_H_DK_C$Gst))
# 
# length(na.omit(genedif_NO_H_DK_H$Gst[genedif_NO_H_DK_H$Gst>c(quantile(na.omit(genedif_NO_H_DK_H$Gst),0.975)) | 
#                                        genedif_NO_H_DK_H$Gst<c(quantile(na.omit(genedif_NO_H_DK_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_NO_H_DK_H$Gst))
# 
# length(na.omit(genedif_NO_H_EE_C$Gst[genedif_NO_H_EE_C$Gst>c(quantile(na.omit(genedif_NO_H_EE_C$Gst),0.975)) | 
#                                        genedif_NO_H_EE_C$Gst<c(quantile(na.omit(genedif_NO_H_EE_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_NO_H_EE_C$Gst))
# 
# length(na.omit(genedif_NO_H_TU_H$Gst[genedif_NO_H_TU_H$Gst>c(quantile(na.omit(genedif_NO_H_TU_H$Gst),0.975)) | 
#                                        genedif_NO_H_TU_H$Gst<c(quantile(na.omit(genedif_NO_H_TU_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_NO_H_TU_H$Gst))
# 
# ## DK_C
# length(na.omit(genedif_DK_C_DK_H$Gst[genedif_DK_C_DK_H$Gst>c(quantile(na.omit(genedif_DK_C_DK_H$Gst),0.975)) | 
#                                        genedif_DK_C_DK_H$Gst<c(quantile(na.omit(genedif_DK_C_DK_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_DK_C_DK_H$Gst))
# 
# length(na.omit(genedif_DK_C_EE_C$Gst[genedif_DK_C_EE_C$Gst>c(quantile(na.omit(genedif_DK_C_EE_C$Gst),0.975)) | 
#                                        genedif_DK_C_EE_C$Gst<c(quantile(na.omit(genedif_DK_C_EE_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_DK_C_EE_C$Gst))
# 
# length(na.omit(genedif_DK_C_TU_H$Gst[genedif_DK_C_TU_H$Gst>c(quantile(na.omit(genedif_DK_C_TU_H$Gst),0.975)) | 
#                                        genedif_DK_C_TU_H$Gst<c(quantile(na.omit(genedif_DK_C_TU_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_DK_C_TU_H$Gst))
# 
# ## DK_H
# length(na.omit(genedif_DK_H_EE_C$Gst[genedif_DK_H_EE_C$Gst>c(quantile(na.omit(genedif_DK_H_EE_C$Gst),0.975)) | 
#                                        genedif_DK_H_EE_C$Gst<c(quantile(na.omit(genedif_DK_H_EE_C$Gst),0.025)) ]))/
#   length(na.omit(genedif_DK_H_EE_C$Gst))
# 
# length(na.omit(genedif_DK_H_TU_H$Gst[genedif_DK_H_TU_H$Gst>c(quantile(na.omit(genedif_DK_H_TU_H$Gst),0.975)) | 
#                                        genedif_DK_H_TU_H$Gst<c(quantile(na.omit(genedif_DK_H_TU_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_DK_H_TU_H$Gst))
# 
# ## EE_C
# length(na.omit(genedif_EE_C_TU_H$Gst[genedif_EE_C_TU_H$Gst>c(quantile(na.omit(genedif_EE_C_TU_H$Gst),0.975)) | 
#                                        genedif_EE_C_TU_H$Gst<c(quantile(na.omit(genedif_EE_C_TU_H$Gst),0.025)) ]))/
#   length(na.omit(genedif_EE_C_TU_H$Gst))
# 


## MDS distances
gst_mat_for_mds<-read.table("Gst_sd_matrix_for_MDS_270521.txt")
pops_for_mds<-read.table("pops_for_MDS_270521.txt")
gst_cmd<-cmdscale(as.dist(gst_mat_for_mds)) #t:turns the matrix around, as.dist:changes it to distance based, [-1,-1]:takes away the korean individual, cmdscale: changes it to a coordinates
gst_cmd
#pdf("gst_cmd_270521.pdf", width = 13, height = 8)
par(mar=c(5,5,1,1))
plot(gst_cmd, ylim=c(-0.075,0.1), xlim=c(-0.25, 0.25), xlab="MDS 1", ylab="MDS 2", pch=19, col="black", cex=1.5, cex.axis=2, cex.lab=2)
text(gst_cmd, pops_for_mds$V1, cex=2.5, font=3, pos=2) #puts in the pop names
#dev.off()


#### Theta and nucleotide diversity PopGenome
vcf_92_0.05_popgenome<-PopGenome::readData(path = "VCF_PopGenome", format="VCF")
vcf_92_0.05_popgenome<-set.populations(vcf_92_0.05_popgenome, list(1:25, 26:37, 38:42, 43:53, 
                                                                   54:56, 57:68, 69:76, 77:78,
                                                                   79:91, 92, 1:92))
names(vcf_92_0.05_popgenome@populations)<-c("IS_C", "NO_C", "DK_H", "DK_C", "EE_C", 
                                            "GL_C", "GL_H", "IS_H", "NO_H", "TU_H", "ALL")

vcf_92_0.05_popgenome<-F_ST.stats(vcf_92_0.05_popgenome)

get.F_ST(vcf_92_0.05_popgenome)

get.diversity(vcf_92_0.05_popgenome)[[11]] #ALL
get.diversity(vcf_92_0.05_popgenome)[[1]] #IS_C
get.diversity(vcf_92_0.05_popgenome)[[2]] #NO_C
get.diversity(vcf_92_0.05_popgenome)[[3]] #DK_H
get.diversity(vcf_92_0.05_popgenome)[[4]] #DK_C
get.diversity(vcf_92_0.05_popgenome)[[5]] #EE_C
get.diversity(vcf_92_0.05_popgenome)[[6]] #GL_C
get.diversity(vcf_92_0.05_popgenome)[[7]] #GL_H
get.diversity(vcf_92_0.05_popgenome)[[8]] #IS_H
get.diversity(vcf_92_0.05_popgenome)[[9]] #NO_H
get.diversity(vcf_92_0.05_popgenome)[[10]] #TU_H

win_vcf_92_0.05_popgenome <- sliding.window.transform(vcf_92_0.05_popgenome, 
                                    width=7, jump=5, 
                                    type=1,
                                    whole.data=T)

win_vcf_92_0.05_popgenome<- F_ST.stats(win_vcf_92_0.05_popgenome)
win_vcf_92_0.05_popgenome@nucleotide.F_ST
win_vcf_92_0.05_popgenome@nuc.diversity.within

get.neutrality(vcf_92_0.05_popgenome,theta=FALSE,stats=TRUE)

##### ROH
ROH_92_relaxed<- read.table("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_newchr26_relaxed_out_altered.hom.indiv" ,header=T)
ggplot(data=ROH_92_relaxed,
       aes(x=KB/1000, y=NSEG, color=CO_TI)) +
  geom_point()

ROH_92_relaxed_310321<-read.table("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_newchr26_relaxed_out_310321_1_altered.hom.indiv", header = T)
# ROH_92_relaxed_130921_plot<-
  ggplot(data=ROH_92_relaxed_310321,
       aes(x=KB/1000, y=NSEG, color=CO_TI,  size=CO_TI, fill=CO_TI)) + # for shape insert shape=CO_TI
  #scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9,10)) + # for shape uncomments 
  geom_point() +
  scale_size_manual(values=c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5), element_blank(), name="") +
  scale_color_manual(values = c("#1F78B4", "#A6CEE3", "#6A3D9A", "#33A02C", "#B2DF8A", "#E31A1C", "#FB9A99", "#FF7F00", "#FDBF6F", "#CAB2D6" ), name="Pop_Time") +
  scale_fill_manual(values = c("#1F78B4", "#A6CEE3", "#6A3D9A", "#33A02C", "#B2DF8A", "#E31A1C", "#FB9A99", "#FF7F00", "#FDBF6F", "#CAB2D6" ),  name="Pop_Time") +
  guides(size = F) + 
  guides(color = guide_legend(override.aes = list(size = 5))) + 
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
    text = element_text(size = 14))

# ggsave(plot = ROH_92_relaxed_130921_plot, filename = "ROHplot_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_newchr26_relaxed_out_310321_1_altered.hom.indiv_130921.png",
#        width = 30, height = 20, units = "cm", device = "png")

## ROH length vs missing data 
  
ROH_92_relaxed_310321_withmissing<-read.table("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_newchr26_relaxed_out_310321_1_altered.hom.indiv_withmissing", header = T)
 ROH_92_relaxed_130921_missingplot<-
ggplot(data=ROH_92_relaxed_310321_withmissing,
       aes(y=KB/1000, x=MISS, color=CO_TI,  size=CO_TI, fill=CO_TI)) + # for shape insert shape=CO_TI
  #scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9,10)) + # for shape uncomments 
  geom_point() +
  scale_size_manual(values=c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5), element_blank(), name="") +
  scale_color_manual(values = c("#1F78B4", "#A6CEE3", "#6A3D9A", "#33A02C", "#B2DF8A", "#E31A1C", "#FB9A99", "#FF7F00", "#FDBF6F", "#CAB2D6" ), name="Pop_Time") +
  scale_fill_manual(values = c("#1F78B4", "#A6CEE3", "#6A3D9A", "#33A02C", "#B2DF8A", "#E31A1C", "#FB9A99", "#FF7F00", "#FDBF6F", "#CAB2D6" ),  name="Pop_Time") +
  guides(size = F) + 
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  #ylim(0,500) +
  labs(y ="Length of ROH (Mb)", x = "Missingness") +
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
    text = element_text(size = 14))

# ggsave(plot = ROH_92_relaxed_130921_missingplot, filename = "ROHplot_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_newchr26_relaxed_out_310321_1_altered.hom.indiv_againstmissing_280621.png",
#        width = 30, height = 20, units = "cm", device = "png")

#### D-statistics
Dstat_92_bald_BBAA<-read.table("Dsuite/all_results_mac1_92ind_and_baldnewname_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_1_BBAA.txt", header=T)
plot(Dstat_92_bald_BBAA$Z.score[Dstat_92_bald_BBAA$p.value<0.05])
length(Dstat_92_bald_BBAA$Z.score)
length(Dstat_92_bald_BBAA$p.value[Dstat_92_bald_BBAA$p.value<0.05])
plot(Dstat_92_bald_BBAA$Dstatistic[Dstat_92_bald_BBAA$p.value<0.05])
Dstat_92_bald_Dmin<-read.table("Dsuite/all_results_mac1_92ind_and_baldnewname_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_1_Dmin.txt", header=T)
plot(Dstat_92_bald_Dmin$Z.score[Dstat_92_bald_Dmin$p.value<0.05])
plot(Dstat_92_bald_Dmin$Dstatistic[Dstat_92_bald_Dmin$p.value<0.05])

#### Treemix
plotting_funcs<-source("plotting_funcs.R")
get_f("Treemix/output/all_results_mac1_92ind_and_baldnewname_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_rootbald_boot_noss_global_k500_1_0")
plot_tree("Treemix/output/all_results_mac1_92ind_and_baldnewname_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_rootbald_boot_noss_global_se_k500_1_5")
plot_tree("Treemix/output/all_results_mac1_92ind_and_baldnewname_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_rootbald_boot_noss_global_se_k500_30_0")

#pdf("treemix_all_results_mac1_92ind_and_baldnewname_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_rootbald_boot_noss_global_se_k500_1_5_010621.pdf", height = 8, width = 12)
par(mar=c(4.2,0.5,0.5,0.5))
plot_tree("Treemix/output/all_results_mac1_92ind_and_baldnewname_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_rootbald_boot_noss_global_se_k500_1_5")
#dev.off()

all_f_explained<-read.table("Treemix/output/f_values_tree_rootbald_boot_k500_noss_global_se_rep50_m0to5.txt",header = T)
all_likely<-read.table("Treemix/output/llikely_values_tree_rootbald_boot_k500_noss_global_se_rep50_m0to5.txt", header = T)

#pdf("treemix_residuals_all_results_mac1_92ind_and_baldnewname_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_rootbald_boot_noss_global_se_k500_1_5_010621.pdf", height = 8, width = 12)
par(mar=c(7,7,0.5,0.5))
plot_resid("Treemix/output/all_results_mac1_92ind_and_baldnewname_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_rootbald_boot_noss_global_se_k500_1_5", "Treemix/output/poporder_paperfriendly.txt")
#dev.off()

#pdf("treemix_explain_like_all_results_mac1_92ind_and_baldnewname_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_rootbald_boot_noss_global_se_k500_1_5_010621.pdf", height = 8, width = 12)
par(mfrow=c(1,2))
par(mar=c(4.2,4,0.5,0.5))
boxplot(x = all_f_explained, names = c("0", "1", "2", "3", "4", "5"),  xlab="m edges", ylab="proportion of explained variance")
boxplot(x = all_likely, names = c("0", "1", "2", "3", "4", "5"),  xlab="m edges", ylab="likelihood")
#dev.off()

#### vcf2sfs
vcf_gt_92<-vcf2gt("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.vcf", "Individuallist_pops_92inds.txt")
samSize(vcf_gt_92)
vcf_gt_92_GL_C_95<-vcf2gt("VCF_pops/GL_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_GL_C_12.txt")
vcf_gt_92_GL_H_95<-vcf2gt("VCF_pops/GL_hist_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_GL_H_8.txt")
vcf_gt_92_IS_C_95<-vcf2gt("VCF_pops/IS_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_IS_C_25.txt")
vcf_gt_92_IS_H_95<-vcf2gt("VCF_pops/IS_hist_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_IS_H_2.txt")
vcf_gt_92_NO_C_12_95<-vcf2gt("VCF_pops/NO_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_NO_C_12.txt")
vcf_gt_92_NO_C_11_95<-vcf2gt("VCF_pops/NO_modern_from92_noBB38_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_NO_C_11.txt")
vcf_gt_92_NO_H_95<-vcf2gt("VCF_pops/NO_hist_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_NO_H_13.txt")
vcf_gt_92_DK_C_95<-vcf2gt("VCF_pops/DK_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_DK_C_11.txt")
vcf_gt_92_DK_H_95<-vcf2gt("VCF_pops/DK_hist_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_DK_H_5.txt")
vcf_gt_92_EE_C_95<-vcf2gt("VCF_pops/EE_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_EE_C_3.txt")
vcf_gt_92_TU_H_95<-vcf2gt("VCF_pops/TU_hist_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_TU_H_1.txt")

vcf_gt_92_GL_C_90<-vcf2gt("VCF_pops/GL_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.90_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_GL_C_12.txt")
vcf_gt_92_GL_H_90<-vcf2gt("VCF_pops/GL_hist_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.90_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_GL_H_8.txt")
vcf_gt_92_IS_C_90<-vcf2gt("VCF_pops/IS_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.90_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_IS_C_25.txt")
vcf_gt_92_IS_H_90<-vcf2gt("VCF_pops/IS_hist_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.90_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_IS_H_2.txt")
vcf_gt_92_NO_C_12_90<-vcf2gt("VCF_pops/NO_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.90_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_NO_C_12.txt")
vcf_gt_92_NO_C_11_90<-vcf2gt("VCF_pops/NO_modern_from92_noBB38_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.90_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_NO_C_11.txt")
vcf_gt_92_NO_H_90<-vcf2gt("VCF_pops/NO_hist_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.90_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_NO_H_13.txt")
vcf_gt_92_DK_C_90<-vcf2gt("VCF_pops/DK_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.90_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_DK_C_11.txt")
vcf_gt_92_DK_H_90<-vcf2gt("VCF_pops/DK_hist_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.90_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_DK_H_5.txt")
vcf_gt_92_EE_C_90<-vcf2gt("VCF_pops/EE_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.90_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_EE_C_3.txt")
vcf_gt_92_TU_H_90<-vcf2gt("VCF_pops/TU_hist_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.90_HetHom_minMQ30_biallelic.recode.vcf", "Individuallist_TU_H_1.txt")

vcf_sfs_92_GL_C_95_raw<-gt2sfs.raw(vcf_gt_92_GL_C_95, "GL_C")
vcf_sfs_92_GL_C_90_raw<-gt2sfs.raw(vcf_gt_92_GL_C_90, "GL_C")
plot(as.data.frame(vcf_sfs_92_GL_C_95_raw)[c(2:24),1],as.data.frame(vcf_sfs_92_GL_C_95_raw)[c(2:24),2])
plot(as.data.frame(vcf_sfs_92_GL_C_90_raw)[c(2:24),1],as.data.frame(vcf_sfs_92_GL_C_90_raw)[c(2:24),2])

vcf_sfs_92_IS_C_95_raw<-gt2sfs.raw(vcf_gt_92_IS_C_95, "IS_C")
vcf_sfs_92_IS_C_90_raw<-gt2sfs.raw(vcf_gt_92_IS_C_90, "IS_C")
plot(as.data.frame(vcf_sfs_92_IS_C_95_raw)[c(2:50),1],as.data.frame(vcf_sfs_92_IS_C_95_raw)[c(2:50),2])
plot(as.data.frame(vcf_sfs_92_IS_C_90_raw)[c(2:50),1],as.data.frame(vcf_sfs_92_IS_C_90_raw)[c(2:50),2])

vcf_sfs_92_NO_C_12_95_raw<-gt2sfs.raw(vcf_gt_92_NO_C_12_95, "NO_C")
vcf_sfs_92_NO_C_12_90_raw<-gt2sfs.raw(vcf_gt_92_NO_C_12_90, "NO_C")
plot(as.data.frame(vcf_sfs_92_NO_C_12_95_raw)[c(2:24),1],as.data.frame(vcf_sfs_92_NO_C_12_95_raw)[c(2:24),2])
plot(as.data.frame(vcf_sfs_92_NO_C_12_90_raw)[c(2:24),1],as.data.frame(vcf_sfs_92_NO_C_12_90_raw)[c(2:24),2])

vcf_sfs_92_NO_C_11_95_raw<-gt2sfs.raw(vcf_gt_92_NO_C_11_95, "NO_C")
vcf_sfs_92_NO_C_11_90_raw<-gt2sfs.raw(vcf_gt_92_NO_C_11_90, "NO_C")
plot(as.data.frame(vcf_sfs_92_NO_C_11_95_raw)[c(2:22),1],as.data.frame(vcf_sfs_92_NO_C_11_95_raw)[c(2:22),2])
plot(as.data.frame(vcf_sfs_92_NO_C_11_90_raw)[c(2:22),1],as.data.frame(vcf_sfs_92_NO_C_11_90_raw)[c(2:22),2])

vcf_sfs_92_DK_C_95_raw<-gt2sfs.raw(vcf_gt_92_DK_C_95, "DK_C")
vcf_sfs_92_DK_C_90_raw<-gt2sfs.raw(vcf_gt_92_DK_C_90, "DK_C")
plot(as.data.frame(vcf_sfs_92_DK_C_95_raw)[c(2:22),1],as.data.frame(vcf_sfs_92_DK_C_95_raw)[c(2:22),2])
plot(as.data.frame(vcf_sfs_92_DK_C_90_raw)[c(2:22),1],as.data.frame(vcf_sfs_92_DK_C_90_raw)[c(2:22),2])

vcf_gt_92_GL_C_95_reduced120K<-vcf2gt("VCF_pops/GL_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic_redused_120Ksnps.recode.vcf", "Individuallist_GL_C_12.txt")
vcf_gt_92_IS_C_95_reduced120K<-vcf2gt("VCF_pops/IS_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic_redused_120Ksnps.recode.vcf", "Individuallist_IS_C_25.txt")
vcf_gt_92_NO_C_11_95_reduced120K<-vcf2gt("VCF_pops/NO_modern_from92_noBB38_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic_redused_120Ksnps.recode.vcf", "Individuallist_NO_C_11.txt")
vcf_gt_92_DK_C_95_reduced120K<-vcf2gt("VCF_pops/DK_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.95_HetHom_minMQ30_biallelic_redused_120Ksnps.recode.vcf", "Individuallist_DK_C_11.txt")

vcf_sfs_92_GL_C_95_reduced120K_raw<-gt2sfs.raw(vcf_gt_92_GL_C_95_reduced120K, "GL_C")
vcf_sfs_92_IS_C_95_reduced120K_raw<-gt2sfs.raw(vcf_gt_92_IS_C_95_reduced120K, "IS_C")
vcf_sfs_92_NO_C_11_95_reduced120K_raw<-gt2sfs.raw(vcf_gt_92_NO_C_11_95_reduced120K, "NO_C")
vcf_sfs_92_DK_C_95_reduced120K_raw<-gt2sfs.raw(vcf_gt_92_DK_C_95_reduced120K, "DK_C")

plot(as.data.frame(vcf_sfs_92_GL_C_95_reduced120K_raw)[c(2:24),1],as.data.frame(vcf_sfs_92_GL_C_95_reduced120K_raw)[c(2:24),2])
plot(as.data.frame(vcf_sfs_92_IS_C_95_reduced120K_raw)[c(2:50),1],as.data.frame(vcf_sfs_92_IS_C_95_reduced120K_raw)[c(2:50),2])
plot(as.data.frame(vcf_sfs_92_NO_C_11_95_reduced120K_raw)[c(2:22),1],as.data.frame(vcf_sfs_92_NO_C_11_95_reduced120K_raw)[c(2:22),2])
plot(as.data.frame(vcf_sfs_92_DK_C_95_reduced120K_raw)[c(2:22),1],as.data.frame(vcf_sfs_92_DK_C_95_reduced120K_raw)[c(2:22),2])

# write.table(vcf_sfs_92_GL_C_95_reduced120K_raw, "sfs_out_R/vcf_sfs_92_GL_C_raw_95_reduced120K_250521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_IS_C_95_reduced120K_raw, "sfs_out_R/vcf_sfs_92_IS_C_raw_95_reduced120K_250521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_NO_C_11_95_reduced120K_raw, "sfs_out_R/vcf_sfs_92_NO_C_11_raw_95_reduced120K_250521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_DK_C_95_reduced120K_raw, "sfs_out_R/vcf_sfs_92_DK_C_raw_95_reduced120K_250521.txt", quote = F, row.names = F, sep = "\t")

# write.table(vcf_sfs_92_GL_C_95_raw, "sfs_out_R/vcf_sfs_92_GL_C_raw_95_070521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_IS_C_95_raw, "sfs_out_R/vcf_sfs_92_IS_C_raw_95_070521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_NO_C_11_95_raw, "sfs_out_R/vcf_sfs_92_NO_C_11_raw_95_070521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_DK_C_95_raw, "sfs_out_R/vcf_sfs_92_DK_C_raw_95_070521.txt", quote = F, row.names = F, sep = "\t")

# write.table(vcf_sfs_92_GL_C_imp, "sfs_out_R/vcf_sfs_92_GL_C_imp_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_GL_H_imp, "sfs_out_R/vcf_sfs_92_GL_H_imp_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_IS_C_imp, "sfs_out_R/vcf_sfs_92_IS_C_imp_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_IS_H_imp, "sfs_out_R/vcf_sfs_92_IS_H_imp_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_NO_C_imp, "sfs_out_R/vcf_sfs_92_NO_C_imp_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_NO_H_imp, "sfs_out_R/vcf_sfs_92_NO_H_imp_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_DK_C_imp, "sfs_out_R/vcf_sfs_92_DK_C_imp_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_DK_H_imp, "sfs_out_R/vcf_sfs_92_DK_H_imp_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_EE_C_imp, "sfs_out_R/vcf_sfs_92_EE_C_imp_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_GL_C_raw, "sfs_out_R/vcf_sfs_92_GL_C_raw_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_GL_H_raw, "sfs_out_R/vcf_sfs_92_GL_H_raw_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_IS_C_raw, "sfs_out_R/vcf_sfs_92_IS_C_raw_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_IS_H_raw, "sfs_out_R/vcf_sfs_92_IS_H_raw_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_NO_C_raw, "sfs_out_R/vcf_sfs_92_NO_C_raw_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_NO_H_raw, "sfs_out_R/vcf_sfs_92_NO_H_raw_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_DK_C_raw, "sfs_out_R/vcf_sfs_92_DK_C_raw_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_DK_H_raw, "sfs_out_R/vcf_sfs_92_DK_H_raw_050521.txt", quote = F, row.names = F, sep = "\t")
# write.table(vcf_sfs_92_EE_C_raw, "sfs_out_R/vcf_sfs_92_EE_C_raw_050521.txt", quote = F, row.names = F, sep = "\t")

##################################
##### Wattersons theta (S) First with pegas - virker ikke! 

#vcf_pegas<-pegas::read.vcf("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.vcf", from= 1, to = 800000)
vcf_92_vcfR<-read.vcfR("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.vcf")
vcf_92_dnabin<-vcfR2DNAbin(vcf_92_vcfR)
names(vcf_92_dnabin)

## No of segregating sites
seg.sites(vcf_92_dnabin)
length(seg.sites(vcf_92_dnabin[1:50,]))
length(seg.sites(mtDNA[sub.eagles_all=='Norway']))
length(seg.sites(mtDNA[sub.eagles_all=='Denmark']))
length(seg.sites(mtDNA[sub.eagles_all=='Estonia']))
length(seg.sites(mtDNA[sub.eagles_all=='Greenland']))

## theta S 
theta.s(124, 89) # All
theta.s(32, 42) # Iceland 
theta.s(30, 21) # Norway 
theta.s(72, 11) # Denmark 
theta.s(7, 3) # Estonia 
theta.s(25, 12) # Greenland 

## theta S sd 
sqrt(theta.s(124, 89, variance = T)[2])
sqrt(theta.s(32, 42, variance = T)[2]) # Iceland 
sqrt(theta.s(30, 21, variance = T)[2]) # Norway 
sqrt(theta.s(72, 11, variance = T)[2]) # Denmark 
sqrt(theta.s(7, 3, variance = T)[2]) # Estonia 
sqrt(theta.s(25, 12, variance = T)[2]) # Greenland 


# source("C:/Users/ccrha/Dropbox/PhD eagle Island/SambaR-master/SambaR-master/SAMBAR_v1.03.txt")
# getpackages(myrepos='http://cran.us.r-project.org',mylib=NULL)
# importdata(inputprefix="all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30",sumstatsfile=FALSE,depthfile=FALSE)
# filterdata(indmiss=0.25,snpmiss=0.1,min_mac=2,dohefilter=TRUE)	


#################################################
### Heterozygosity per site 

Hobs=function(x){
  x=x[!is.na(x)]
  length(x[x=="1"])/length(x)
}

Obs_het_GL_C<-Hobs(vcf_gt_92_GL_C_95$genotype[])


####################################################
### Stairway

GL_C_9<- read.table("Stairway_out/GL_C_raw_95_250521_T15.6_mut2.3_10-9.final.summary", head=T)
IS_C_9<- read.table("Stairway_out/IS_C_raw_95_250521_T15.6_mut2.3_10-9.final.summary", head=T)
NO_C_9<- read.table("Stairway_out/NO_C_raw_95_250521_T15.6_mut2.3_10-9.final.summary", head=T)
DK_C_9<- read.table("Stairway_out/DK_C_raw_95_250521_T15.6_mut2.3_10-9.final.summary", head=T)

#pdf("stairway_unfolded_12ind_T2.pdf")
####
plot(GL_C_9$year, GL_C_9$Ne_median, col="green", type="l", xlim=c(0,1000000), ylim=c(0,1000000))
lines(GL_C_9$year, GL_C_9$Ne_2.5., col="green", lty=2, lwd=0.5)
lines(GL_C_9$year, GL_C_9$Ne_97.5., col="green", lty=2, lwd=0.5)
plot(IS_C_9$year, IS_C_9$Ne_median, type="l", xlim=c(0,1000000), ylim=c(0,1000000))
lines(IS_C_9$year, IS_C_9$Ne_2.5., col="green", lty=2, lwd=0.5)
lines(IS_C_9$year, IS_C_9$Ne_97.5., col="green", lty=2, lwd=0.5)
plot(NO_C_9$year, NO_C_9$Ne_median, type="l", xlim=c(0,1000000), ylim=c(0,1000000))
lines(NO_C_9$year, NO_C_9$Ne_2.5., col="green", lty=2, lwd=0.5)
lines(NO_C_9$year, NO_C_9$Ne_97.5., col="green", lty=2, lwd=0.5)
plot(DK_C_9$year, DK_C_9$Ne_median, type="l", xlim=c(0,1000000), ylim=c(0,1000000))
lines(DK_C_9$year, DK_C_9$Ne_2.5., col="green", lty=2, lwd=0.5)
lines(DK_C_9$year, DK_C_9$Ne_97.5., col="green", lty=2, lwd=0.5)

### All in one
plot(GL_C_9$year, GL_C_9$Ne_median, col="green", type="l", xlim=c(0,1000000), ylim=c(0,1000000))
lines(GL_C_9$year, GL_C_9$Ne_2.5., col="green", lty=2, lwd=0.5)
lines(GL_C_9$year, GL_C_9$Ne_97.5., col="green", lty=2, lwd=0.5)

lines(IS_C_9$year, IS_C_9$Ne_median, col="red", lwd=2)
lines(IS_C_9$year, IS_C_9$Ne_2.5., col="red", lty=2, lwd=0.5)
lines(IS_C_9$year, IS_C_9$Ne_97.5., col="red", lty=2, lwd=0.5)

lines(NO_C_9$year, NO_C_9$Ne_median, col="orange", lwd=2)
lines(NO_C_9$year, NO_C_9$Ne_2.5., col="orange", lty=2, lwd=0.5)
lines(NO_C_9$year, NO_C_9$Ne_97.5., col="orange", lty=2, lwd=0.5)

lines(DK_C_9$year, DK_C_9$Ne_median, col="blue", lwd=2)
lines(DK_C_9$year, DK_C_9$Ne_2.5., col="blue", lty=2, lwd=0.5)
lines(DK_C_9$year, DK_C_9$Ne_97.5., col="blue", lty=2, lwd=0.5)

### All in one - 100K years
#pdf("plot_stairway_4pops_raw_95_250521_T15.6_mut2.3_10-9_100Kyears_010621.pdf", height = 8, width = 12)
par(mar=c(6,4.5,0.5,0.5))
plot(GL_C_9$year, GL_C_9$Ne_median/1000, col="green", type="l", xlim=c(0,100000), ylim=c(0,150), cex.axis=1.5, cex.lab=1.5, ylab="Ne/1000", xlab = "years back in time")
lines(GL_C_9$year, GL_C_9$Ne_2.5./1000, col="green", lty=2, lwd=1)
lines(GL_C_9$year, GL_C_9$Ne_97.5./1000, col="green", lty=2, lwd=1)

lines(IS_C_9$year, IS_C_9$Ne_median/1000, col="red", lwd=1)
lines(IS_C_9$year, IS_C_9$Ne_2.5./1000, col="red", lty=2, lwd=1)
lines(IS_C_9$year, IS_C_9$Ne_97.5./1000, col="red", lty=2, lwd=1)

lines(NO_C_9$year, NO_C_9$Ne_median/1000, col="orange", lwd=1)
lines(NO_C_9$year, NO_C_9$Ne_2.5./1000, col="orange", lty=2, lwd=1)
lines(NO_C_9$year, NO_C_9$Ne_97.5./1000, col="orange", lty=2, lwd=1)

lines(DK_C_9$year, DK_C_9$Ne_median/1000, col="blue", lwd=1)
lines(DK_C_9$year, DK_C_9$Ne_2.5./1000, col="blue", lty=2, lwd=1)
lines(DK_C_9$year, DK_C_9$Ne_97.5./1000, col="blue", lty=2, lwd=1)
legend("topleft",        # Add legend to plot
       legend = c("GL_C", "IS_C", "NO_C", "DK_C"),
       col = c("green", "red", "orange", "blue"),
       pch = 16,
       cex = 1.5)
#dev.off()

### All in one - 10K years
#pdf("plot_stairway_4pops_raw_95_250521_T15.6_mut2.3_10-9_10Kyears_010621.pdf", height = 8, width = 12)
par(mar=c(6,4.5,0.5,0.5))
plot(GL_C_9$year, GL_C_9$Ne_median/1000, col="green", type="l", xlim=c(0,10000), ylim=c(0,30), cex.axis=1.5, cex.lab=1.5, ylab="Ne/1000", xlab = "years back in time")
lines(GL_C_9$year, GL_C_9$Ne_2.5./1000, col="green", lty=2, lwd=1)
lines(GL_C_9$year, GL_C_9$Ne_97.5./1000, col="green", lty=2, lwd=1)

lines(IS_C_9$year, IS_C_9$Ne_median/1000, col="red", lwd=1)
lines(IS_C_9$year, IS_C_9$Ne_2.5./1000, col="red", lty=2, lwd=1)
lines(IS_C_9$year, IS_C_9$Ne_97.5./1000, col="red", lty=2, lwd=1)

lines(NO_C_9$year, NO_C_9$Ne_median/1000, col="orange", lwd=1)
lines(NO_C_9$year, NO_C_9$Ne_2.5./1000, col="orange", lty=2, lwd=1)
lines(NO_C_9$year, NO_C_9$Ne_97.5./1000, col="orange", lty=2, lwd=1)

lines(DK_C_9$year, DK_C_9$Ne_median/1000, col="blue", lwd=1)
lines(DK_C_9$year, DK_C_9$Ne_2.5./1000, col="blue", lty=2, lwd=1)
lines(DK_C_9$year, DK_C_9$Ne_97.5./1000, col="blue", lty=2, lwd=1)
legend("topleft",        # Add legend to plot
       legend = c("GL_C", "IS_C", "NO_C", "DK_C"),
       col = c("green", "red", "orange", "blue"),
       pch = 16,
       cex = 1.5)
#dev.off()

### All in one - 3K years
plot(GL_C_9$year, GL_C_9$Ne_median/1000, col="green", type="l", xlim=c(0,3000), ylim=c(0,20))
lines(GL_C_9$year, GL_C_9$Ne_2.5./1000, col="green", lty=2, lwd=0.5)
lines(GL_C_9$year, GL_C_9$Ne_97.5./1000, col="green", lty=2, lwd=0.5)

lines(IS_C_9$year, IS_C_9$Ne_median/1000, col="red", lwd=2)
lines(IS_C_9$year, IS_C_9$Ne_2.5./1000, col="red", lty=2, lwd=0.5)
lines(IS_C_9$year, IS_C_9$Ne_97.5./1000, col="red", lty=2, lwd=0.5)

lines(NO_C_9$year, NO_C_9$Ne_median/1000, col="orange", lwd=2)
lines(NO_C_9$year, NO_C_9$Ne_2.5./1000, col="orange", lty=2, lwd=0.5)
lines(NO_C_9$year, NO_C_9$Ne_97.5./1000, col="orange", lty=2, lwd=0.5)

lines(DK_C_9$year, DK_C_9$Ne_median/1000, col="blue", lwd=2)
lines(DK_C_9$year, DK_C_9$Ne_2.5./1000, col="blue", lty=2, lwd=0.5)
lines(DK_C_9$year, DK_C_9$Ne_97.5./1000, col="blue", lty=2, lwd=0.5)

##### 10-8 all sites
GL_C_8_all<- read.table("Stairway_out/GL_C_raw_95_180521_T15.6_mut10-8.final.summary", head=T)
IS_C_8_all<- read.table("Stairway_out/IS_C_raw_95_180521_T15.6_mut10-8.final.summary", head=T)
NO_C_8_all<- read.table("Stairway_out/NO_C_11_raw_95_180521_T15.6_mut10-8.final.summary", head=T)
DK_C_8_all<- read.table("Stairway_out/DK_C_raw_95_180521_T15.6_mut10-8.final.summary", head=T)

plot(GL_C_8_all$year, GL_C_8_all$Ne_median, col="green", type="l", xlim=c(0,100000), ylim=c(0,10000))
lines(IS_C_8_all$year, IS_C_8_all$Ne_median, col="red")#, type="l",  xlim=c(0,100000), ylim=c(0,10000))
lines(NO_C_8_all$year, NO_C_8_all$Ne_median, col="orange")#, type="l",  xlim=c(0,100000), ylim=c(0,10000))
lines(DK_C_8_all$year, DK_C_8_all$Ne_median, col="blue")#, type="l", col="blue")#, xlim=c(0,100000), ylim=c(0,10000))

#### 5*10-9 

GL_C_5_109<- read.table("Stairway_out/GL_C_raw_95_290621_120Ksnps_T15.6_mut5_10-9.final.summary", head=T)
IS_C_5_109<- read.table("Stairway_out/IS_C_raw_95_290621_120Ksnps_T15.6_mut5_10-9.final.summary", head=T)
NO_C_5_109<- read.table("Stairway_out/NO_C_11_raw_95_290621_120Ksnps_T15.6_mut5_10-9.final.summary", head=T)
DK_C_5_109<- read.table("Stairway_out/DK_C_raw_95_290621_T15.6_mut5_10-9.final.summary", head=T)

### All in one 5*10-9 - 200K years
#pdf("plot_stairway_4pops_raw_95_290621_T15.6_mut5_10-9_200Kyears_290621.pdf", height = 8, width = 12)
par(mar=c(6,4.5,0.5,0.5))
plot(GL_C_5_109$year, GL_C_5_109$Ne_median/1000, col="green", type="l", xlim=c(0,200000), ylim=c(0,200), cex.axis=1.5, cex.lab=1.5, ylab="Ne/1000", xlab = "years back in time")
lines(GL_C_5_109$year, GL_C_5_109$Ne_2.5./1000, col="green", lty=2, lwd=1)
lines(GL_C_5_109$year, GL_C_5_109$Ne_97.5./1000, col="green", lty=2, lwd=1)

lines(IS_C_5_109$year, IS_C_5_109$Ne_median/1000, col="red", lwd=1)
lines(IS_C_5_109$year, IS_C_5_109$Ne_2.5./1000, col="red", lty=2, lwd=1)
lines(IS_C_5_109$year, IS_C_5_109$Ne_97.5./1000, col="red", lty=2, lwd=1)

lines(NO_C_5_109$year, NO_C_5_109$Ne_median/1000, col="orange", lwd=1)
lines(NO_C_5_109$year, NO_C_5_109$Ne_2.5./1000, col="orange", lty=2, lwd=1)
lines(NO_C_5_109$year, NO_C_5_109$Ne_97.5./1000, col="orange", lty=2, lwd=1)

lines(DK_C_5_109$year, DK_C_5_109$Ne_median/1000, col="blue", lwd=1)
lines(DK_C_5_109$year, DK_C_5_109$Ne_2.5./1000, col="blue", lty=2, lwd=1)
lines(DK_C_5_109$year, DK_C_5_109$Ne_97.5./1000, col="blue", lty=2, lwd=1)
legend("topleft",        # Add legend to plot
       legend = c("GL_C", "IS_C", "NO_C", "DK_C"),
       col = c("green", "red", "orange", "blue"),
       pch = 16,
       cex = 1.5)
#dev.off()

### All in one 1*10-9 - 10K years
#pdf("plot_stairway_4pops_raw_95_290621_T15.6_mut5_10-9_10Kyears_290621.pdf", height = 8, width = 12)
par(mar=c(6,4.5,0.5,0.5))
plot(GL_C_5_109$year, GL_C_5_109$Ne_median/1000, col="green", type="l", xlim=c(0,10000), ylim=c(0,30), cex.axis=1.5, cex.lab=1.5, ylab="Ne/1000", xlab = "years back in time")
lines(GL_C_5_109$year, GL_C_5_109$Ne_2.5./1000, col="green", lty=2, lwd=1)
lines(GL_C_5_109$year, GL_C_5_109$Ne_97.5./1000, col="green", lty=2, lwd=1)

lines(IS_C_5_109$year, IS_C_5_109$Ne_median/1000, col="red", lwd=1)
lines(IS_C_5_109$year, IS_C_5_109$Ne_2.5./1000, col="red", lty=2, lwd=1)
lines(IS_C_5_109$year, IS_C_5_109$Ne_97.5./1000, col="red", lty=2, lwd=1)

lines(NO_C_5_109$year, NO_C_5_109$Ne_median/1000, col="orange", lwd=1)
lines(NO_C_5_109$year, NO_C_5_109$Ne_2.5./1000, col="orange", lty=2, lwd=1)
lines(NO_C_5_109$year, NO_C_5_109$Ne_97.5./1000, col="orange", lty=2, lwd=1)

lines(DK_C_5_109$year, DK_C_5_109$Ne_median/1000, col="blue", lwd=1)
lines(DK_C_5_109$year, DK_C_5_109$Ne_2.5./1000, col="blue", lty=2, lwd=1)
lines(DK_C_5_109$year, DK_C_5_109$Ne_97.5./1000, col="blue", lty=2, lwd=1)
legend("topleft",        # Add legend to plot
       legend = c("GL_C", "IS_C", "NO_C", "DK_C"),
       col = c("green", "red", "orange", "blue"),
       pch = 16,
       cex = 1.5)
#dev.off()


#### 1*10-9 

GL_C_1_109<- read.table("Stairway_out/GL_C_raw_95_250521_T15.6_mut10-9.final.summary", head=T)
IS_C_1_109<- read.table("Stairway_out/IS_C_raw_95_250521_T15.6_mut10-9.final.summary", head=T)
NO_C_1_109<- read.table("Stairway_out/NO_C_raw_95_250521_T15.6_mut10-9.final.summary", head=T)
DK_C_1_109<- read.table("Stairway_out/DK_C_raw_95_250521_T15.6_mut10-9.final.summary", head=T)

### All in one 1*10-9 - 200K years
#pdf("plot_stairway_4pops_raw_95_290621_T15.6_mut5_10-9_200Kyears_290621.pdf", height = 8, width = 12)
par(mar=c(6,4.5,0.5,0.5))
plot(GL_C_1_109$year, GL_C_1_109$Ne_median/1000, col="green", type="l", xlim=c(0,100000), ylim=c(0,200), cex.axis=1.5, cex.lab=1.5, ylab="Ne/1000", xlab = "years back in time")
lines(GL_C_1_109$year, GL_C_1_109$Ne_2.5./1000, col="green", lty=2, lwd=1)
lines(GL_C_1_109$year, GL_C_1_109$Ne_97.5./1000, col="green", lty=2, lwd=1)

lines(IS_C_1_109$year, IS_C_1_109$Ne_median/1000, col="red", lwd=1)
lines(IS_C_1_109$year, IS_C_1_109$Ne_2.5./1000, col="red", lty=2, lwd=1)
lines(IS_C_1_109$year, IS_C_1_109$Ne_97.5./1000, col="red", lty=2, lwd=1)

lines(NO_C_1_109$year, NO_C_1_109$Ne_median/1000, col="orange", lwd=1)
lines(NO_C_1_109$year, NO_C_1_109$Ne_2.5./1000, col="orange", lty=2, lwd=1)
lines(NO_C_1_109$year, NO_C_1_109$Ne_97.5./1000, col="orange", lty=2, lwd=1)

lines(DK_C_1_109$year, DK_C_1_109$Ne_median/1000, col="blue", lwd=1)
lines(DK_C_1_109$year, DK_C_1_109$Ne_2.5./1000, col="blue", lty=2, lwd=1)
lines(DK_C_1_109$year, DK_C_1_109$Ne_97.5./1000, col="blue", lty=2, lwd=1)
legend("topleft",        # Add legend to plot
       legend = c("GL_C", "IS_C", "NO_C", "DK_C"),
       col = c("green", "red", "orange", "blue"),
       pch = 16,
       cex = 1.5)
#dev.off()

### All in one 1*10-9 - 10K years
#pdf("plot_stairway_4pops_raw_95_250521_T15.6_mut1_10-9_10Kyears_290621.pdf", height = 8, width = 12)
par(mar=c(6,4.5,0.5,0.5))
plot(GL_C_1_109$year, GL_C_1_109$Ne_median/1000, col="green", type="l", xlim=c(0,10000), ylim=c(0,40), cex.axis=1.5, cex.lab=1.5, ylab="Ne/1000", xlab = "years back in time")
lines(GL_C_1_109$year, GL_C_1_109$Ne_2.5./1000, col="green", lty=2, lwd=1)
lines(GL_C_1_109$year, GL_C_1_109$Ne_97.5./1000, col="green", lty=2, lwd=1)

lines(IS_C_1_109$year, IS_C_1_109$Ne_median/1000, col="red", lwd=1)
lines(IS_C_1_109$year, IS_C_1_109$Ne_2.5./1000, col="red", lty=2, lwd=1)
lines(IS_C_1_109$year, IS_C_1_109$Ne_97.5./1000, col="red", lty=2, lwd=1)

lines(NO_C_1_109$year, NO_C_1_109$Ne_median/1000, col="orange", lwd=1)
lines(NO_C_1_109$year, NO_C_1_109$Ne_2.5./1000, col="orange", lty=2, lwd=1)
lines(NO_C_1_109$year, NO_C_1_109$Ne_97.5./1000, col="orange", lty=2, lwd=1)

lines(DK_C_1_109$year, DK_C_1_109$Ne_median/1000, col="blue", lwd=1)
lines(DK_C_1_109$year, DK_C_1_109$Ne_2.5./1000, col="blue", lty=2, lwd=1)
lines(DK_C_1_109$year, DK_C_1_109$Ne_97.5./1000, col="blue", lty=2, lwd=1)
legend("topleft",        # Add legend to plot
       legend = c("GL_C", "IS_C", "NO_C", "DK_C"),
       col = c("green", "red", "orange", "blue"),
       pch = 16,
       cex = 1.5)
#dev.off()



#### 1*10-9 ANGSD out 100000000

GL_C_1_109_ang<- read.table("Stairway_out/GL_C_angsd_080721_10boot_100miosites_T15.6_mut10-9.final.summary", head=T)
IS_C_1_109_ang<- read.table("Stairway_out/IS_C_angsd_080721_10boot_100miosites_T15.6_mut10-9.final.summary", head=T)
NO_C_1_109_ang<- read.table("Stairway_out/NO_C_11_angsd_080721_10boot_100miosites_T15.6_mut10-9.final.summary", head=T)
DK_C_1_109_ang<- read.table("Stairway_out/DK_C_angsd_080721_10boot_100miosites_T15.6_mut10-9.final.summary", head=T)

### All in one 1*10-9 - 100K years
#pdf("plot_stairway_4pops_raw_95_290621_T15.6_mut5_10-9_200Kyears_290621.pdf", height = 8, width = 12)
par(mar=c(6,4.5,0.5,0.5))
plot(GL_C_1_109_ang$year, GL_C_1_109_ang$Ne_median/1000, col="green", type="l", xlim=c(0,100000), ylim=c(0,200), cex.axis=1.5, cex.lab=1.5, ylab="Ne/1000", xlab = "years back in time")
lines(GL_C_1_109_ang$year, GL_C_1_109_ang$Ne_2.5./1000, col="green", lty=2, lwd=1)
lines(GL_C_1_109_ang$year, GL_C_1_109_ang$Ne_97.5./1000, col="green", lty=2, lwd=1)

lines(IS_C_1_109_ang$year, IS_C_1_109_ang$Ne_median/1000, col="red", lwd=1)
lines(IS_C_1_109_ang$year, IS_C_1_109_ang$Ne_2.5./1000, col="red", lty=2, lwd=1)
lines(IS_C_1_109_ang$year, IS_C_1_109_ang$Ne_97.5./1000, col="red", lty=2, lwd=1)

lines(NO_C_1_109_ang$year, NO_C_1_109_ang$Ne_median/1000, col="orange", lwd=1)
lines(NO_C_1_109_ang$year, NO_C_1_109_ang$Ne_2.5./1000, col="orange", lty=2, lwd=1)
lines(NO_C_1_109_ang$year, NO_C_1_109_ang$Ne_97.5./1000, col="orange", lty=2, lwd=1)

lines(DK_C_1_109_ang$year, DK_C_1_109_ang$Ne_median/1000, col="blue", lwd=1)
lines(DK_C_1_109_ang$year, DK_C_1_109_ang$Ne_2.5./1000, col="blue", lty=2, lwd=1)
lines(DK_C_1_109_ang$year, DK_C_1_109_ang$Ne_97.5./1000, col="blue", lty=2, lwd=1)
legend("topleft",        # Add legend to plot
       legend = c("GL_C", "IS_C", "NO_C", "DK_C"),
       col = c("green", "red", "orange", "blue"),
       pch = 16,
       cex = 1.5)
#dev.off()

### All in one 1*10-9 - 10K years
#pdf("plot_stairway_4pops_raw_95_250521_T15.6_mut1_10-9_10Kyears_290621.pdf", height = 8, width = 12)
par(mar=c(6,4.5,0.5,0.5))
plot(GL_C_1_109_ang$year, GL_C_1_109_ang$Ne_median/1000, col="green", type="l", xlim=c(0,50000), ylim=c(0,100), cex.axis=1.5, cex.lab=1.5, ylab="Ne/1000", xlab = "years back in time", options(scipen=5))
lines(GL_C_1_109_ang$year, GL_C_1_109_ang$Ne_2.5./1000, col="green", lty=2, lwd=1)
lines(GL_C_1_109_ang$year, GL_C_1_109_ang$Ne_97.5./1000, col="green", lty=2, lwd=1)

lines(IS_C_1_109_ang$year, IS_C_1_109_ang$Ne_median/1000, col="red", lwd=1)
lines(IS_C_1_109_ang$year, IS_C_1_109_ang$Ne_2.5./1000, col="red", lty=2, lwd=1)
lines(IS_C_1_109_ang$year, IS_C_1_109_ang$Ne_97.5./1000, col="red", lty=2, lwd=1)

lines(NO_C_1_109_ang$year, NO_C_1_109_ang$Ne_median/1000, col="orange", lwd=1)
lines(NO_C_1_109_ang$year, NO_C_1_109_ang$Ne_2.5./1000, col="orange", lty=2, lwd=1)
lines(NO_C_1_109_ang$year, NO_C_1_109_ang$Ne_97.5./1000, col="orange", lty=2, lwd=1)

lines(DK_C_1_109_ang$year, DK_C_1_109_ang$Ne_median/1000, col="blue", lwd=1)
lines(DK_C_1_109_ang$year, DK_C_1_109_ang$Ne_2.5./1000, col="blue", lty=2, lwd=1)
lines(DK_C_1_109_ang$year, DK_C_1_109_ang$Ne_97.5./1000, col="blue", lty=2, lwd=1)
legend("topleft",        # Add legend to plot
       legend = c("GL_C", "IS_C", "NO_C", "DK_C"),
       col = c("green", "red", "orange", "blue"),
       pch = 16,
       cex = 1.5)
#dev.off()


##### stairway ANGSD only auto par(mar=c(6,4.5,0.5,0.5))
### All in one 2.3*10-9 - 100K years
GL_C_2.3_109_ang_auto<- read.table("Stairway_out/angsd_sites_auto_120721/GL_C_angsd_120721_10boot_100miosites_T15.6_mut2.3_10-9.final.summary", head=T)
IS_C_2.3_109_ang_auto<- read.table("Stairway_out/angsd_sites_auto_120721/IS_C_angsd_120721_10boot_100miosites_T15.6_mut2.3_10-9.final.summary", head=T)
NO_C_2.3_109_ang_auto<- read.table("Stairway_out/angsd_sites_auto_120721/NO_C_11_angsd_120721_10boot_100miosites_T15.6_mut2.3_10-9.final.summary", head=T)
DK_C_2.3_109_ang_auto<- read.table("Stairway_out/angsd_sites_auto_120721/DK_C_angsd_120721_10boot_100miosites_T15.6_mut2.3_10-9.final.summary", head=T)

options(scipen=5)

#pdf("plot_stairway_4pops_raw_95_150721_T15.6_mut2.3_10-9_100Kyears_angsd_auto.pdf", height = 8, width = 12)
par(mar=c(6,4.5,0.5,0.5))
plot(GL_C_2.3_109_ang_auto$year, GL_C_2.3_109_ang_auto$Ne_median/1000, col="green", type="l", xlim=c(0,100000), ylim=c(0,200), cex.axis=1.5, cex.lab=1.5, ylab="Ne/1000", xlab = "years back in time")
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
legend("topleft",        # Add legend to plot
       legend = c("GL_C", "IS_C", "NO_C", "DK_C"),
       col = c("green", "red", "orange", "blue"),
       pch = 16,
       cex = 1.5)
#dev.off()



### All in one 2.3*10-9 - 10K years
#pdf("plot_stairway_4pops_raw_95_150721_T15.6_mut2.3_10-9_10Kyears_angsd_auto.pdf", height = 8, width = 12)
par(mar=c(6,4.5,0.5,0.5))
plot(GL_C_2.3_109_ang_auto$year, GL_C_2.3_109_ang_auto$Ne_median/1000, col="green", type="l", xlim=c(0,10000), ylim=c(0,50), cex.axis=1.5, cex.lab=1.5, ylab="Ne/1000", xlab = "years back in time")
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
legend("topleft",        # Add legend to plot
       legend = c("GL_C", "IS_C", "NO_C", "DK_C"),
       col = c("green", "red", "orange", "blue"),
       pch = 16,
       cex = 1.5)
#dev.off()

### All in one 1*10-9 - 100K years
GL_C_1_109_ang_auto<- read.table("Stairway_out/angsd_sites_auto_120721/GL_C_angsd_120721_10boot_100miosites_T15.6_mut10-9.final.summary", head=T)
IS_C_1_109_ang_auto<- read.table("Stairway_out/angsd_sites_auto_120721/IS_C_angsd_120721_10boot_100miosites_T15.6_mut10-9.final.summary", head=T)
NO_C_1_109_ang_auto<- read.table("Stairway_out/angsd_sites_auto_120721/NO_C_11_angsd_120721_10boot_100miosites_T15.6_mut10-9.final.summary", head=T)
DK_C_1_109_ang_auto<- read.table("Stairway_out/angsd_sites_auto_120721/DK_C_angsd_120721_10boot_100miosites_T15.6_mut10-9.final.summary", head=T)

#pdf("plot_stairway_4pops_raw_95_140721_T15.6_mut1_10-9_1000Kyears_angsd_auto.pdf", height = 8, width = 12)
par(mar=c(6,4.5,0.5,0.5))
plot(GL_C_1_109_ang_auto$year, GL_C_1_109_ang_auto$Ne_median/1000, col="green", type="l", xlim=c(0,1000000), ylim=c(0,2500), cex.axis=1.5, cex.lab=1.5, ylab="Ne/1000", xlab = "years back in time")
lines(GL_C_1_109_ang_auto$year, GL_C_1_109_ang_auto$Ne_2.5./1000, col="green", lty=2, lwd=1)
lines(GL_C_1_109_ang_auto$year, GL_C_1_109_ang_auto$Ne_97.5./1000, col="green", lty=2, lwd=1)

lines(IS_C_1_109_ang_auto$year, IS_C_1_109_ang_auto$Ne_median/1000, col="red", lwd=1)
lines(IS_C_1_109_ang_auto$year, IS_C_1_109_ang_auto$Ne_2.5./1000, col="red", lty=2, lwd=1)
lines(IS_C_1_109_ang_auto$year, IS_C_1_109_ang_auto$Ne_97.5./1000, col="red", lty=2, lwd=1)

lines(NO_C_1_109_ang_auto$year, NO_C_1_109_ang_auto$Ne_median/1000, col="orange", lwd=1)
lines(NO_C_1_109_ang_auto$year, NO_C_1_109_ang_auto$Ne_2.5./1000, col="orange", lty=2, lwd=1)
lines(NO_C_1_109_ang_auto$year, NO_C_1_109_ang_auto$Ne_97.5./1000, col="orange", lty=2, lwd=1)

lines(DK_C_1_109_ang_auto$year, DK_C_1_109_ang_auto$Ne_median/1000, col="blue", lwd=1)
lines(DK_C_1_109_ang_auto$year, DK_C_1_109_ang_auto$Ne_2.5./1000, col="blue", lty=2, lwd=1)
lines(DK_C_1_109_ang_auto$year, DK_C_1_109_ang_auto$Ne_97.5./1000, col="blue", lty=2, lwd=1)
legend("topleft",        # Add legend to plot
       legend = c("GL_C", "IS_C", "NO_C", "DK_C"),
       col = c("green", "red", "orange", "blue"),
       pch = 16,
       cex = 1.5)
#dev.off()

### All in one 1*10-9 - 10K years
#pdf("plot_stairway_4pops_raw_95_140721_T15.6_mut1_10-9_30Kyears_angsd_auto.pdf", height = 8, width = 12)
par(mar=c(6,4.5,0.5,0.5))
plot(GL_C_1_109_ang_auto$year, GL_C_1_109_ang_auto$Ne_median/1000, col="green", type="l", xlim=c(0,30000), ylim=c(0,100), cex.axis=1.5, cex.lab=1.5, ylab="Ne/1000", xlab = "years back in time")
lines(GL_C_1_109_ang_auto$year, GL_C_1_109_ang_auto$Ne_2.5./1000, col="green", lty=2, lwd=1)
lines(GL_C_1_109_ang_auto$year, GL_C_1_109_ang_auto$Ne_97.5./1000, col="green", lty=2, lwd=1)

lines(IS_C_1_109_ang_auto$year, IS_C_1_109_ang_auto$Ne_median/1000, col="red", lwd=1)
lines(IS_C_1_109_ang_auto$year, IS_C_1_109_ang_auto$Ne_2.5./1000, col="red", lty=2, lwd=1)
lines(IS_C_1_109_ang_auto$year, IS_C_1_109_ang_auto$Ne_97.5./1000, col="red", lty=2, lwd=1)

lines(NO_C_1_109_ang_auto$year, NO_C_1_109_ang_auto$Ne_median/1000, col="orange", lwd=1)
lines(NO_C_1_109_ang_auto$year, NO_C_1_109_ang_auto$Ne_2.5./1000, col="orange", lty=2, lwd=1)
lines(NO_C_1_109_ang_auto$year, NO_C_1_109_ang_auto$Ne_97.5./1000, col="orange", lty=2, lwd=1)

lines(DK_C_1_109_ang_auto$year, DK_C_1_109_ang_auto$Ne_median/1000, col="blue", lwd=1)
lines(DK_C_1_109_ang_auto$year, DK_C_1_109_ang_auto$Ne_2.5./1000, col="blue", lty=2, lwd=1)
lines(DK_C_1_109_ang_auto$year, DK_C_1_109_ang_auto$Ne_97.5./1000, col="blue", lty=2, lwd=1)       
legend("topleft",        # Add legend to plot
       legend = c("GL_C", "IS_C", "NO_C", "DK_C"),
       col = c("green", "red", "orange", "blue"),
       pch = 16,
       cex = 1.5)
#dev.off()



#####################
### Gst along the genome
### 
#fst<-read.table("genedif_GL_C_GL_H", header=TRUE)
genedif_GL_C_GL_H<-read.table("Gst_vcfR/Gst_genedif_GL_C_GL_H.txt", header=TRUE)
Gstsubset<-genedif_GL_C_GL_H[complete.cases(genedif_GL_C_GL_H),]
mydf_meangst<-tapply(mydf$Gst, list(mydf$CHROM, mydf$bin), mean)
SNP<-c(1:(nrow(Gstsubset)))
mydf<-data.frame(SNP,Gstsubset)
mydf$CHROM <- revalue(mydf$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
mydf$CHROM <- as.numeric(mydf$CHROM)
mydf$POS <- as.numeric(mydf$POS)
mydf$bin <- mydf$POS%/%100000
manhattan(mydf, chr="CHROM", bp="POS", p="Gst", snp="SNP", logp=FALSE, ylab="Nei´s Gst")

# Gst along the genome GL_C_GL_H
genedif_GL_C_GL_H_sub<-genedif_GL_C_GL_H[complete.cases(genedif_GL_C_GL_H),]
genedif_GL_C_GL_H_sub$bin<-genedif_GL_C_GL_H_sub$bin <- as.numeric(genedif_GL_C_GL_H_sub$POS)%/%100000
genedif_GL_C_GL_H_sub_meangst<-tapply(genedif_GL_C_GL_H_sub$Gst, list(genedif_GL_C_GL_H_sub$CHROM, genedif_GL_C_GL_H_sub$bin), mean)
new.genedif_GL_C_GL_H_sub_meangst=as.data.frame(cbind(rep(row.names(genedif_GL_C_GL_H_sub_meangst),each=855),rep(0:854,each=nrow(genedif_GL_C_GL_H_sub_meangst)),as.numeric(t(genedif_GL_C_GL_H_sub_meangst))))
GL_C_GL_H_sub_meangst<-rep(0:854,length.out=22230)
table_GL_C_GL_H_sub_meangs<-new.genedif_GL_C_GL_H_sub_meangst
table_GL_C_GL_H_sub_meangs$V2<-bins
SNP_table_GL_C_GL_H_sub_meangs<-c(1:(nrow(table_GL_C_GL_H_sub_meangs)))
df_table_GL_C_GL_H_sub_meangs<-data.frame(SNP_table_GL_C_GL_H_sub_meangs,table_GL_C_GL_H_sub_meangs)
colnames(df_table_GL_C_GL_H_sub_meangs) <- c("SNP","CHROM", "bin","Gst_mean")
df_table_GL_C_GL_H_sub_meangs$CHROM <- revalue(df_table_GL_C_GL_H_sub_meangs$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
df_table_GL_C_GL_H_sub_meangs$CHROM <- as.numeric(df_table_GL_C_GL_H_sub_meangs$CHROM)
df_table_GL_C_GL_H_sub_meangs$Gst_mean[df_table_GL_C_GL_H_sub_meangs$Gst_mean == NA] <- NA
df_table_GL_C_GL_H_sub_meangs_sub<-df_table_GL_C_GL_H_sub_meangs[complete.cases(df_table_GL_C_GL_H_sub_meangs),]
df_table_GL_C_GL_H_sub_meangs_sub$Gst_mean <- as.numeric(df_table_GL_C_GL_H_sub_meangs_sub$Gst_mean)
manhattan(df_table_GL_C_GL_H_sub_meangs_sub, chr="CHROM", bp="bin", p="Gst_mean", snp="SNP", logp=FALSE, ylab="Nei´s Gst")

# Gst along the genome GL_C_IS_C
genedif_GL_C_IS_C<-read.table("Gst_vcfR/Gst_genedif_GL_C_IS_C.txt", header=TRUE)
genedif_GL_C_IS_C_sub<-genedif_GL_C_IS_C[complete.cases(genedif_GL_C_IS_C),]
genedif_GL_C_IS_C_sub$bin<-genedif_GL_C_IS_C_sub$bin <- as.numeric(genedif_GL_C_IS_C_sub$POS)%/%100000
genedif_GL_C_IS_C_sub_meangst<-tapply(genedif_GL_C_IS_C_sub$Gst, list(genedif_GL_C_IS_C_sub$CHROM, genedif_GL_C_IS_C_sub$bin), mean)
new.genedif_GL_C_IS_C_sub_meangst=as.data.frame(cbind(rep(row.names(genedif_GL_C_IS_C_sub_meangst),each=855),rep(0:854,each=nrow(genedif_GL_C_IS_C_sub_meangst)),as.numeric(t(genedif_GL_C_IS_C_sub_meangst))))
GL_C_IS_C_sub_meangst<-rep(0:854,length.out=22230)
table_GL_C_IS_C_sub_meangs<-new.genedif_GL_C_IS_C_sub_meangst
table_GL_C_IS_C_sub_meangs$V2<-bins
SNP_table_GL_C_IS_C_sub_meangs<-c(1:(nrow(table_GL_C_IS_C_sub_meangs)))
df_table_GL_C_IS_C_sub_meangs<-data.frame(SNP_table_GL_C_IS_C_sub_meangs,table_GL_C_IS_C_sub_meangs)
colnames(df_table_GL_C_IS_C_sub_meangs) <- c("SNP","CHROM", "bin","Gst_mean")
df_table_GL_C_IS_C_sub_meangs$CHROM <- revalue(df_table_GL_C_IS_C_sub_meangs$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
df_table_GL_C_IS_C_sub_meangs$CHROM <- as.numeric(df_table_GL_C_IS_C_sub_meangs$CHROM)
df_table_GL_C_IS_C_sub_meangs$Gst_mean[df_table_GL_C_IS_C_sub_meangs$Gst_mean == NA] <- NA
df_table_GL_C_IS_C_sub_meangs_sub<-df_table_GL_C_IS_C_sub_meangs[complete.cases(df_table_GL_C_IS_C_sub_meangs),]
df_table_GL_C_IS_C_sub_meangs_sub$Gst_mean <- as.numeric(df_table_GL_C_IS_C_sub_meangs_sub$Gst_mean)
manhattan(df_table_GL_C_IS_C_sub_meangs_sub, chr="CHROM", bp="bin", p="Gst_mean", snp="SNP", logp=FALSE, ylab="Nei´s Gst")

# Gst along the genome GL_C_NO_C
genedif_GL_C_NO_C<-read.table("Gst_vcfR/Gst_genedif_GL_C_NO_C.txt", header=TRUE)
genedif_GL_C_NO_C_sub<-genedif_GL_C_NO_C[complete.cases(genedif_GL_C_NO_C),]
genedif_GL_C_NO_C_sub$bin<-genedif_GL_C_NO_C_sub$bin <- as.numeric(genedif_GL_C_NO_C_sub$POS)%/%100000
genedif_GL_C_NO_C_sub_meangst<-tapply(genedif_GL_C_NO_C_sub$Gst, list(genedif_GL_C_NO_C_sub$CHROM, genedif_GL_C_NO_C_sub$bin), mean)
new.genedif_GL_C_NO_C_sub_meangst=as.data.frame(cbind(rep(row.names(genedif_GL_C_NO_C_sub_meangst),each=855),rep(0:854,each=nrow(genedif_GL_C_NO_C_sub_meangst)),as.numeric(t(genedif_GL_C_NO_C_sub_meangst))))
GL_C_NO_C_sub_meangst<-rep(0:854,length.out=22230)
table_GL_C_NO_C_sub_meangs<-new.genedif_GL_C_NO_C_sub_meangst
table_GL_C_NO_C_sub_meangs$V2<-bins
SNP_table_GL_C_NO_C_sub_meangs<-c(1:(nrow(table_GL_C_NO_C_sub_meangs)))
df_table_GL_C_NO_C_sub_meangs<-data.frame(SNP_table_GL_C_NO_C_sub_meangs,table_GL_C_NO_C_sub_meangs)
colnames(df_table_GL_C_NO_C_sub_meangs) <- c("SNP","CHROM", "bin","Gst_mean")
df_table_GL_C_NO_C_sub_meangs$CHROM <- revalue(df_table_GL_C_NO_C_sub_meangs$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
df_table_GL_C_NO_C_sub_meangs$CHROM <- as.numeric(df_table_GL_C_NO_C_sub_meangs$CHROM)
df_table_GL_C_NO_C_sub_meangs$Gst_mean[df_table_GL_C_NO_C_sub_meangs$Gst_mean == NA] <- NA
df_table_GL_C_NO_C_sub_meangs_sub<-df_table_GL_C_NO_C_sub_meangs[complete.cases(df_table_GL_C_NO_C_sub_meangs),]
df_table_GL_C_NO_C_sub_meangs_sub$Gst_mean <- as.numeric(df_table_GL_C_NO_C_sub_meangs_sub$Gst_mean)
manhattan(df_table_GL_C_NO_C_sub_meangs_sub, chr="CHROM", bp="bin", p="Gst_mean", snp="SNP", logp=FALSE, ylab="Nei´s Gst")

# Gst along the genome IS_C_IS_H
genedif_GL_C_GL_H<-read.table("Gst_vcfR/Gst_genedif_IS_C_IS_H.txt", header=TRUE)
genedif_IS_C_IS_H_sub<-genedif_IS_C_IS_H[complete.cases(genedif_IS_C_IS_H),]
genedif_IS_C_IS_H_sub$bin<-genedif_IS_C_IS_H_sub$bin <- as.numeric(genedif_IS_C_IS_H_sub$POS)%/%100000
genedif_IS_C_IS_H_sub_meangst<-tapply(genedif_IS_C_IS_H_sub$Gst, list(genedif_IS_C_IS_H_sub$CHROM, genedif_IS_C_IS_H_sub$bin), mean)
new.genedif_IS_C_IS_H_sub_meangst=as.data.frame(cbind(rep(row.names(genedif_IS_C_IS_H_sub_meangst),each=855),rep(0:854,each=nrow(genedif_IS_C_IS_H_sub_meangst)),as.numeric(t(genedif_IS_C_IS_H_sub_meangst))))
IS_C_IS_H_sub_meangst<-rep(0:854,length.out=22230)
table_IS_C_IS_H_sub_meangs<-new.genedif_IS_C_IS_H_sub_meangst
table_IS_C_IS_H_sub_meangs$V2<-bins
SNP_table_IS_C_IS_H_sub_meangs<-c(1:(nrow(table_IS_C_IS_H_sub_meangs)))
df_table_IS_C_IS_H_sub_meangs<-data.frame(SNP_table_IS_C_IS_H_sub_meangs,table_IS_C_IS_H_sub_meangs)
colnames(df_table_IS_C_IS_H_sub_meangs) <- c("SNP","CHROM", "bin","Gst_mean")
df_table_IS_C_IS_H_sub_meangs$CHROM <- revalue(df_table_IS_C_IS_H_sub_meangs$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
df_table_IS_C_IS_H_sub_meangs$CHROM <- as.numeric(df_table_IS_C_IS_H_sub_meangs$CHROM)
df_table_IS_C_IS_H_sub_meangs$Gst_mean[df_table_IS_C_IS_H_sub_meangs$Gst_mean == NA] <- NA
df_table_IS_C_IS_H_sub_meangs_sub<-df_table_IS_C_IS_H_sub_meangs[complete.cases(df_table_IS_C_IS_H_sub_meangs),]
df_table_IS_C_IS_H_sub_meangs_sub$Gst_mean <- as.numeric(df_table_IS_C_IS_H_sub_meangs_sub$Gst_mean)
manhattan(df_table_IS_C_IS_H_sub_meangs_sub, chr="CHROM", bp="bin", p="Gst_mean", snp="SNP", logp=FALSE, ylab="Nei´s Gst")

tail(sort(df_table_IS_C_IS_H_sub_meangs_sub$Gst_mean), n=5)
df_table_IS_C_IS_H_sub_meangs_sub %>% filter(Gst_mean>0.6695735)
genedif_IS_C_IS_H_sub %>% filter(CHROM=="LR606182.1") %>% filter(bin == 312 | bin == 320) 
# 312 contains NUBPL, nucleotide binding protein like
# 320 contains PRKD1, protein kinase D1 
genedif_IS_C_IS_H_sub %>% filter(CHROM=="LR606196.1") %>% filter(bin == 161 | bin == 162 | bin == 181) 
# 161 contains CALCB, calcitonin related polypeptide alpha and CYP2R1, cytochrome P450 family 2 subfamily R member 1
# 162 just behind 161 contains PDE3B, phosphodiesterase 3B
# 181 contains ADM, adrenomedullin

# Gst along the genome IS_C_NO_C
genedif_GL_C_GL_H<-read.table("Gst_vcfR/Gst_genedif_IS_C_NO_C.txt", header=TRUE)
genedif_IS_C_NO_C_sub<-genedif_IS_C_NO_C[complete.cases(genedif_IS_C_NO_C),]
genedif_IS_C_NO_C_sub$bin<-genedif_IS_C_NO_C_sub$bin <- as.numeric(genedif_IS_C_NO_C_sub$POS)%/%100000
genedif_IS_C_NO_C_sub_meangst<-tapply(genedif_IS_C_NO_C_sub$Gst, list(genedif_IS_C_NO_C_sub$CHROM, genedif_IS_C_NO_C_sub$bin), mean)
new.genedif_IS_C_NO_C_sub_meangst=as.data.frame(cbind(rep(row.names(genedif_IS_C_NO_C_sub_meangst),each=855),rep(0:854,each=nrow(genedif_IS_C_NO_C_sub_meangst)),as.numeric(t(genedif_IS_C_NO_C_sub_meangst))))
IS_C_NO_C_sub_meangst<-rep(0:854,length.out=22230)
table_IS_C_NO_C_sub_meangs<-new.genedif_IS_C_NO_C_sub_meangst
table_IS_C_NO_C_sub_meangs$V2<-bins
SNP_table_IS_C_NO_C_sub_meangs<-c(1:(nrow(table_IS_C_NO_C_sub_meangs)))
df_table_IS_C_NO_C_sub_meangs<-data.frame(SNP_table_IS_C_NO_C_sub_meangs,table_IS_C_NO_C_sub_meangs)
colnames(df_table_IS_C_NO_C_sub_meangs) <- c("SNP","CHROM", "bin","Gst_mean")
df_table_IS_C_NO_C_sub_meangs$CHROM <- revalue(df_table_IS_C_NO_C_sub_meangs$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
df_table_IS_C_NO_C_sub_meangs$CHROM <- as.numeric(df_table_IS_C_NO_C_sub_meangs$CHROM)
df_table_IS_C_NO_C_sub_meangs$Gst_mean[df_table_IS_C_NO_C_sub_meangs$Gst_mean == NA] <- NA
df_table_IS_C_NO_C_sub_meangs_sub<-df_table_IS_C_NO_C_sub_meangs[complete.cases(df_table_IS_C_NO_C_sub_meangs),]
df_table_IS_C_NO_C_sub_meangs_sub$Gst_mean <- as.numeric(df_table_IS_C_NO_C_sub_meangs_sub$Gst_mean)
manhattan(df_table_IS_C_NO_C_sub_meangs_sub, chr="CHROM", bp="bin", p="Gst_mean", snp="SNP", logp=FALSE, ylab="Nei´s Gst")

# Gst along the genome NO_C_NO_H
genedif_GL_C_GL_H<-read.table("Gst_vcfR/Gst_genedif_NO_C_NO_H.txt", header=TRUE)
genedif_NO_C_NO_H_sub<-genedif_NO_C_NO_H[complete.cases(genedif_NO_C_NO_H),]
genedif_NO_C_NO_H_sub$bin<-genedif_NO_C_NO_H_sub$bin <- as.numeric(genedif_NO_C_NO_H_sub$POS)%/%100000
genedif_NO_C_NO_H_sub_meangst<-tapply(genedif_NO_C_NO_H_sub$Gst, list(genedif_NO_C_NO_H_sub$CHROM, genedif_NO_C_NO_H_sub$bin), mean)
new.genedif_NO_C_NO_H_sub_meangst=as.data.frame(cbind(rep(row.names(genedif_NO_C_NO_H_sub_meangst),each=855),rep(0:854,each=nrow(genedif_NO_C_NO_H_sub_meangst)),as.numeric(t(genedif_NO_C_NO_H_sub_meangst))))
NO_C_NO_H_sub_meangst<-rep(0:854,length.out=22230)
table_NO_C_NO_H_sub_meangs<-new.genedif_NO_C_NO_H_sub_meangst
table_NO_C_NO_H_sub_meangs$V2<-bins
SNP_table_NO_C_NO_H_sub_meangs<-c(1:(nrow(table_NO_C_NO_H_sub_meangs)))
df_table_NO_C_NO_H_sub_meangs<-data.frame(SNP_table_NO_C_NO_H_sub_meangs,table_NO_C_NO_H_sub_meangs)
colnames(df_table_NO_C_NO_H_sub_meangs) <- c("SNP","CHROM", "bin","Gst_mean")
df_table_NO_C_NO_H_sub_meangs$CHROM <- revalue(df_table_NO_C_NO_H_sub_meangs$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
df_table_NO_C_NO_H_sub_meangs$CHROM <- as.numeric(df_table_NO_C_NO_H_sub_meangs$CHROM)
df_table_NO_C_NO_H_sub_meangs$Gst_mean[df_table_NO_C_NO_H_sub_meangs$Gst_mean == NA] <- NA
df_table_NO_C_NO_H_sub_meangs_sub<-df_table_NO_C_NO_H_sub_meangs[complete.cases(df_table_NO_C_NO_H_sub_meangs),]
df_table_NO_C_NO_H_sub_meangs_sub$Gst_mean <- as.numeric(df_table_NO_C_NO_H_sub_meangs_sub$Gst_mean)
manhattan(df_table_NO_C_NO_H_sub_meangs_sub, chr="CHROM", bp="bin", p="Gst_mean", snp="SNP", logp=FALSE, ylab="Nei´s Gst")

# Gst manhattens, all in one 
#pdf("Gst_mean_along_chromosome_6pops_IS_GL_NO_030621.pdf", width=13, height=8)
par(mfrow=c(6,1))
par(mar=c(2.5,3,0.5,0.5))
manhattan(df_table_GL_C_GL_H_sub_meangs_sub, chr="CHROM", bp="bin", p="Gst_mean", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5)
title("A", adj=0.02, line = -1)
manhattan(df_table_GL_C_IS_C_sub_meangs_sub, chr="CHROM", bp="bin", p="Gst_mean", snp="SNP", logp=FALSE,ylab="", xlab="", cex.axis=1.5)
title("B", adj=0.02, line = -1)
manhattan(df_table_GL_C_NO_C_sub_meangs_sub, chr="CHROM", bp="bin", p="Gst_mean", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5)
title("C", adj=0.02, line = -1)
manhattan(df_table_IS_C_IS_H_sub_meangs_sub, chr="CHROM", bp="bin", p="Gst_mean", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5)
title("D", adj=0.02, line = -1)
manhattan(df_table_IS_C_NO_C_sub_meangs_sub, chr="CHROM", bp="bin", p="Gst_mean", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5)
title("E", adj=0.02, line = -1)
par(mar=c(2.5,3,0.5,0.5))
manhattan(df_table_NO_C_NO_H_sub_meangs_sub, chr="CHROM", bp="bin", p="Gst_mean", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5)
title("F", adj=0.02, line = -1)
#dev.off()

tail(sort(df_table_GL_C_GL_H_sub_meangs_sub$Gst_mean), n=10)
mean(df_table_GL_C_GL_H_sub_meangs_sub$Gst_mean)
quantile(sort(df_table_GL_C_GL_H_sub_meangs_sub$Gst_mean), 0.99)
quantile(sort(df_table_GL_C_GL_H_sub_meangs_sub$Gst_mean), 0.999)
df_table_GL_C_GL_H_sub_meangs_sub %>% filter(Gst_mean>0.2178788)
df_table_GL_C_GL_H_sub_meangs_sub %>% filter(Gst_mean>0.233976)
genedif_GL_C_GL_H_sub %>% filter(CHROM=="LR606182.1") %>% filter(bin == 621 )# LR606182.1 sites 62118628 - 62190432 
genedif_GL_C_GL_H_sub %>% filter(CHROM=="LR606182.1") %>% filter(bin == 704 )# LR606182.1 sites 70403430 - 70491063
genedif_GL_C_GL_H_sub %>% filter(CHROM=="LR606182.1") %>% filter(bin == 726 )# LR606182.1 sites 72652540 - 72664623 
genedif_GL_C_GL_H_sub %>% filter(CHROM=="LR606187.1") %>% filter(bin == 230 )# LR606187.1 sites 23001409 - 23094257 
genedif_GL_C_GL_H_sub %>% filter(CHROM=="LR606206.1") %>% filter(bin == 153 )# LR606206.1 sites 15300071 - 15377626
# [1] 0.2178788 0.2204283 0.2218441 0.2221453 0.2292619 0.2339763 0.2450859 0.2479435 0.2533147 0.2769024
# 0.03424935
# 0.1507834, 109 windows are in the top 1% percentile
# 0.2168228, 11 
tail(sort(df_table_GL_C_IS_C_sub_meangs_sub$Gst_mean), n=10)
mean(df_table_GL_C_IS_C_sub_meangs_sub$Gst_mean)
quantile(sort(df_table_GL_C_IS_C_sub_meangs_sub$Gst_mean), 0.99)
df_table_GL_C_IS_C_sub_meangs_sub %>% filter(Gst_mean>0.2290198)
df_table_GL_C_IS_C_sub_meangs_sub %>% filter(Gst_mean>0.252007)
genedif_GL_C_IS_C_sub %>% filter(CHROM=="LR606183.1") %>% filter(bin == 182 )# LR606183.1 sites 18212785 - 18215171
genedif_GL_C_IS_C_sub %>% filter(CHROM=="LR606183.1") %>% filter(bin == 228 )# LR606183.1 sites 22816578 - 22892621
genedif_GL_C_IS_C_sub %>% filter(CHROM=="LR606183.1") %>% filter(bin == 229 )# LR606183.1 sites 22901488 - 22994858
genedif_GL_C_IS_C_sub %>% filter(CHROM=="LR606183.1") %>% filter(bin == 231 )# LR606183.1 sites 23130506 - 23182772
genedif_GL_C_IS_C_sub %>% filter(CHROM=="LR606193.1") %>% filter(bin == 289 )# LR606193.1 sites 28908817 - 28992270
# [1] 0.2290198 0.2302330 0.2324133 0.2345182 0.2444082 0.2520074 0.2558147 0.2593828 0.2678723 0.2749365
# [1] 0.03745914
# 0.1406353, 110 windows are in the top 1% percentile
tail(sort(df_table_GL_C_NO_C_sub_meangs_sub$Gst_mean), n=10)
mean(df_table_GL_C_NO_C_sub_meangs_sub$Gst_mean)
quantile(sort(df_table_GL_C_NO_C_sub_meangs_sub$Gst_mean), 0.99)
df_table_GL_C_NO_C_sub_meangs_sub %>% filter(Gst_mean>0.5597770)
df_table_GL_C_NO_C_sub_meangs_sub %>% filter(Gst_mean>0.615926)
genedif_GL_C_NO_C_sub %>% filter(CHROM=="LR606181.1") %>% filter(bin == 582 )# LR606181.1 sites 58200511 - 58298449 
genedif_GL_C_NO_C_sub %>% filter(CHROM=="LR606187.1") %>% filter(bin == 247 )# LR606187.1 sites 24706860 - 24797653 
genedif_GL_C_NO_C_sub %>% filter(CHROM=="LR606193.1") %>% filter(bin == 209 )# LR606193.1 sites 20911001 - 20998861
genedif_GL_C_NO_C_sub %>% filter(CHROM=="LR606193.1") %>% filter(bin == 385 )# LR606193.1 sites 38500440 - 38597429
genedif_GL_C_NO_C_sub %>% filter(CHROM=="LR606196.1") %>% filter(bin == 91  )# LR606196.1 sites 9115417  - 9198779 
# [1] 0.5597770 0.5630246 0.5951097 0.6107685 0.6118879 0.6159263 0.6184662 0.6396692 0.6808996 0.7839295
# 0.1747237
# 0.4035339, 111 windows are in the top 1% percentile
tail(sort(df_table_IS_C_IS_H_sub_meangs_sub$Gst_mean), n=10)
mean(df_table_IS_C_IS_H_sub_meangs_sub$Gst_mean)
quantile(sort(df_table_IS_C_IS_H_sub_meangs_sub$Gst_mean), 0.99)
df_table_IS_C_IS_H_sub_meangs_sub %>% filter(Gst_mean>0.5112795)
df_table_IS_C_IS_H_sub_meangs_sub %>% filter(Gst_mean>0.669573)
genedif_IS_C_IS_H_sub %>% filter(CHROM=="LR606182.1") %>% filter(bin == 312 )# LR606182.1 sites 31249530 - 31282461
genedif_IS_C_IS_H_sub %>% filter(CHROM=="LR606182.1") %>% filter(bin == 320 )# LR606182.1 sites 32013738 - 32085671
genedif_IS_C_IS_H_sub %>% filter(CHROM=="LR606196.1") %>% filter(bin == 161 )# LR606196.1 sites 16112794 - 16156396
genedif_IS_C_IS_H_sub %>% filter(CHROM=="LR606196.1") %>% filter(bin == 162 )# LR606196.1 sites 16248560 - 16281711
genedif_IS_C_IS_H_sub %>% filter(CHROM=="LR606196.1") %>% filter(bin == 181 )# LR606196.1 sites 18124721 - 18171544 
# [1] 0.5112795 0.5307581 0.5400703 0.5766690 0.6527778 0.6695736 0.7448980 0.7959984 0.8191862 0.8392328
# 0.04466201
# 0.417714, 106 windows are in the top 1% percentile
tail(sort(df_table_IS_C_NO_C_sub_meangs_sub$Gst_mean), n=10)
mean(df_table_IS_C_NO_C_sub_meangs_sub$Gst_mean)
quantile(sort(df_table_IS_C_NO_C_sub_meangs_sub$Gst_mean), 0.99)
df_table_IS_C_NO_C_sub_meangs_sub %>% filter(Gst_mean>0.6041078)
df_table_IS_C_NO_C_sub_meangs_sub %>% filter(Gst_mean>0.635071)
genedif_IS_C_NO_C_sub %>% filter(CHROM=="LR606185.1") %>% filter(bin == 367 )# LR606185.1 sites 36702960 - 36796730
genedif_IS_C_NO_C_sub %>% filter(CHROM=="LR606187.1") %>% filter(bin == 247 )# LR606187.1 sites 24706860 - 24797653
genedif_IS_C_NO_C_sub %>% filter(CHROM=="LR606192.1") %>% filter(bin == 256 )# LR606192.1 sites 25601336 - 25690140
genedif_IS_C_NO_C_sub %>% filter(CHROM=="LR606196.1") %>% filter(bin == 91  )# LR606196.1 sites 9115417  - 9198779
genedif_IS_C_NO_C_sub %>% filter(CHROM=="LR606203.1") %>% filter(bin == 2   )# LR606203.1 sites 203475   - 296826
# [1] 0.6041078 0.6080262 0.6161049 0.6306242 0.6328931 0.6350713 0.6355762 0.6454940 0.6528318 0.6581610
# 0.181222
# 0.441, 150 windows are in the top 1% percentile
tail(sort(df_table_NO_C_NO_H_sub_meangs_sub$Gst_mean), n=10)
mean(df_table_NO_C_NO_H_sub_meangs_sub$Gst_mean)
quantile(sort(df_table_NO_C_NO_H_sub_meangs_sub$Gst_mean), 0.99)
df_table_NO_C_NO_H_sub_meangs_sub %>% filter(Gst_mean>0.1463457)
df_table_NO_C_NO_H_sub_meangs_sub %>% filter(Gst_mean>0.167242)
genedif_NO_C_NO_H_sub %>% filter(CHROM=="LR606181.1") %>% filter(bin == 582 )# LR606181.1 sites 58200511 - 58298449 
genedif_NO_C_NO_H_sub %>% filter(CHROM=="LR606182.1") %>% filter(bin == 267 )# LR606182.1 sites 26703706 - 26799250 
genedif_NO_C_NO_H_sub %>% filter(CHROM=="LR606183.1") %>% filter(bin == 557 )# LR606183.1 sites 55700166 - 55785664
genedif_NO_C_NO_H_sub %>% filter(CHROM=="LR606185.1") %>% filter(bin == 159 )# LR606185.1 sites 15906196 - 15990426 
genedif_NO_C_NO_H_sub %>% filter(CHROM=="LR606193.1") %>% filter(bin == 224 )# LR606193.1 sites 22405370 - 22472063 
# [1] 0.1463457 0.1544603 0.1571986 0.1606028 0.1667754 0.1672422 0.1690192 0.1859624 0.2229830 0.2251594
# 0.04048088
# 0.09510602, 111 windows are in the top 1% percentile

# for(i in 1:tail(mydf$bin)[1]){
#   mean(mydf$Gst[mydf$bin==[i]])
# }
# 
# mydf$meanGst <- mean(mydf$Gst[mydf$bin])
# mydf$sdGst <- sd(mydf$Gst[mydf$bin==177])
# manhattan(mydf, chr="CHROM", bp="POS", p="meanGst", snp="SNP", logp=FALSE, ylab="Nei´s Gst")


c(quantile(na.omit(genedif_GL_C_GL_H$Gst),0.025))
c(quantile(na.omit(genedif_GL_C_GL_H$Gst),0.5))
c(quantile(na.omit(genedif_GL_C_GL_H$Gst),0.975))
length(na.omit(genedif_GL_C_GL_H$Gst))
length(na.omit(genedif_GL_C_GL_H$Gst[genedif_GL_C_GL_H$Gst>0.1628959 | genedif_GL_C_GL_H$Gst<0]))
length(na.omit(genedif_GL_C_GL_H$Gst[genedif_GL_C_GL_H$Gst>0.1628959 | genedif_GL_C_GL_H$Gst<0]))/
  length(na.omit(genedif_GL_C_GL_H$Gst))

c(quantile(na.omit(genedif_GL_C_NO_H$Gst),0.025))
c(quantile(na.omit(genedif_GL_C_NO_H$Gst),0.5))
c(quantile(na.omit(genedif_GL_C_NO_H$Gst),0.975))
length(na.omit(genedif_GL_C_GL_H$Gst))
length(na.omit(genedif_GL_C_NO_H$Gst[genedif_GL_C_NO_H$Gst>0.7272727  | genedif_GL_C_GL_H$Gst<0.0002020202 ]))
length(na.omit(genedif_GL_C_NO_H$Gst[genedif_GL_C_NO_H$Gst>0.7272727  | genedif_GL_C_GL_H$Gst<0.0002020202 ]))/
  length(na.omit(genedif_GL_C_NO_H$Gst))

length(na.omit(genedif_GL_C_NO_H$Gst[genedif_GL_C_NO_H$Gst>0.7272727  | genedif_GL_C_GL_H$Gst<0.0002020202 ]))/
  +   length(na.omit(genedif_GL_C_NO_H$Gst))

quantile(genedif_GL_C_GL_H$Gst, c(0.025, 0.975),na.rm=T)

?fisher.test()         

#### Fst in windows 100K
fst_win_GL_C_GL_H<-read.table("Fst_windows_vcftools/fst_weir_GLH_GLM_windows100k_230821.windowed.weir.fst", header=T)
fst_win_GL_C_IS_C<-read.table("Fst_windows_vcftools/fst_weir_GLM_ISM_windows100k_230821.windowed.weir.fst", header=T)
fst_win_GL_C_NO_C<-read.table("Fst_windows_vcftools/fst_weir_GLM_NOM_windows100k_230821.windowed.weir.fst", header=T)
fst_win_IS_C_IS_H<-read.table("Fst_windows_vcftools/fst_weir_ISM_ISH_windows100k_230821.windowed.weir.fst", header=T)
fst_win_IS_C_NO_C<-read.table("Fst_windows_vcftools/fst_weir_ISM_NOM_windows100k_230821.windowed.weir.fst", header=T)
fst_win_NO_C_NO_H<-read.table("Fst_windows_vcftools/fst_weir_NOM_NOH_windows100k_230821.windowed.weir.fst", header=T)

fst_win_GL_C_GL_H$CHROM <- revalue(fst_win_GL_C_GL_H$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
fst_win_GL_C_IS_C$CHROM <- revalue(fst_win_GL_C_IS_C$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
fst_win_GL_C_NO_C$CHROM <- revalue(fst_win_GL_C_NO_C$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
fst_win_IS_C_IS_H$CHROM <- revalue(fst_win_IS_C_IS_H$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
fst_win_IS_C_NO_C$CHROM <- revalue(fst_win_IS_C_NO_C$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))
fst_win_NO_C_NO_H$CHROM <- revalue(fst_win_NO_C_NO_H$CHROM, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26"))

fst_win_GL_C_GL_H$CHROM <- as.numeric(fst_win_GL_C_GL_H$CHROM)
fst_win_GL_C_IS_C$CHROM <- as.numeric(fst_win_GL_C_IS_C$CHROM)
fst_win_GL_C_NO_C$CHROM <- as.numeric(fst_win_GL_C_NO_C$CHROM)
fst_win_IS_C_IS_H$CHROM <- as.numeric(fst_win_IS_C_IS_H$CHROM)
fst_win_IS_C_NO_C$CHROM <- as.numeric(fst_win_IS_C_NO_C$CHROM)
fst_win_NO_C_NO_H$CHROM <- as.numeric(fst_win_NO_C_NO_H$CHROM)

fst_win_GL_C_GL_H$SNP<-seq.int(nrow(fst_win_GL_C_GL_H))
fst_win_GL_C_IS_C$SNP<-seq.int(nrow(fst_win_GL_C_IS_C))
fst_win_GL_C_NO_C$SNP<-seq.int(nrow(fst_win_GL_C_NO_C))
fst_win_IS_C_IS_H$SNP<-seq.int(nrow(fst_win_IS_C_IS_H))
fst_win_IS_C_NO_C$SNP<-seq.int(nrow(fst_win_IS_C_NO_C))
fst_win_NO_C_NO_H$SNP<-seq.int(nrow(fst_win_NO_C_NO_H))

pdf("Fst_mean_along_chromosome_6pops_IS_GL_NO_240821.pdf", width=13, height=8)
par(mfrow=c(6,1))
par(mar=c(2.5,3,0.5,0.5))
manhattan(fst_win_GL_C_GL_H, chr="CHROM", bp="BIN_START", p="WEIGHTED_FST", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5)
title("A", adj=0.02, line = -1)
manhattan(fst_win_GL_C_IS_C, chr="CHROM", bp="BIN_START", p="WEIGHTED_FST", snp="SNP", logp=FALSE,ylab="", xlab="", cex.axis=1.5)
title("B", adj=0.02, line = -1)
manhattan(fst_win_GL_C_NO_C, chr="CHROM", bp="BIN_START", p="WEIGHTED_FST", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5)
title("C", adj=0.02, line = -1)
manhattan(fst_win_IS_C_IS_H, chr="CHROM", bp="BIN_START", p="WEIGHTED_FST", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5)
title("D", adj=0.02, line = -1)
manhattan(fst_win_IS_C_NO_C, chr="CHROM", bp="BIN_START", p="WEIGHTED_FST", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5)
title("E", adj=0.02, line = -1)
par(mar=c(2.5,3,0.5,0.5))
manhattan(fst_win_NO_C_NO_H, chr="CHROM", bp="BIN_START", p="WEIGHTED_FST", snp="SNP", logp=FALSE, ylab="", xlab="", cex.axis=1.5)
title("F", adj=0.02, line = -1)
dev.off()

mean(fst_win_GL_C_GL_H$WEIGHTED_FST)
mean(fst_win_GL_C_IS_C$WEIGHTED_FST)
mean(fst_win_GL_C_NO_C$WEIGHTED_FST)
mean(fst_win_IS_C_IS_H$WEIGHTED_FST)
mean(fst_win_IS_C_NO_C$WEIGHTED_FST)
mean(fst_win_NO_C_NO_H$WEIGHTED_FST)


fst_win_IS_C_IS_H %>% filter(WEIGHTED_FST>0.9285)
# 2	32000001	32100000 # PRKD
# 16	16100001	16200000
# 16	16200001	16300000 # CYP2R1, CALCB, PDE3B
# 15	21600001	21700000
# 7	21100001	21200000
# 17	9700001	9800000
# 16	18100001	18200000
# 4	43700001	43800000
# 5	17500001	17600000
# 6	12000001	12100000
# 6	12400001	12500000
# 6	13300001	13400000
# 6	20000001	20100000
# 6	20200001	20300000
# 6	21100001	21200000
# 7	16800001	16900000 ##  GJA8, HTR1F, PDZK1, POUI1F1, VGLL3, CHMP2B, CADM 
# 7	16900001	17000000
# 7	17500001	17600000
# 7	18200001	18300000
# 7	18300001	18400000

fst_win_GL_C_IS_C %>% filter(WEIGHTED_FST>0.90)
# LR606192.1:4600001-4700000
# LR606192.1:4800001-4900000
# LR606196.1:16200001-16300000
# LR606187.1:20400001-23100000 # super mange ROBO1 and ROBO2, and , RBM11, HSPA13, SAMSN1, NRIP1,
# LR606182.1:32000001-32100000

# 7	20500001	20600000
# 7	20900001	21000000
# 7	21600001	21700000
# 7	21300001	21400000
# 7	21000001	21100000
# 7	21400001	21500000
# 7	21100001	21200000
# 7	21200001	21300000
# 16	16200001	16300000
# 7	23000001	23100000
# 7	20600001	20700000
# 7	20400001	20500000
# 7	21700001	21800000
# 12	4800001	4900000
# 7	20700001	20800000
# 7	22000001	22100000
# 7	22900001	23000000
# 7	20800001	20900000
# 2	32000001	32100000
# 7	21900001	22000000
# 7	21800001	21900000
# 12	4600001	4700000
# 7	22200001	22300000

  
### heterozygosity per pos
hardy_GL_C<-read.table(text = gsub("/", "\t", readLines("hardy_het_per_site/hardy_GL_C_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe")), header = T)
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


hardy_GL_C<-read.table(text = gsub("/", "\t", readLines("hardy_het_per_site/hardy_GL_C_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe")), header = T)
hardy_GL_C$hetero<-hardy_GL_C$HET/(hardy_GL_C$OBS.HOM1+hardy_GL_C$HET+hardy_GL_C$HOM2.)
hardy_GL_C$hetero_exp<-hardy_GL_C$HET.1/(hardy_GL_C$E.HOM1+hardy_GL_C$HET.1+hardy_GL_C$HOM2..1)
hardy_GL_C$CHR <- as.numeric(revalue(hardy_GL_C$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_GL_H<-read.table(text = gsub("/", "\t", readLines("hardy_het_per_site/hardy_GL_H_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe")), header = T)
hardy_GL_H$hetero<-hardy_GL_H$HET/(hardy_GL_H$OBS.HOM1+hardy_GL_H$HET+hardy_GL_H$HOM2.)
hardy_GL_H$hetero_exp<-hardy_GL_H$HET.1/(hardy_GL_H$E.HOM1+hardy_GL_H$HET.1+hardy_GL_H$HOM2..1)
hardy_GL_H$CHR <- as.numeric(revalue(hardy_GL_H$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_IS_C<-read.table(text = gsub("/", "\t", readLines("hardy_het_per_site/hardy_IS_C_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe")), header = T)
hardy_IS_C$hetero<-hardy_IS_C$HET/(hardy_IS_C$OBS.HOM1+hardy_IS_C$HET+hardy_IS_C$HOM2.)
hardy_IS_C$hetero_exp<-hardy_IS_C$HET.1/(hardy_IS_C$E.HOM1+hardy_IS_C$HET.1+hardy_IS_C$HOM2..1)
hardy_IS_C$CHR <- as.numeric(revalue(hardy_IS_C$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_IS_H<-read.table(text = gsub("/", "\t", readLines("hardy_het_per_site/hardy_IS_H_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe")), header = T)
hardy_IS_H$hetero<-hardy_IS_H$HET/(hardy_IS_H$OBS.HOM1+hardy_IS_H$HET+hardy_IS_H$HOM2.)
hardy_IS_H$hetero_exp<-hardy_IS_H$HET.1/(hardy_IS_H$E.HOM1+hardy_IS_H$HET.1+hardy_IS_H$HOM2..1)
hardy_IS_H$CHR <- as.numeric(revalue(hardy_IS_H$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_NO_C<-read.table(text = gsub("/", "\t", readLines("hardy_het_per_site/hardy_NO_C_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe")), header = T)
hardy_NO_C$hetero<-hardy_NO_C$HET/(hardy_NO_C$OBS.HOM1+hardy_NO_C$HET+hardy_NO_C$HOM2.)
hardy_NO_C$hetero_exp<-hardy_NO_C$HET.1/(hardy_NO_C$E.HOM1+hardy_NO_C$HET.1+hardy_NO_C$HOM2..1)
hardy_NO_C$CHR <- as.numeric(revalue(hardy_NO_C$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_NO_H<-read.table(text = gsub("/", "\t", readLines("hardy_het_per_site/hardy_NO_H_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe")), header = T)
hardy_NO_H$hetero<-hardy_NO_H$HET/(hardy_NO_H$OBS.HOM1+hardy_NO_H$HET+hardy_NO_H$HOM2.)
hardy_NO_H$hetero_exp<-hardy_NO_H$HET.1/(hardy_NO_H$E.HOM1+hardy_NO_H$HET.1+hardy_NO_H$HOM2..1)
hardy_NO_H$CHR <- as.numeric(revalue(hardy_NO_H$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_DK_C<-read.table(text = gsub("/", "\t", readLines("hardy_het_per_site/hardy_DK_C_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe")), header = T)
hardy_DK_C$hetero<-hardy_DK_C$HET/(hardy_DK_C$OBS.HOM1+hardy_DK_C$HET+hardy_DK_C$HOM2.)
hardy_DK_C$hetero_exp<-hardy_DK_C$HET.1/(hardy_DK_C$E.HOM1+hardy_DK_C$HET.1+hardy_DK_C$HOM2..1)
hardy_DK_C$CHR <- as.numeric(revalue(hardy_DK_C$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_DK_H<-read.table(text = gsub("/", "\t", readLines("hardy_het_per_site/hardy_DK_H_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe")), header = T)
hardy_DK_H$hetero<-hardy_DK_H$HET/(hardy_DK_H$OBS.HOM1+hardy_DK_H$HET+hardy_DK_H$HOM2.)
hardy_DK_H$hetero_exp<-hardy_DK_H$HET.1/(hardy_DK_H$E.HOM1+hardy_DK_H$HET.1+hardy_DK_H$HOM2..1)
hardy_DK_H$CHR <- as.numeric(revalue(hardy_DK_H$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_EE_C<-read.table(text = gsub("/", "\t", readLines("hardy_het_per_site/hardy_EE_C_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe")), header = T)
hardy_EE_C$hetero<-hardy_EE_C$HET/(hardy_EE_C$OBS.HOM1+hardy_EE_C$HET+hardy_EE_C$HOM2.)
hardy_EE_C$hetero_exp<-hardy_EE_C$HET.1/(hardy_EE_C$E.HOM1+hardy_EE_C$HET.1+hardy_EE_C$HOM2..1)
hardy_EE_C$CHR <- as.numeric(revalue(hardy_EE_C$CHR, c("LR606181.1"="1", "LR606182.1"="2", "LR606183.1"="3", "LR606184.1"="4", "LR606185.1"="5", "LR606186.1"="6", "LR606187.1"="7", "LR606188.1"="8", "LR606189.1"="9", "LR606190.1"="10", "LR606191.1"="11", "LR606192.1"="12", "LR606193.1"="13", "LR606194.1"="14", "LR606195.1"="15", "LR606196.1"="16", "LR606197.1"="17", "LR606198.1"="18", "LR606199.1"="19", "LR606200.1"="20", "LR606201.1"="21", "LR606202.1"="22", "LR606203.1"="23", "LR606204.1"="24", "LR606205.1"="25", "LR606206.1"="26")))
hardy_TU_H<-read.table(text = gsub("/", "\t", readLines("hardy_het_per_site/hardy_TU_H_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe")), header = T)
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
#dev.off()

plot(ecdf(sort(log(na.omit(hardy_IS_C[,"P_HWE"])))), do.points = FALSE)
lines(ecdf(log(hardy_IS_C[,"P_HWE"])), do.points = FALSE)

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

pdf("pvalues_hwe_heterozygosity_allpops_250821.pdf", width = 13, height = 8)
plot(x=sort(log(hardy_GL_C$P_HWE)), y=cumsum_hardy_GL_C/length(cumsum_hardy_GL_C), type="l", lwd = 2, col="black", xlab="log(HWE p-value)", ylab="cumulative proportion", cex.lab=1.5, cex.axis=1.5)
lines(x=sort(log(hardy_GL_H$P_HWE)), y=cumsum_hardy_GL_H/length(cumsum_hardy_GL_H), col="black", lwd = 2, lty=2)
lines(x=sort(log(hardy_IS_C$P_HWE)), y=cumsum_hardy_IS_C/length(cumsum_hardy_IS_C), col="grey", lwd = 2)
lines(x=sort(log(hardy_IS_H$P_HWE)), y=cumsum_hardy_IS_H/length(cumsum_hardy_IS_H), col="grey", lwd = 2 , lty=2)
lines(x=sort(log(hardy_NO_C$P_HWE)), y=cumsum_hardy_NO_C/length(cumsum_hardy_NO_C), col="green", lwd = 2)
lines(x=sort(log(hardy_NO_H$P_HWE)), y=cumsum_hardy_NO_H/length(cumsum_hardy_NO_H), col="green", lwd = 2 , lty=2)
lines(x=sort(log(hardy_DK_C$P_HWE)), y=cumsum_hardy_DK_C/length(cumsum_hardy_DK_C), col="red", lwd = 2)
lines(x=sort(log(hardy_DK_H$P_HWE)), y=cumsum_hardy_DK_H/length(cumsum_hardy_DK_H), col="red", lwd = 2 , lty=2)
lines(x=sort(log(hardy_EE_C$P_HWE)), y=cumsum_hardy_EE_C/length(cumsum_hardy_EE_C), col="blue", lwd = 2)
lines(x=sort(log(hardy_TU_H$P_HWE)), y=cumsum_hardy_TU_H/length(cumsum_hardy_TU_H), col="yellow", lwd = 2, lty=2)
legend("topleft",        # Add legend to plot
       legend = c("GL", "IS", "NO", "DK", "EE", "TU"),
       col = c("black", "grey", "green","red", "blue", "yellow"),
       pch = 16,
       cex = 1.5)
dev.off()

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

hardy_IS_C_clean<-completeFun(hardy_IS_C, "hetero")
hardy_IS_H_clean<-completeFun(hardy_IS_H, "hetero")
hardy_DK_H_clean<-completeFun(hardy_DK_H, "hetero")
#for (i in 1:26) print(mean(hardy_DK_H_clean$hetero[hardy_DK_H_clean$CHR==2]))

# GL_C
mean_df <- NULL;
sd_df <- NULL;
for (i in 1:26)
{ 
  tmp_mean <- mean(na.omit(hardy_GL_C$hetero[hardy_GL_C$CHR==i]))
  tmp_sd <- sd(na.omit(hardy_GL_C$hetero[hardy_GL_C$CHR==i]))
  mean_df <- rbind(mean_df, tmp_mean)
  sd_df <- rbind(sd_df, tmp_sd)
}
hardy_GL_C_mean_sd <- as.data.frame(t(rbind(mean_df[, 1], sd_df[, 1])))
row.names(hardy_GL_C_mean_sd) <- c(1:26)
colnames(hardy_GL_C_mean_sd) <- c("mean", "sd")

# GL_H
mean_df <- NULL;
sd_df <- NULL;
for (i in 1:26)
{ 
  tmp_mean <- mean(na.omit(hardy_GL_H$hetero[hardy_GL_H$CHR==i]))
  tmp_sd <- sd(na.omit(hardy_GL_H$hetero[hardy_GL_H$CHR==i]))
  mean_df <- rbind(mean_df, tmp_mean)
  sd_df <- rbind(sd_df, tmp_sd)
}
hardy_GL_H_mean_sd <- as.data.frame(t(rbind(mean_df[, 1], sd_df[, 1])))
row.names(hardy_GL_H_mean_sd) <- c(1:26)
colnames(hardy_GL_H_mean_sd) <- c("mean", "sd")

# IS_C
mean_df <- NULL;
sd_df <- NULL;
for (i in 1:26)
{ 
  tmp_mean <- mean(na.omit(hardy_IS_C$hetero[hardy_IS_C$CHR==i]))
  tmp_sd <- sd(na.omit(hardy_IS_C$hetero[hardy_IS_C$CHR==i]))
  mean_df <- rbind(mean_df, tmp_mean)
  sd_df <- rbind(sd_df, tmp_sd)
}
hardy_IS_C_mean_sd <- as.data.frame(t(rbind(mean_df[, 1], sd_df[, 1])))
row.names(hardy_IS_C_mean_sd) <- c(1:26)
colnames(hardy_IS_C_mean_sd) <- c("mean", "sd")

# IS_H
mean_df <- NULL;
sd_df <- NULL;
for (i in 1:26)
{ 
  tmp_mean <- mean(na.omit(hardy_IS_H$hetero[hardy_IS_H$CHR==i]))
  tmp_sd <- sd(na.omit(hardy_IS_H$hetero[hardy_IS_H$CHR==i]))
  mean_df <- rbind(mean_df, tmp_mean)
  sd_df <- rbind(sd_df, tmp_sd)
}
hardy_IS_H_mean_sd <- as.data.frame(t(rbind(mean_df[, 1], sd_df[, 1])))
row.names(hardy_IS_H_mean_sd) <- c(1:26)
colnames(hardy_IS_H_mean_sd) <- c("mean", "sd")

# NO_C
mean_df <- NULL;
sd_df <- NULL;
for (i in 1:26)
{ 
  tmp_mean <- mean(na.omit(hardy_NO_C$hetero[hardy_NO_C$CHR==i]))
  tmp_sd <- sd(na.omit(hardy_NO_C$hetero[hardy_NO_C$CHR==i]))
  mean_df <- rbind(mean_df, tmp_mean)
  sd_df <- rbind(sd_df, tmp_sd)
}
hardy_NO_C_mean_sd <- as.data.frame(t(rbind(mean_df[, 1], sd_df[, 1])))
row.names(hardy_NO_C_mean_sd) <- c(1:26)
colnames(hardy_NO_C_mean_sd) <- c("mean", "sd")

# NO_H
mean_df <- NULL;
sd_df <- NULL;
for (i in 1:26)
{ 
  tmp_mean <- mean(na.omit(hardy_NO_H$hetero[hardy_NO_H$CHR==i]))
  tmp_sd <- sd(na.omit(hardy_NO_H$hetero[hardy_NO_H$CHR==i]))
  mean_df <- rbind(mean_df, tmp_mean)
  sd_df <- rbind(sd_df, tmp_sd)
}
hardy_NO_H_mean_sd <- as.data.frame(t(rbind(mean_df[, 1], sd_df[, 1])))
row.names(hardy_NO_H_mean_sd) <- c(1:26)
colnames(hardy_NO_H_mean_sd) <- c("mean", "sd")

# DK_C
mean_df <- NULL;
sd_df <- NULL;
for (i in 1:26)
{ 
  tmp_mean <- mean(na.omit(hardy_DK_C$hetero[hardy_DK_C$CHR==i]))
  tmp_sd <- sd(na.omit(hardy_DK_C$hetero[hardy_DK_C$CHR==i]))
  mean_df <- rbind(mean_df, tmp_mean)
  sd_df <- rbind(sd_df, tmp_sd)
}
hardy_DK_C_mean_sd <- as.data.frame(t(rbind(mean_df[, 1], sd_df[, 1])))
row.names(hardy_DK_C_mean_sd) <- c(1:26)
colnames(hardy_DK_C_mean_sd) <- c("mean", "sd")

# DK_H
mean_df <- NULL;
sd_df <- NULL;
for (i in 1:26)
{ 
  tmp_mean <- mean(na.omit(hardy_DK_H$hetero[hardy_DK_H$CHR==i]))
  tmp_sd <- sd(na.omit(hardy_DK_H$hetero[hardy_DK_H$CHR==i]))
  mean_df <- rbind(mean_df, tmp_mean)
  sd_df <- rbind(sd_df, tmp_sd)
}
hardy_DK_H_mean_sd <- as.data.frame(t(rbind(mean_df[, 1], sd_df[, 1])))
row.names(hardy_DK_H_mean_sd) <- c(1:26)
colnames(hardy_DK_H_mean_sd) <- c("mean", "sd")

# EE_C
mean_df <- NULL;
sd_df <- NULL;
for (i in 1:26)
{ 
  tmp_mean <- mean(na.omit(hardy_EE_C$hetero[hardy_EE_C$CHR==i]))
  tmp_sd <- sd(na.omit(hardy_EE_C$hetero[hardy_EE_C$CHR==i]))
  mean_df <- rbind(mean_df, tmp_mean)
  sd_df <- rbind(sd_df, tmp_sd)
}
hardy_EE_C_mean_sd <- as.data.frame(t(rbind(mean_df[, 1], sd_df[, 1])))
row.names(hardy_EE_C_mean_sd) <- c(1:26)
colnames(hardy_EE_C_mean_sd) <- c("mean", "sd")

# TU_H
mean_df <- NULL;
sd_df <- NULL;
for (i in 1:26)
{ 
  tmp_mean <- mean(na.omit(hardy_TU_H$hetero[hardy_TU_H$CHR==i]))
  tmp_sd <- sd(na.omit(hardy_TU_H$hetero[hardy_TU_H$CHR==i]))
  mean_df <- rbind(mean_df, tmp_mean)
  sd_df <- rbind(sd_df, tmp_sd)
}
hardy_TU_H_mean_sd <- as.data.frame(t(rbind(mean_df[, 1], sd_df[, 1])))
row.names(hardy_TU_H_mean_sd) <- c(1:26)
colnames(hardy_TU_H_mean_sd) <- c("mean", "sd")


# mean and sd
#pdf("Heterozygosity_mean_sd_along_chromosome_10pops_180821.pdf", width=13, height=8)
par(mfrow=c(5,2))
par(mar=c(2.5,2.5,1.5,0.5))
plot(hardy_GL_C_mean_sd$mean, ylim=c(-0.3,0.8), ylab="", xlab="", cex.axis=1.5)
lines((hardy_GL_C_mean_sd$mean+hardy_GL_C_mean_sd$sd))
lines((hardy_GL_C_mean_sd$mean-hardy_GL_C_mean_sd$sd))
title("A", adj=0.02, line = 0.5)
plot(hardy_GL_H_mean_sd$mean, ylim=c(-0.3,0.8), ylab="", xlab="", cex.axis=1.5)
lines((hardy_GL_H_mean_sd$mean+hardy_GL_H_mean_sd$sd))
lines((hardy_GL_H_mean_sd$mean-hardy_GL_H_mean_sd$sd))
title("B", adj=0.02, line = 0.5)
plot(hardy_IS_C_mean_sd$mean, ylim=c(-0.3,0.8), ylab="", xlab="", cex.axis=1.5)
lines((hardy_IS_C_mean_sd$mean+hardy_IS_C_mean_sd$sd))
lines((hardy_IS_C_mean_sd$mean-hardy_IS_C_mean_sd$sd))
title("C", adj=0.02, line = 0.5)
plot(hardy_IS_H_mean_sd$mean, ylim=c(-0.3,0.8), ylab="", xlab="", cex.axis=1.5)
lines((hardy_IS_H_mean_sd$mean+hardy_IS_H_mean_sd$sd))
lines((hardy_IS_H_mean_sd$mean-hardy_IS_H_mean_sd$sd))
title("D", adj=0.02, line = 0.5)
plot(hardy_NO_C_mean_sd$mean, ylim=c(-0.3,0.8), ylab="", xlab="", cex.axis=1.5)
lines((hardy_NO_C_mean_sd$mean+hardy_NO_C_mean_sd$sd))
lines((hardy_NO_C_mean_sd$mean-hardy_NO_C_mean_sd$sd))
title("E", adj=0.02, line = 0.5)
plot(hardy_NO_H_mean_sd$mean, ylim=c(-0.3,0.8), ylab="", xlab="", cex.axis=1.5)
lines((hardy_NO_H_mean_sd$mean+hardy_NO_H_mean_sd$sd))
lines((hardy_NO_H_mean_sd$mean-hardy_NO_H_mean_sd$sd))
title("F", adj=0.02, line = 0.5)
plot(hardy_DK_C_mean_sd$mean, ylim=c(-0.3,0.8), ylab="", xlab="", cex.axis=1.5)
lines((hardy_DK_C_mean_sd$mean+hardy_DK_C_mean_sd$sd))
lines((hardy_DK_C_mean_sd$mean-hardy_DK_C_mean_sd$sd))
title("G", adj=0.02, line = 0.5)
plot(hardy_DK_H_mean_sd$mean, ylim=c(-0.3,0.8), ylab="", xlab="", cex.axis=1.5)
lines((hardy_DK_H_mean_sd$mean+hardy_DK_H_mean_sd$sd))
lines((hardy_DK_H_mean_sd$mean-hardy_DK_H_mean_sd$sd))
title("H", adj=0.02, line = 0.5)
plot(hardy_EE_C_mean_sd$mean, ylim=c(-0.3,0.8), ylab="", xlab="", cex.axis=1.5)
lines((hardy_EE_C_mean_sd$mean+hardy_EE_C_mean_sd$sd))
lines((hardy_EE_C_mean_sd$mean-hardy_EE_C_mean_sd$sd))
title("I", adj=0.02, line = 0.5)
plot(hardy_TU_H_mean_sd$mean, ylim=c(-0.3,0.8), ylab="", xlab="", cex.axis=1.5)
lines((hardy_TU_H_mean_sd$mean+hardy_TU_H_mean_sd$sd))
lines((hardy_TU_H_mean_sd$mean-hardy_TU_H_mean_sd$sd))
title("J", adj=0.02, line = 0.5)
#dev.off()

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




## Tree from Sina
#vcf_path <- "all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newsamplename.recode.vcf"

#sometimes spits the dummy, then execute
showfile.gds(closeall=TRUE)

#turn vcf into gds
#snpgdsVCF2GDS(vcf_path,"all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newsamplename.recode.gds",method ="biallelic.only")

#format data correctly to create a Identity by State matrix
genofile<-snpgdsOpen("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newsamplename.recode.gds", readonly = T)
pop_code <- scan("popcode_for_tree_040621.txt", what=character())
#samp_code <- scan("popcode_for_tree_040621_2.txt", what=character())
ibs.hc<-snpgdsHCluster(snpgdsIBS(genofile,num.thread=2, autosome.only=FALSE))

#test<-snpgdsIBS(genofile,num.thread=2, autosome.only=FALSE)
#write.table(test$ibs, "all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newsamplename_IBS_dist_MAT.txt")
#write.table(ibs.hc$dist, "all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_maf0.05_newsamplename_IBS_dist_clustMAT.txt")

#take clustering results from before and turn the numerical values into a dendrogram
rv <- snpgdsCutTree(ibs.hc)
#rv$sample.id=c(read.table("popcode_for_tree_040621.txt"))
#samp.annot<- data.frame(pop.group = c("IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "DK_H", "DK_H", "DK_H", "DK_H", "DK_H", "DK_C", "DK_C", "DK_C", "DK_C", "DK_C", "DK_C", "DK_C", "DK_C", "DK_C", "DK_C", "DK_C", "EE_C", "EE_C", "EE_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_H", "GL_H", "GL_H", "GL_H", "GL_H", "GL_H", "GL_H", "GL_H", "IS_H", "IS_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "TU_H"))
#add.gdsn(genofile, "sample.annot", samp.annot)
#pop_code <- read.gdsn(index.gdsn(genofile, path="popcode_for_tree_040621.txt"))
rv2 <- snpgdsCutTree(ibs.hc, samp.group = as.factor(pop_code), label.Z=T, n.perm=1000)

#plot dendrogram
#pdf("dendogram_all10_pops_040621.pdf", width = 16, height = 8)
par(mar=c(4,2,0.2,0))
plot(rv2$dendrogram, main="", col=c(brewer.pal(10,"Paired")))
#dev.off()
plot(rv$dendrogram, main="Dendogram from IBS", col=c(brewer.pal(10,"Paired")))

## Tree from Sina ONLY IS
#vcf_path_IS <- "all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss1_HetHom_minMQ30_onlyISall.recode.vcf"

#sometimes spits the dummy, then execute
showfile.gds(closeall=TRUE)

#turn vcf into gds
#snpgdsVCF2GDS(vcf_path_IS,"all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss1_HetHom_minMQ30_onlyISall.recode.gds",method ="biallelic.only")

#format data correctly to create a Identity by State matrix
genofile_IS<-snpgdsOpen("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss1_HetHom_minMQ30_onlyISall.recode.gds", readonly = T)
#pop_code_IS <- scan("popcode_for_tree_040621.txt", what=character())
#samp_code_IS <- scan("popcode_for_tree_040621_2.txt", what=character())
set.seed(1000)
ibs.hc_IS<-snpgdsHCluster(snpgdsIBS(genofile_IS,num.thread=2, autosome.only=FALSE))

#take clustering results from before and turn the numerical values into a dendrogram
rv_IS <- snpgdsCutTree(ibs.hc_IS)
#rv$sample.id=c(read.table("popcode_for_tree_040621.txt"))
#samp.annot<- data.frame(pop.group = c("IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "IS_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "NO_C", "DK_H", "DK_H", "DK_H", "DK_H", "DK_H", "DK_C", "DK_C", "DK_C", "DK_C", "DK_C", "DK_C", "DK_C", "DK_C", "DK_C", "DK_C", "DK_C", "EE_C", "EE_C", "EE_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_C", "GL_H", "GL_H", "GL_H", "GL_H", "GL_H", "GL_H", "GL_H", "GL_H", "IS_H", "IS_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "NO_H", "TU_H"))
#add.gdsn(genofile, "sample.annot", samp.annot)
#pop_code <- read.gdsn(index.gdsn(genofile, path="popcode_for_tree_040621.txt"))
#rv2 <- snpgdsCutTree(ibs.hc_IS, samp.group = as.factor(pop_code), label.Z=T, n.perm=1000)

#plot dendrogram
#pdf("dendogram_all10_pops_040621.pdf", width = 16, height = 8)
par(mar=c(4,2,0.2,0))
plot(rv2$dendrogram, main="", col=c(brewer.pal(10,"Paired")))
#dev.off()
plot(rv_IS$dendrogram, main="Dendrogram from IBS", col=c(brewer.pal(10,"Paired")))

showfile.gds(closeall=TRUE)
#vcf_path_GL <- "all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss1_HetHom_minMQ30_onlyGLall.recode.vcf"
#snpgdsVCF2GDS(vcf_path_GL,"all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss1_HetHom_minMQ30_onlyGLall.recode.gds",method ="biallelic.only")
genofile_GL<-snpgdsOpen("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss1_HetHom_minMQ30_onlyGLall.recode.gds", readonly = T)
set.seed(1000)
ibs.hc_GL<-snpgdsHCluster(snpgdsIBS(genofile_GL,num.thread=2, autosome.only=FALSE))

rv_GL <- snpgdsCutTree(ibs.hc_GL)

plot(rv_GL$dendrogram, main="Dendrogram from IBS", col=c(brewer.pal(10,"Paired")))


#### het vs depth
mean_depth_true<-read.table("Truemeandepth_onlysitespresentinindividual.txt", header = T)
missingness<-read.table("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.imiss", header = T)
het_data_withdepth<-het_data
het_data_withdepth$mean_depth<-mean_depth_true$MEAN_DEPTH
het_data_withdepth$F_MISS<-missingness$F_MISS

png("HeteroVsMeandepth_130921.png", width = 1000, height = 1000)
par(mar=c(5,5,0.2,0.2))
plot(het_data_withdepth$het_frac, het_data_withdepth$mean_depth, xlab=c("Heterozygosity"), ylab=c("Mean depth"), pch=20, cex=3.5, cex.lab=2)
dev.off()

png("HeteroVsMissingness_130921.png", width = 1000, height = 1000)
par(mar=c(5,5,0.2,0.2))
plot(het_data_withdepth$het_frac, het_data_withdepth$F_MISS, xlab=c("Heterozygosity"), ylab=c("Missingness"), pch=20, cex=3.5, cex.lab=2)
dev.off()

cor(het_data_withdepth$het_frac, het_data_withdepth$mean_depth)
#-0.08692925
cor(het_data_withdepth$F_MISS, het_data_withdepth$het_frac)
#0.3317026

### ROH vs depth
ROH_92_relaxed_310321_withmissing$meandepth<-het_data_withdepth$mean_depth
png("LenghtOfROHVsMeandepth_130921.png", width = 1000, height = 1000)
par(mar=c(5,5,0.2,0.2))
plot(ROH_92_relaxed_310321_withmissing$KB/1000, ROH_92_relaxed_310321_withmissing$meandepth, xlab=c("Length of ROH (Mb)"), ylab=c("Mean depth"), pch=20, cex=3.5, cex.lab=2)
dev.off()

cor(ROH_92_relaxed_310321_withmissing$KB, ROH_92_relaxed_310321_withmissing$meandepth)
#-0.2314151
cor(ROH_92_relaxed_310321_withmissing$KB, ROH_92_relaxed_310321_withmissing$MISS)
#0.5740728

### Wrights inbreeding - Snæbjörn formula
hardy_GL_C_basic<-read.table("hardy_het_per_site/hardy_GL_C_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe", header = T)
hardy_GL_H_basic<-read.table("hardy_het_per_site/hardy_GL_H_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe", header = T)
hardy_IS_C_basic<-read.table("hardy_het_per_site/hardy_IS_C_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe", header = T)
hardy_IS_H_basic<-read.table("hardy_het_per_site/hardy_IS_H_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe", header = T)
hardy_NO_C_basic<-read.table("hardy_het_per_site/hardy_NO_C_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe", header = T)
hardy_NO_H_basic<-read.table("hardy_het_per_site/hardy_NO_H_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe", header = T)
hardy_DK_C_basic<-read.table("hardy_het_per_site/hardy_DK_C_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe", header = T)
hardy_DK_H_basic<-read.table("hardy_het_per_site/hardy_DK_H_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe", header = T)
hardy_EE_C_basic<-read.table("hardy_het_per_site/hardy_EE_C_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe", header = T)
hardy_TU_H_basic<-read.table("hardy_het_per_site/hardy_TU_H_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe", header = T)

hardy_all_basic<-read.table("hardy_het_per_site/hardy_all_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.hwe", header = T)
  
Fis=function(obs,exp,P)
{
  #calculates Fis from HW per loci in vcftools
  counts =as.numeric(unlist(strsplit(obs,"/")))
  n=sum(counts)
  HO=(counts[2]/n)*2*n/(2*n-1)
  HE=(as.numeric(unlist(strsplit(exp,"/"))[2])/n)*2*n/(2*n-1)
  P=P
  Fis=1-HO/HE
  return(c(Fis,HO,HE,P))
}


Fis.GLC=matrix(NA,nrow=nrow(hardy_GL_C_basic),ncol=4)
for(i in 1:nrow(Fis.GLC)) Fis.GLC[i,]=  
  Fis(hardy_GL_C_basic[i,3],hardy_GL_C_basic[i,4],hardy_GL_C_basic[i,6])
1-mean(na.omit(Fis.GLC[!is.na(Fis.GLC[,3]),2]))/mean(na.omit(Fis.GLC[!is.na(Fis.GLC[,3]),3]))
mean_Fis.GLC<-apply(Fis.GLC, 2, mean, na.rm=T)
mean_Fis.GLC[1]<-1-mean(na.omit(Fis.GLC[!is.na(Fis.GLC[,3]),2]))/mean(na.omit(Fis.GLC[!is.na(Fis.GLC[,3]),3]))
mean_Fis.GLC_var<-apply(Fis.GLC, 2, var, na.rm=T)
#mean_Fis.GLC_sd[1]<-1-mean(na.omit(Fis.GLC[!is.na(Fis.GLC[,3]),2]))/mean(na.omit(Fis.GLC[!is.na(Fis.GLC[,3]),3]))
#[1] -0.08885921

Fis.GLH=matrix(NA,nrow=nrow(hardy_GL_H_basic),ncol=4)
for(i in 1:nrow(Fis.GLH)) Fis.GLH[i,]=  
  Fis(hardy_GL_H_basic[i,3],hardy_GL_H_basic[i,4],hardy_GL_H_basic[i,6])
1-mean(na.omit(Fis.GLH[!is.na(Fis.GLH[,3]),2]))/mean(na.omit(Fis.GLH[!is.na(Fis.GLH[,3]),3]))
mean_Fis.GLH<-apply(Fis.GLH, 2, mean, na.rm=T)
mean_Fis.GLH[1]<-1-mean(na.omit(Fis.GLH[!is.na(Fis.GLH[,3]),2]))/mean(na.omit(Fis.GLH[!is.na(Fis.GLH[,3]),3]))
mean_Fis.GLH_var<-apply(Fis.GLH, 2, var, na.rm=T)
#[1] -0.160006

Fis.ISC=matrix(NA,nrow=778735,ncol=4)
Fis.ISC=matrix(NA,nrow(hardy_IS_C_basic),ncol=4)
for(i in 1:nrow(Fis.ISC)) Fis.ISC[i,]=
  Fis(hardy_IS_C_basic[i,3],hardy_IS_C_basic[i,4],hardy_IS_C_basic[i,6])
1-mean(na.omit(Fis.ISC[!is.na(Fis.ISC[,3]),2]))/mean(na.omit(Fis.ISC[!is.na(Fis.ISC[,3]),3]))
mean_Fis.ISC<-apply(Fis.ISC, 2, mean, na.rm=T)
mean_Fis.ISC[1]<-1-mean(na.omit(Fis.ISC[!is.na(Fis.ISC[,3]),2]))/mean(na.omit(Fis.ISC[!is.na(Fis.ISC[,3]),3]))
mean_Fis.ISC_var<-apply(Fis.ISC, 2, var, na.rm=T)
#[1] -0.0488635

Fis.ISH=matrix(NA,nrow(hardy_IS_H_basic),ncol=4)
for(i in 1:nrow(Fis.ISH)) Fis.ISH[i,]=
  Fis(hardy_IS_H_basic[i,3],hardy_IS_H_basic[i,4],hardy_IS_H_basic[i,6])
1-mean(na.omit(Fis.ISH[!is.na(Fis.ISH[,3]),2]))/mean(na.omit(Fis.ISH[!is.na(Fis.ISH[,3]),3]))
mean_Fis.ISH<-apply(Fis.ISH, 2, mean, na.rm=T)
mean_Fis.ISH[1]<-1-mean(na.omit(Fis.ISH[!is.na(Fis.ISH[,3]),2]))/mean(na.omit(Fis.ISH[!is.na(Fis.ISH[,3]),3]))
mean_Fis.ISH_var<-apply(Fis.ISH, 2, var, na.rm=T)
#[1] 0.3165549

Fis.NOC=matrix(NA,nrow(hardy_NO_C_basic),ncol=4)
for(i in 1:nrow(Fis.NOC)) Fis.NOC[i,]=
  Fis(hardy_NO_C_basic[i,3],hardy_NO_C_basic[i,4],hardy_NO_C_basic[i,6])
1-mean(na.omit(Fis.NOC[!is.na(Fis.NOC[,3]),2]))/mean(na.omit(Fis.NOC[!is.na(Fis.NOC[,3]),3]))
mean_Fis.NOC<-apply(Fis.NOC, 2, mean, na.rm=T)
mean_Fis.NOC[1]<-1-mean(na.omit(Fis.NOC[!is.na(Fis.NOC[,3]),2]))/mean(na.omit(Fis.NOC[!is.na(Fis.NOC[,3]),3]))
mean_Fis.NOC_var<-apply(Fis.NOC, 2, var, na.rm=T)
#[1] -0.07270307

Fis.NOH=matrix(NA,nrow(hardy_NO_H_basic),ncol=4)
for(i in 1:nrow(Fis.NOH)) Fis.NOH[i,]=
  Fis(hardy_NO_H_basic[i,3],hardy_NO_H_basic[i,4],hardy_NO_H_basic[i,6])
1-mean(na.omit(Fis.NOH[!is.na(Fis.NOH[,3]),2]))/mean(na.omit(Fis.NOH[!is.na(Fis.NOH[,3]),3]))
mean_Fis.NOH<-apply(Fis.NOH, 2, mean, na.rm=T)
mean_Fis.NOH[1]<-1-mean(na.omit(Fis.NOH[!is.na(Fis.NOH[,3]),2]))/mean(na.omit(Fis.NOH[!is.na(Fis.NOH[,3]),3]))
mean_Fis.NOH_var<-apply(Fis.NOH, 2, var, na.rm=T)
#[1] -0.1091536

Fis.DKC=matrix(NA,nrow(hardy_DK_C_basic),ncol=4)
for(i in 1:nrow(Fis.DKC)) Fis.DKC[i,]=
  Fis(hardy_DK_C_basic[i,3],hardy_DK_C_basic[i,4],hardy_DK_C_basic[i,6])
1-mean(na.omit(Fis.DKC[!is.na(Fis.DKC[,3]),2]))/mean(na.omit(Fis.DKC[!is.na(Fis.DKC[,3]),3]))
mean_Fis.DKC<-apply(Fis.DKC, 2, mean, na.rm=T)
mean_Fis.DKC[1]<-1-mean(na.omit(Fis.DKC[!is.na(Fis.DKC[,3]),2]))/mean(na.omit(Fis.DKC[!is.na(Fis.DKC[,3]),3]))
mean_Fis.DKC_var<-apply(Fis.DKC, 2, var, na.rm=T)
#[1] -0.0431765

Fis.DKH=matrix(NA,nrow(hardy_DK_H_basic),ncol=4)
for(i in 1:nrow(Fis.DKH)) Fis.DKH[i,]=
  Fis(hardy_DK_H_basic[i,3],hardy_DK_H_basic[i,4],hardy_DK_H_basic[i,6])
1-mean(na.omit(Fis.DKH[!is.na(Fis.DKH[,3]),2]))/mean(na.omit(Fis.DKH[!is.na(Fis.DKH[,3]),3]))
mean_Fis.DKH<-apply(Fis.DKH, 2, mean, na.rm=T)
mean_Fis.DKH[1]<-1-mean(na.omit(Fis.DKH[!is.na(Fis.DKH[,3]),2]))/mean(na.omit(Fis.DKH[!is.na(Fis.DKH[,3]),3]))
mean_Fis.DKH_var<-apply(Fis.DKH, 2, var, na.rm=T)
#[1] -0.2711167

Fis.EEC=matrix(NA,nrow(hardy_EE_C_basic),ncol=4)
for(i in 1:nrow(Fis.EEC)) Fis.EEC[i,]=
  Fis(hardy_EE_C_basic[i,3],hardy_EE_C_basic[i,4],hardy_EE_C_basic[i,6])
1-mean(na.omit(Fis.EEC[!is.na(Fis.EEC[,3]),2]))/mean(na.omit(Fis.EEC[!is.na(Fis.EEC[,3]),3]))
mean_Fis.EEC<-apply(Fis.EEC, 2, mean, na.rm=T)
mean_Fis.EEC[1]<-1-mean(na.omit(Fis.EEC[!is.na(Fis.EEC[,3]),2]))/mean(na.omit(Fis.EEC[!is.na(Fis.EEC[,3]),3]))
mean_Fis.EEC_var<-apply(Fis.EEC, 2, var, na.rm=T)
#[1] -0.2254058

Fis.TUH=matrix(NA,nrow(hardy_TU_H_basic),ncol=4)
for(i in 1:nrow(Fis.TUH)) Fis.TUH[i,]=
  Fis(hardy_TU_H_basic[i,3],hardy_TU_H_basic[i,4],hardy_TU_H_basic[i,6])
1-mean(na.omit(Fis.TUH[!is.na(Fis.TUH[,3]),2]))/mean(na.omit(Fis.TUH[!is.na(Fis.TUH[,3]),3]))
mean_Fis.TUH<-apply(Fis.TUH, 2, mean, na.rm=T)
mean_Fis.TUH[1]<-1-mean(na.omit(Fis.TUH[!is.na(Fis.TUH[,3]),2]))/mean(na.omit(Fis.TUH[!is.na(Fis.TUH[,3]),3]))
mean_Fis.TUH_var<-apply(Fis.TUH, 2, var, na.rm=T)
#[1] -0.8554759

Fis.all=matrix(NA,nrow(hardy_all_basic),ncol=4)
for(i in 1:nrow(Fis.all)) Fis.all[i,]=
  Fis(hardy_all_basic[i,3],hardy_all_basic[i,4],hardy_all_basic[i,6])
1-mean(na.omit(Fis.all[!is.na(Fis.all[,3]),2]))/mean(na.omit(Fis.all[!is.na(Fis.all[,3]),3]))
mean_Fis.all<-apply(Fis.all, 2, mean, na.rm=T)
mean_Fis.all[1]<-1-mean(na.omit(Fis.all[!is.na(Fis.all[,3]),2]))/mean(na.omit(Fis.all[!is.na(Fis.all[,3]),3]))
mean_Fis.all_var<-apply(Fis.all, 2, var, na.rm=T)  

Fis.all<-as.data.frame(round(rbind(mean_Fis.GLC, mean_Fis.GLH, mean_Fis.ISC, mean_Fis.ISH, mean_Fis.NOC, mean_Fis.NOH, 
            mean_Fis.DKC, mean_Fis.DKH, mean_Fis.EEC, mean_Fis.TUH ), 3))

colnames(Fis.all)<-c("Fis", "Het_obs", "Het_exp", "P")
#write.table(Fis.all, "Fis_Het_wrights_041021.txt", quote = F)

Fis.all_var<-as.data.frame(round(rbind(mean_Fis.GLC_var, mean_Fis.GLH_var, mean_Fis.ISC_var, mean_Fis.ISH_var, mean_Fis.NOC_var, mean_Fis.NOH_var, 
                                   mean_Fis.DKC_var, mean_Fis.DKH_var, mean_Fis.EEC_var, mean_Fis.TUH_var ), 3))

colnames(Fis.all_var)<-c("Fis", "Het_obs", "Het_exp", "P")

### shitty het calculated per pop
het_data_new<-read.table("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_all_individualpopcalc.het.txt", h = T)
#colnames(het_data_new)[1]<-c("ORDER")
head(het_data_new)
het_data_new$O.HET<-het_data_new$N_SITES-het_data_new$O.HOM.
het_data_new$E.HET<-het_data_new$N_SITES-het_data_new$E.HOM.
het_data_new$het_frac<-het_data_new$O.HET/het_data_new$N_SITES
het_data_new$CO_TI<-factor(het_data_new$CO_TI, levels = c("GL_C","GL_H","IS_C","IS_H","NO_C","NO_H","DK_C","DK_H","EE_C","TU_H"))
het_data_new[order(het_data_new$O.HET),]

het_obs_tmp_new<-het_data_new[,c(1,2,3,4,7,9)]
het_data_new_obs_exp_for_ggplot<-het_obs_tmp_new[order(het_obs_tmp_new$ORDER),]
het_data_new_obs_exp_for_ggplot$obs_or_exp<-rep("obs", 92)
colnames(het_data_new_obs_exp_for_ggplot)[6]<-c("HET")

het_exp_tmp_new<-het_data_new[,c(1,2,3,4,7,10)]
het_exp_tmp_new2<-het_exp_tmp_new[order(het_exp_tmp_new$ORDER),]
het_exp_tmp_new2$obs_or_exp<-rep("exp", 92)
colnames(het_exp_tmp_new2)[6]<-c("HET")

het_inb_tmp_new<-het_data_new[,c(1,2,3,4,7,8)]
het_inb_tmp_new2<-het_inb_tmp_new[order(het_inb_tmp_new$ORDER),]
het_inb_tmp_new2$obs_or_exp<-rep("inb", 92)
colnames(het_inb_tmp_new2)[6]<-c("HET")

het_data_new_obs_exp_for_ggplot_true<-rbind(het_data_new_obs_exp_for_ggplot, het_exp_tmp_new2)

ggplot_het_bothexpandobsandinb_300921<-
  ggplot(het_data_new_obs_exp_for_ggplot_true, aes(x=CO_TI, y=HET/N_SITES, color = obs_or_exp)) +
  geom_boxplot() +
  geom_boxplot(data = het_inb_tmp_new2, aes(x=CO_TI, y=HET, color = obs_or_exp), fill="#636363") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        title = element_text(size = 16), legend.text = element_text(size = 14)) +
  #annotate(geom= "text", x=seq_len(unique(het_data_new_obs_exp_for_ggplot_true$CO_TI)), y=10, label=het_data_new_obs_exp_for_ggplot_true$obs_or_exp) +
  scale_color_manual(values = c( "#aaaaaa", "#000000", "#000000")) + 
  labs(x = "Population_time", 
       y = "Heterozygosity & Inbreeding coefficient",
       color="") + 
  theme(legend.position = "none")

#ggsave(plot = ggplot_het_bothexpandobsandinb_300921, filename = "het_plot_bothhetandobsandinb_300921_perpopF.png", width = 30, height = 20, units = "cm", device = "png")

#roh92_<-all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_newchr26_relaxed_out_310321_1.hom.inv

Froh<-read.table("Froh_300921.txt", header = T)
missing<-read.table("all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.imiss", header=T)
head(Froh)
Froh$Fh<-het_data$F
Froh$imiss_F<-missing$F_MISS
plot(Froh$Froh, Froh$Fh)
cor(Froh$Froh[Froh$imiss_F<0.5], Froh$Fh[Froh$imiss_F<0.5])
#[1] 0.9591388
sd(Froh$Froh[Froh$imiss_F<0.5]-Froh$Fh[Froh$imiss_F<0.5])
max(Froh$Froh[Froh$imiss_F<0.5]-Froh$Fh[Froh$imiss_F<0.5])
min(Froh$Froh[Froh$imiss_F<0.5]-Froh$Fh[Froh$imiss_F<0.5])
#[1] 0.08057715
#png("Rplot_inbreeding_Fh_Froh_78indfrom_92ind_miss75_061021.png", height = 800, width = 800)
par(mar=c(4.2,4.5,0.2,0.2))
plot(Froh$Froh[Froh$imiss_F<0.5], Froh$Fh[Froh$imiss_F<0.5], cex=1.5, pch=19, xlab=expression('F'[ROH]), ylab= expression('F'[H]), cex.axis = 1.5, cex.lab = 1.5)
#dev.off()


##LD testing
# LD_no_200k<-read.table("LD_test_200K_NO_modern_from92_noBB38_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_nonan.txt.geno.ld", header=T, fill=T)
# LD_no_200k$dist<-LD_no_200k$POS2-LD_no_200k$POS1
# LD_no_200k$dist_kb<-LD_no_200k$dist%/%1000
# head(LD_no_200k)
# LD_no_200k_NoNaN<-na.omit(LD_no_200k)
# length(LD_no_200k_NoNaN$CHR)
# LD_no_200k_NoNaN_tapply_mean<-as.data.frame(tapply(LD_no_200k_NoNaN$R.2, LD_no_200k_NoNaN$dist_kb, mean))
# write.table(LD_no_200k_NoNaN_tapply_mean, "LD_NO_200k_NoNaN_tapply_mean_LD_test_200k_NO_modern_from92_noBB38_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.txt")
LD_no_200k_NoNaN_tapply_mean<-read.table("LD_NO_200k_NoNaN_tapply_mean_LD_test_200k_NO_modern_from92_noBB38_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.txt", header =T)
colnames(LD_no_200k_NoNaN_tapply_mean)<-c("mean_LD")
plot(rownames(LD_no_200k_NoNaN_tapply_mean), LD_no_200k_NoNaN_tapply_mean$mean_LD)
head(LD_no_200k_NoNaN_tapply_mean)
#0.3919585/2
#View(LD_no_200k_NoNaN_tapply_mean)
#34 0.1952584

# LD_is_200k<-read.table("LD_test_200K_IS_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_nonan.txt.geno.ld", header=T, fill=T)
# LD_is_200k$dist<-LD_is_200k$POS2-LD_is_200k$POS1
# LD_is_200k$dist_kb<-LD_is_200k$dist%/%1000
# head(LD_is_200k)
# LD_is_200k_NoNaN<-na.omit(LD_is_200k)
# length(LD_is_200k_NoNaN$CHR)
# LD_is_200k_NoNaN_tapply_mean<-as.data.frame(tapply(LD_is_200k_NoNaN$R.2, LD_is_200k_NoNaN$dist_kb, mean))
# write.table(LD_is_200k_NoNaN_tapply_mean, "LD_IS_200k_NoNaN_tapply_mean_LD_test_200k_IS_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.txt")
LD_is_200k_NoNaN_tapply_mean<-read.table("LD_IS_200k_NoNaN_tapply_mean_LD_test_200k_IS_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30.txt", header =T)
colnames(LD_is_200k_NoNaN_tapply_mean)<-c("mean_LD")
plot(rownames(LD_is_200k_NoNaN_tapply_mean), LD_is_200k_NoNaN_tapply_mean$mean_LD)
head(LD_is_200k_NoNaN_tapply_mean)
#0.6654564/2
#View(LD_is_200k_NoNaN_tapply_mean)
# 134 0.3316085

# LD_gl_200k<-read.table("LD_test_200K_GL_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_nonan.txt.geno.ld", header=T, fill=T)
# LD_gl_200k$dist<-LD_gl_200k$POS2-LD_gl_200k$POS1
# LD_gl_200k$dist_kb<-LD_gl_200k$dist%/%1000
# head(LD_gl_200k)
# LD_gl_200k_NoNaN<-na.omit(LD_gl_200k)
# length(LD_gl_200k_NoNaN$CHR)
# LD_gl_200k_NoNaN_tapply_mean<-as.data.frame(tapply(LD_gl_200k_NoNaN$R.2, LD_gl_200k_NoNaN$dist_kb, mean))
# write.table(LD_gl_200k_NoNaN_tapply_mean, "LD_GL_200k_NoNaN_tapply_mean_LD_test_200k_GL_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_mgls0.75_HetHom_minMQ30.txt")
LD_gl_200k_NoNaN_tapply_mean<-read.table("LD_GL_200k_NoNaN_tapply_mean_LD_test_200k_GL_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_mgls0.75_HetHom_minMQ30.txt", header =T)
colnames(LD_gl_200k_NoNaN_tapply_mean)<-c("mean_LD")
plot(rownames(LD_gl_200k_NoNaN_tapply_mean), LD_gl_200k_NoNaN_tapply_mean$mean_LD)
head(LD_gl_200k_NoNaN_tapply_mean)
#0.6719751/2
#View(LD_gl_200k_NoNaN_tapply_mean)
# 56 0.3351957

# LD_dk_200k<-read.table("LD_test_200K_DK_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_miss0.75_HetHom_minMQ30_nonan.txt.geno.ld", header=T, fill=T)
# LD_dk_200k$dist<-LD_dk_200k$POS2-LD_dk_200k$POS1
# LD_dk_200k$dist_kb<-LD_dk_200k$dist%/%1000
# head(LD_dk_200k)
# LD_dk_200k_NoNaN<-na.omit(LD_dk_200k)
# length(LD_dk_200k_NoNaN$CHR)
# LD_dk_200k_NoNaN_tapply_mean<-as.data.frame(tapply(LD_dk_200k_NoNaN$R.2, LD_dk_200k_NoNaN$dist_kb, mean))
# write.table(LD_dk_200k_NoNaN_tapply_mean, "LD_DK_200k_NoNaN_tapply_mean_LD_test_200k_DK_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_mdks0.75_HetHom_minMQ30.txt")
LD_dk_200k_NoNaN_tapply_mean<-read.table("LD_DK_200k_NoNaN_tapply_mean_LD_test_200k_DK_modern_from92_all_results_mac1_92ind_Q1000_GQ20_DP8_autosomes_mdks0.75_HetHom_minMQ30.txt", header =T)
colnames(LD_dk_200k_NoNaN_tapply_mean)<-c("mean_LD")
plot(rownames(LD_dk_200k_NoNaN_tapply_mean), LD_dk_200k_NoNaN_tapply_mean$mean_LD)
head(LD_dk_200k_NoNaN_tapply_mean)
# 0.4063505/2
# View(LD_dk_200k_NoNaN_tapply_mean)
# 46 0.2023231

LD_is_200k_NoNaN_tapply_mean[ "kb_dist" ] <- rownames(LD_is_200k_NoNaN_tapply_mean)
LD_all_200k_NoNaN_tapply_mean<-as.data.frame(LD_is_200k_NoNaN_tapply_mean$kb_dist)
LD_all_200k_NoNaN_tapply_mean$is_mean_LD<-LD_is_200k_NoNaN_tapply_mean$mean_LD
LD_all_200k_NoNaN_tapply_mean$no_mean_LD<-LD_no_200k_NoNaN_tapply_mean$mean_LD
LD_all_200k_NoNaN_tapply_mean$gl_mean_LD<-LD_gl_200k_NoNaN_tapply_mean$mean_LD
LD_all_200k_NoNaN_tapply_mean$dk_mean_LD<-LD_dk_200k_NoNaN_tapply_mean$mean_LD
colnames(LD_all_200k_NoNaN_tapply_mean)[1]<-c("kb_dist")

ggplot_LD_mean_allpops<-
ggplot(LD_all_200k_NoNaN_tapply_mean, aes(x=order(kb_dist))) +#, y=mean_LD)) + #, color = aes(red))) +
  geom_line(aes(y=is_mean_LD[order(kb_dist)]), colour = "red") +
  geom_point(aes(x=134, y=0.331), colour="red", size=2) + 
  geom_line(aes(y=gl_mean_LD[order(kb_dist)]), colour = "green") +
  geom_point(aes(x=56, y=0.335), colour="green", size=2) + 
  geom_line(aes(y=no_mean_LD[order(kb_dist)]), colour = "orange") +
  geom_point(aes(x=34, y=0.195), colour="orange", size=2) +
  geom_line(aes(y=dk_mean_LD[order(kb_dist)]), colour = "blue") +  
  geom_point(aes(x=46, y=0.202), colour="blue", size=2) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        title = element_text(size = 16), legend.text = element_text(size = 14)) +
  labs(x = "kb distance", 
       y = "Mean R2",
       color="") + 
  theme(legend.position = "none")

#ggsave(plot = ggplot_LD_mean_allpops, filename = "ggplot_LD_mean_allpops_130421.png", width = 30, height = 20, units = "cm", device = "png")
#ggsave(plot = ggplot_LD_mean_allpops, filename = "ggplot_LD_mean_allpops_130421.pdf", width = 30, height = 20, units = "cm", device = "pdf")
