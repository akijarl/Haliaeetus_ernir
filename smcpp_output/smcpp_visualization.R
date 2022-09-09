setwd("/Users/akijarl/Documents/UI/rannsoknir/Haliaeetus/Haliaeetus_ernir/smcpp_output/")
require(ggplot2)

GC1<-read.csv("Greenland_con1.csv")
GC2<-read.csv("Greenland_con2.csv")
GC3<-read.csv("Greenland_con3.csv")
GC4<-read.csv("Greenland_con4.csv")
GC5<-read.csv("Greenland_con5.csv")
GC6<-read.csv("Greenland_con6.csv")
GC7<-read.csv("Greenland_con7.csv")
GC8<-read.csv("Greenland_con8.csv")
GC9<-read.csv("Greenland_con9.csv")
GC10<-read.csv("Greenland_con10.csv")

GC1$replicate<-"A"
GC2$replicate<-"B"
GC3$replicate<-"C"
GC4$replicate<-"D"
GC5$replicate<-"E"
GC6$replicate<-"F"
GC7$replicate<-"G"
GC8$replicate<-"H"
GC9$replicate<-"I"
GC10$replicate<-"J"

GC<-rbind(GC1,GC2,GC3,GC4,GC5,GC6,GC7,GC8,GC9,GC10)

IC1<-read.csv("Iceland_con1.csv")
IC2<-read.csv("Iceland_con2.csv")
IC3<-read.csv("Iceland_con3.csv")
IC4<-read.csv("Iceland_con4.csv")
IC5<-read.csv("Iceland_con5.csv")
IC6<-read.csv("Iceland_con6.csv")
IC7<-read.csv("Iceland_con7.csv")
IC8<-read.csv("Iceland_con8.csv")
IC9<-read.csv("Iceland_con9.csv")
IC10<-read.csv("Iceland_con10.csv")

IC1$replicate<-"A"
IC2$replicate<-"B"
IC3$replicate<-"C"
IC4$replicate<-"D"
IC5$replicate<-"E"
IC6$replicate<-"F"
IC7$replicate<-"G"
IC8$replicate<-"H"
IC9$replicate<-"I"
IC10$replicate<-"J"

IC<-rbind(IC1,IC2,IC3,IC4,IC5,IC6,IC7,IC8,IC9,IC10)

ggplot(IC[IC$x>100,])+
  geom_line(aes(x,log10(y),colour=replicate))+
  coord_cartesian(xlim=c(0,100000),ylim=c(0,5))+
  scale_color_manual(values = rep("orange",10),guide="none")+
  theme_classic()

NC1<-read.csv("Norway_con1.csv")
NC2<-read.csv("Norway_con2.csv")
NC3<-read.csv("Norway_con3.csv")
NC4<-read.csv("Norway_con4.csv")
NC5<-read.csv("Norway_con5.csv")
NC6<-read.csv("Norway_con6.csv")
NC7<-read.csv("Norway_con7.csv")
NC8<-read.csv("Norway_con8.csv")
NC9<-read.csv("Norway_con9.csv")
NC10<-read.csv("Norway_con10.csv")

NC1$replicate<-"A"
NC2$replicate<-"B"
NC3$replicate<-"C"
NC4$replicate<-"D"
NC5$replicate<-"E"
NC6$replicate<-"F"
NC7$replicate<-"G"
NC8$replicate<-"H"
NC9$replicate<-"I"
NC10$replicate<-"J"

NC<-rbind(NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8,NC9,NC10)

DC1<-read.csv("Denmark_con1.csv")
DC2<-read.csv("Denmark_con2.csv")
DC3<-read.csv("Denmark_con3.csv")
DC4<-read.csv("Denmark_con4.csv")
DC5<-read.csv("Denmark_con5.csv")
DC6<-read.csv("Denmark_con6.csv")
DC7<-read.csv("Denmark_con7.csv")
DC8<-read.csv("Denmark_con8.csv")
DC9<-read.csv("Denmark_con9.csv")
DC10<-read.csv("Denmark_con10.csv")

DC1$replicate<-"A"
DC2$replicate<-"B"
DC3$replicate<-"C"
DC4$replicate<-"D"
DC5$replicate<-"E"
DC6$replicate<-"F"
DC7$replicate<-"G"
DC8$replicate<-"H"
DC9$replicate<-"I"
DC10$replicate<-"J"

DC<-rbind(DC1,DC2,DC3,DC4,DC5,DC6,DC7,DC8,DC9,DC10)

Con <- rbind(IC,GC,NC,DC)

ggplot(Con[Con$x>100,])+
  geom_line(aes(x,y/1000,colour=label,lty=replicate),lwd=1.5,alpha=0.1)+
  xlab("Years ago")+
  ylab("Ne/1000")+
  coord_cartesian(xlim=c(0,300000))+
  scale_linetype_manual(values = rep("solid",10),guide="none")+
  scale_color_manual(name="Sample",values = c("blue","green","red","orange"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme_classic()

setwd("/Users/akijarl/Documents/UI/rannsoknir/Haliaeetus/Haliaeetus_ernir/smcpp_output/nocubic/")

GC1<-read.csv("Greenland_con1.csv")
GC2<-read.csv("Greenland_con2.csv")
GC3<-read.csv("Greenland_con3.csv")
GC4<-read.csv("Greenland_con4.csv")
GC5<-read.csv("Greenland_con5.csv")
GC6<-read.csv("Greenland_con6.csv")
GC7<-read.csv("Greenland_con7.csv")
GC8<-read.csv("Greenland_con8.csv")
GC9<-read.csv("Greenland_con9.csv")
GC10<-read.csv("Greenland_con10.csv")

GC1$replicate<-"A"
GC2$replicate<-"B"
GC3$replicate<-"C"
GC4$replicate<-"D"
GC5$replicate<-"E"
GC6$replicate<-"F"
GC7$replicate<-"G"
GC8$replicate<-"H"
GC9$replicate<-"I"
GC10$replicate<-"J"

GC<-rbind(GC1,GC2,GC3,GC4,GC5,GC6,GC7,GC8,GC9,GC10)

IC1<-read.csv("Iceland_con1.csv")
IC2<-read.csv("Iceland_con2.csv")
IC3<-read.csv("Iceland_con3.csv")
IC4<-read.csv("Iceland_con4.csv")
IC5<-read.csv("Iceland_con5.csv")
IC6<-read.csv("Iceland_con6.csv")
IC7<-read.csv("Iceland_con7.csv")
IC8<-read.csv("Iceland_con8.csv")
IC9<-read.csv("Iceland_con9.csv")
IC10<-read.csv("Iceland_con10.csv")

IC1$replicate<-"A"
IC2$replicate<-"B"
IC3$replicate<-"C"
IC4$replicate<-"D"
IC5$replicate<-"E"
IC6$replicate<-"F"
IC7$replicate<-"G"
IC8$replicate<-"H"
IC9$replicate<-"I"
IC10$replicate<-"J"

IC<-rbind(IC1,IC2,IC3,IC4,IC5,IC6,IC7,IC8,IC9,IC10)
#IC<-rbind(IC1,IC2,IC4,IC5,IC6,IC9,IC10)

NC1<-read.csv("Norway_con1.csv")
NC2<-read.csv("Norway_con2.csv")
NC3<-read.csv("Norway_con3.csv")
NC4<-read.csv("Norway_con4.csv")
NC5<-read.csv("Norway_con5.csv")
NC6<-read.csv("Norway_con6.csv")
NC7<-read.csv("Norway_con7.csv")
NC8<-read.csv("Norway_con8.csv")
NC9<-read.csv("Norway_con9.csv")
NC10<-read.csv("Norway_con10.csv")

NC1$replicate<-"A"
NC2$replicate<-"B"
NC3$replicate<-"C"
NC4$replicate<-"D"
NC5$replicate<-"E"
NC6$replicate<-"F"
NC7$replicate<-"G"
NC8$replicate<-"H"
NC9$replicate<-"I"
NC10$replicate<-"J"

NC<-rbind(NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8,NC9,NC10)

DC1<-read.csv("Denmark_con1.csv")
DC2<-read.csv("Denmark_con2.csv")
DC3<-read.csv("Denmark_con3.csv")
DC4<-read.csv("Denmark_con4.csv")
DC5<-read.csv("Denmark_con5.csv")
DC6<-read.csv("Denmark_con6.csv")
DC7<-read.csv("Denmark_con7.csv")
DC8<-read.csv("Denmark_con8.csv")
DC9<-read.csv("Denmark_con9.csv")
DC10<-read.csv("Denmark_con10.csv")

DC1$replicate<-"A"
DC2$replicate<-"B"
DC3$replicate<-"C"
DC4$replicate<-"D"
DC5$replicate<-"E"
DC6$replicate<-"F"
DC7$replicate<-"G"
DC8$replicate<-"H"
DC9$replicate<-"I"
DC10$replicate<-"J"

DC<-rbind(DC1,DC2,DC3,DC4,DC5,DC6,DC7,DC8,DC9,DC10)

EC1<-read.csv("Estonia_con1.csv")
EC2<-read.csv("Estonia_con2.csv")
EC3<-read.csv("Estonia_con3.csv")
EC4<-read.csv("Estonia_con4.csv")
EC5<-read.csv("Estonia_con5.csv")
EC6<-read.csv("Estonia_con6.csv")
EC7<-read.csv("Estonia_con7.csv")
EC8<-read.csv("Estonia_con8.csv")
EC9<-read.csv("Estonia_con9.csv")
EC10<-read.csv("Estonia_con10.csv")

EC1$replicate<-"A"
EC2$replicate<-"B"
EC3$replicate<-"C"
EC4$replicate<-"D"
EC5$replicate<-"E"
EC6$replicate<-"F"
EC7$replicate<-"G"
EC8$replicate<-"H"
EC9$replicate<-"I"
EC10$replicate<-"J"

EC<-rbind(EC1,EC2,EC3,EC4,EC5,EC6,EC7,EC8,EC9,EC10)

Con <- rbind(IC,GC,NC,DC,EC)

ggplot(Con[Con$x>100,])+
  geom_line(aes(x*1.56,(y/1000),colour=label,lty=replicate),lwd=1.5,alpha=0.1)+
  xlab("Years ago")+
  ylab("Ne/1000")+
  coord_cartesian(xlim=c(0,50000))+#,ylim=c(0,23))+
  scale_linetype_manual(values = rep("solid",10),guide="none")+
  scale_color_manual(name="Sample",values = c("blue","black","green","red","orange"), labels=c("Denmark","Estonia","Greenland","Iceland","Norway"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  scale_x_continuous(breaks = c(seq(1000,50000,2000)))+
  geom_vline(xintercept=10000,lty=1)+
  geom_vline(xintercept=25000,lty=2)+
  geom_vline(xintercept=110000,lty=3)+
  geom_vline(xintercept=150000,lty=3)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))

#Historic
setwd("/Users/akijarl/Documents/UI/rannsoknir/Haliaeetus/Haliaeetus_ernir/smcpp_output/nocubic/His/")

GH1<-read.csv("Greenland_his1.csv")
GH2<-read.csv("Greenland_his2.csv")
GH3<-read.csv("Greenland_his3.csv")
GH4<-read.csv("Greenland_his4.csv")
GH5<-read.csv("Greenland_his5.csv")
GH6<-read.csv("Greenland_his6.csv")
GH7<-read.csv("Greenland_his7.csv")
GH8<-read.csv("Greenland_his8.csv")
GH9<-read.csv("Greenland_his9.csv")
GH10<-read.csv("Greenland_his10.csv")

GH1$replicate<-"A"
GH2$replicate<-"B"
GH3$replicate<-"C"
GH4$replicate<-"D"
GH5$replicate<-"E"
GH6$replicate<-"F"
GH7$replicate<-"G"
GH8$replicate<-"H"
GH9$replicate<-"I"
GH10$replicate<-"J"

GH<-rbind(GH1,GH2,GH3,GH4,GH5,GH6,GH7,GH8,GH9,GH10)

IH1<-read.csv("Iceland_his1.csv")
IH2<-read.csv("Iceland_his2.csv")
IH3<-read.csv("Iceland_his3.csv")
IH4<-read.csv("Iceland_his4.csv")
IH5<-read.csv("Iceland_his5.csv")
IH6<-read.csv("Iceland_his6.csv")
IH7<-read.csv("Iceland_his7.csv")
IH8<-read.csv("Iceland_his8.csv")
IH9<-read.csv("Iceland_his9.csv")
IH10<-read.csv("Iceland_his10.csv")

IH1$replicate<-"A"
IH2$replicate<-"B"
IH3$replicate<-"C"
IH4$replicate<-"D"
IH5$replicate<-"E"
IH6$replicate<-"F"
IH7$replicate<-"G"
IH8$replicate<-"H"
IH9$replicate<-"I"
IH10$replicate<-"J"

IH<-rbind(IH1,IH2,IH3,IH4,IH5,IH6,IH7,IH8,IH9,IH10)

NH1<-read.csv("Norway_his1.csv")
NH2<-read.csv("Norway_his2.csv")
NH3<-read.csv("Norway_his3.csv")
NH4<-read.csv("Norway_his4.csv")
NH5<-read.csv("Norway_his5.csv")
NH6<-read.csv("Norway_his6.csv")
NH7<-read.csv("Norway_his7.csv")
NH8<-read.csv("Norway_his8.csv")
NH9<-read.csv("Norway_his9.csv")
NH10<-read.csv("Norway_his10.csv")

NH1$replicate<-"A"
NH2$replicate<-"B"
NH3$replicate<-"C"
NH4$replicate<-"D"
NH5$replicate<-"E"
NH6$replicate<-"F"
NH7$replicate<-"G"
NH8$replicate<-"H"
NH9$replicate<-"I"
NH10$replicate<-"J"

NH<-rbind(NH1,NH2,NH3,NH4,NH5,NH6,NH7,NH8,NH9,NH10)

DH1<-read.csv("Denmark_his1.csv") 
DH2<-read.csv("Denmark_his2.csv")
DH3<-read.csv("Denmark_his3.csv")
DH4<-read.csv("Denmark_his4.csv")
DH5<-read.csv("Denmark_his5.csv")
DH6<-read.csv("Denmark_his6.csv")
DH7<-read.csv("Denmark_his7.csv")
DH8<-read.csv("Denmark_his8.csv")
DH9<-read.csv("Denmark_his9.csv")
DH10<-read.csv("Denmark_his10.csv")

DH1$replicate<-"A"
DH2$replicate<-"B"
DH3$replicate<-"C"
DH4$replicate<-"D"
DH5$replicate<-"E"
DH6$replicate<-"F"
DH7$replicate<-"G"
DH8$replicate<-"H"
DH9$replicate<-"I"
DH10$replicate<-"J"

DH<-rbind(DH1,DH2,DH3,DH4,DH5,DH6,DH7,DH8,DH9,DH10)

TH1<-read.csv("Turkey_his1.csv") 
TH2<-read.csv("Turkey_his2.csv")
TH3<-read.csv("Turkey_his3.csv")
TH4<-read.csv("Turkey_his4.csv")
TH5<-read.csv("Turkey_his5.csv")
TH6<-read.csv("Turkey_his6.csv")
TH7<-read.csv("Turkey_his7.csv")
TH8<-read.csv("Turkey_his8.csv")
TH9<-read.csv("Turkey_his9.csv")
TH10<-read.csv("Turkey_his10.csv")

TH1$replicate<-"A"
TH2$replicate<-"B"
TH3$replicate<-"C"
TH4$replicate<-"D"
TH5$replicate<-"E"
TH6$replicate<-"F"
TH7$replicate<-"G"
TH8$replicate<-"H"
TH9$replicate<-"I"
TH10$replicate<-"J"

TH<-rbind(TH1,TH2,TH3,TH4,TH5,TH6,TH7,TH8,TH9,TH10)

His <- rbind(IH,GH,NH,DH,TH)

ggplot(His[His$x>100,])+
  geom_line(aes(x,y/1000,colour=label,lty=replicate),lwd=1.5,alpha=0.1)+
  xlab("Years ago")+
  ylab("Ne/1000")+
  #coord_cartesian(xlim=c(0,200000),ylim=c(0,10))+
  scale_linetype_manual(values = rep("solid",10),guide="none")+
  scale_color_manual(name="Sample",values = c("blue","green","red","orange","black"), labels=c("Denmark","Greenland","Iceland","Norway","Turkey"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme_classic()

ggplot(His[His$x>100,])+
  geom_line(aes(x,(y/1000),colour=label,lty=replicate),lwd=1.5,alpha=0.1)+
  xlab("Years ago")+
  ylab("Ne/1000")+
  coord_cartesian(xlim=c(0,350000),ylim=c(0,30))+
  scale_linetype_manual(values = rep("solid",10),guide="none")+
  scale_color_manual(name="Sample",values = c("blue","green","red","orange","black"), labels=c("Denmark","Greenland","Iceland","Norway","Turkey"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  scale_x_continuous(breaks = c(seq(1000,351000,50000)))+
  geom_vline(xintercept=10000,lty=1)+
  geom_vline(xintercept=25000,lty=2)+
  geom_vline(xintercept=110000,lty=3)+
  geom_vline(xintercept=150000,lty=3)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))

ggplot(His[His$x>100,])+
  geom_line(aes(x*1.56,(y/1000),colour=label,lty=replicate),lwd=1.5,alpha=0.1)+
  xlab("Years ago")+
  ylab("Ne/1000")+
  coord_cartesian(xlim=c(0,150000),ylim=c(0,30))+
  scale_linetype_manual(values = rep("solid",10),guide="none")+
  scale_color_manual(name="Sample",values = c("blue","green","red","orange","black"), labels=c("Denmark","Greenland","Iceland","Norway","Turkey"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  scale_x_continuous(breaks = c(seq(1000,151000,5000)))+
  geom_vline(xintercept=10000,lty=1)+
  geom_vline(xintercept=25000,lty=2)+
  geom_vline(xintercept=110000,lty=3)+
  geom_vline(xintercept=150000,lty=3)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))


setwd("/home/aki/Documents/Rannsóknir/Haaliaeetus/Haliaeetus_ernir/smcpp_output/nocubic/Split/")
require(ggplot2)

IM1<-read.csv("Mall_IC_split_g10_1.csv")
IM2<-read.csv("Mall_IC_split_g10_2.csv")
IM3<-read.csv("Mall_IC_split_g10_3.csv")
IM4<-read.csv("Mall_IC_split_g10_4.csv")
IM5<-read.csv("Mall_IC_split_g10_5.csv")
IM6<-read.csv("Mall_IC_split_g10_6.csv")
IM7<-read.csv("Mall_IC_split_g10_7.csv")
IM8<-read.csv("Mall_IC_split_g10_8.csv")
IM9<-read.csv("Mall_IC_split_g10_9.csv")
IM10<-read.csv("Mainland_IC_split_g10.csv")

IM1$replicate<-"A"
IM2$replicate<-"B"
IM3$replicate<-"C"
IM4$replicate<-"D"
IM5$replicate<-"E"
IM6$replicate<-"F"
IM7$replicate<-"G"
IM8$replicate<-"H"
IM9$replicate<-"I"
IM10$replicate<-"J"

IM<-rbind(IM1,IM2,IM3,IM4,IM5,IM6,IM7,IM8,IM9,IM10)

ggplot(IM[IM$x>100,])+
  geom_line(aes(x,(y/1000),colour=label,lty=replicate),lwd=1.5,alpha=0.1)+
  xlab("Years ago")+
  ylab("Ne/1000")+
  coord_cartesian(xlim=c(0,35000),ylim=c(0,30))+
  scale_linetype_manual(values = rep("solid",10),guide="none")+
  scale_color_manual(name="Sample",values = c("red","black"), labels=c("Iceland","Mainland (all)"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  scale_x_continuous(breaks = c(100,seq(1000,36000,1000)))+
  geom_vline(xintercept=max(IM1[IM1$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(IM2[IM2$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(IM3[IM3$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(IM4[IM4$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(IM5[IM5$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(IM6[IM6$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(IM7[IM7$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(IM8[IM8$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(IM9[IM9$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(IM10[IM10$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))


ID1<-read.csv("ECDHTH_IC_split_g10.2.csv")
ID2<-read.csv("ECDHTH_IC_split1_g10.csv")
ID3<-read.csv("ECDHTH_IC_split2_g10.csv")
ID4<-read.csv("ECDHTH_IC_split3_g10.csv")
ID5<-read.csv("ECDHTH_IC_split4_g10.csv")
ID6<-read.csv("ECDHTH_IC_split5_g10.csv")
ID7<-read.csv("ECDHTH_IC_split6_g10.csv")
ID8<-read.csv("ECDHTH_IC_split7_g10.csv")
ID9<-read.csv("ECDHTH_IC_split8_g10.csv")
ID10<-read.csv("ECDHTH_IC_split9_g10.csv")

ID1$replicate<-"A"
ID2$replicate<-"B"
ID3$replicate<-"C"
ID4$replicate<-"D"
ID5$replicate<-"E"
ID6$replicate<-"F"
ID7$replicate<-"G"
ID8$replicate<-"H"
ID9$replicate<-"I"
ID10$replicate<-"J"

ID<-rbind(ID1,ID2,ID3,ID4,ID5,ID6,ID7,ID8,ID9,ID10)

ggplot(ID[ID$x>100,])+
  geom_line(aes(x,(y/1000),colour=label,lty=replicate),lwd=1.5,alpha=0.1)+
  xlab("Years ago")+
  ylab("Ne/1000")+
  coord_cartesian(xlim=c(0,35000),ylim=c(0,30))+
  scale_linetype_manual(values = rep("solid",10),guide="none")+
  scale_color_manual(name="Sample",values = c("red","black"), labels=c("Iceland","Central Europe (EE,DH,TH)"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  scale_x_continuous(breaks = c(100,seq(1000,36000,1000)))+
  geom_vline(xintercept=max(ID1[ID1$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(ID2[ID2$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(ID3[ID3$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(ID4[ID4$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(ID5[ID5$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(ID6[ID6$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(ID7[ID7$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(ID8[ID8$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(ID9[ID9$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(ID10[ID10$label=="Iceland_con",]$x),lwd=1.5,alpha=0.1,color="blue")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))




GD1<-read.csv("ECDHTH_GC_split_g10_1.csv")
GD2<-read.csv("ECDHTH_GC_split_g10_2.csv")
GD3<-read.csv("ECDHTH_GC_split_g10_3.csv")
GD4<-read.csv("ECDHTH_GC_split_g10_4.csv")
GD5<-read.csv("ECDHTH_GC_split_g10_5.csv")
GD6<-read.csv("ECDHTH_GC_split_g10_6.csv")
GD7<-read.csv("ECDHTH_GC_split_g10_7.csv")
GD8<-read.csv("ECDHTH_GC_split_g10_8.csv")
GD9<-read.csv("ECDHTH_GC_split_g10_9.csv")
GD10<-read.csv("ECDHTH_GC_split_g10_10.csv")

GD1$replicate<-"A"
GD2$replicate<-"B"
GD3$replicate<-"C"
GD4$replicate<-"D"
GD5$replicate<-"E"
GD6$replicate<-"F"
GD7$replicate<-"G"
GD8$replicate<-"H"
GD9$replicate<-"I"
GD10$replicate<-"J"

GD<-rbind(GD1,GD2,GD3,GD4,GD5,GD6,GD7,GD8,GD9,GD10)

ggplot(GD[GD$x>100,])+
  geom_line(aes(x,(y/1000),colour=label,lty=replicate),lwd=1.5,alpha=0.1)+
  xlab("Years ago")+
  ylab("Ne/1000")+
  coord_cartesian(xlim=c(0,35000),ylim=c(0,30))+
  scale_linetype_manual(values = rep("solid",10),guide="none")+
  scale_color_manual(name="Sample",values = c("green","black"), labels=c("Greenland","Central Europe (EE,DH,TH)"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  scale_x_continuous(breaks = c(100,seq(1000,36000,1000)))+
  geom_vline(xintercept=max(GD1[GD1$label=="Greenland_contemporary",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(GD2[GD2$label=="Greenland_contemporary",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(GD3[GD3$label=="Greenland_contemporary",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(GD4[GD4$label=="Greenland_contemporary",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(GD5[GD5$label=="Greenland_contemporary",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(GD6[GD6$label=="Greenland_contemporary",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(GD7[GD7$label=="Greenland_contemporary",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(GD8[GD8$label=="Greenland_contemporary",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(GD9[GD9$label=="Greenland_contemporary",]$x),lwd=1.5,alpha=0.1,color="blue")+
  geom_vline(xintercept=max(GD10[GD10$label=="Greenland_contemporary",]$x),lwd=1.5,alpha=0.1,color="blue")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))
