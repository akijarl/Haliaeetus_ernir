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

# plot(GC1$x,log10(GC1$y),type="l", col=alpha(rgb(0,0,0), 0.2), xlim = c(100,400000), ylim = c(2,5),lwd=2, xlab = "Years", ylab = "Effective population size")
# lines(GC2$x,log10(GC2$y), lwd=2, col=alpha(rgb(0,0,0), 0.2))
# lines(GC3$x,log10(GC3$y), lwd=2, col=alpha(rgb(0,0,0), 0.2))
# lines(GC4$x,log10(GC4$y), lwd=2, col=alpha(rgb(0,0,0), 0.2))
# lines(GC5$x,log10(GC5$y), lwd=2, col=alpha(rgb(0,0,0), 0.2))
# lines(GC6$x,log10(GC6$y), lwd=2, col=alpha(rgb(0,0,0), 0.2))
# lines(GC7$x,log10(GC7$y), lwd=2, col=alpha(rgb(0,0,0), 0.2))
# lines(GC8$x,log10(GC8$y), lwd=2, col=alpha(rgb(0,0,0), 0.2))
# lines(GC9$x,log10(GC9$y), lwd=2, col=alpha(rgb(0,0,0), 0.2))
# lines(GC10$x,log10(GC10$y), lwd=2, col=alpha(rgb(0,0,0), 0.2))

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
#GC<-rbind(GC2,GC3,GC6)
#GC<-rbind(GC1,GC4,GC5,GC7,GC8,GC9,GC10)

#summary(GC$x)
#[GC$x>100,]
ggplot(GC)+
  geom_line(aes(x,log10(y),colour=replicate))+
  #coord_cartesian(xlim=c(0,100000),ylim=c(0,5))+
  scale_color_manual(values = rep("orange",10),guide="none")+
  theme_classic()

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

