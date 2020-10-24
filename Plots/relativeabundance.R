library(readxl)
library(ggplot2)
library(ggpubr)
GXcommspecific_abd <- read_excel("GXcommspecific_abd.xlsx")
GXminimal_abd <- read_excel("GXminimal_abd.xlsx")
GXexcess_abd <- read_excel("GXexcess_abd.xlsx")
View(GXcustom_abd)

Minimal <- ggplot(GXminimal_abd, aes(x=Species,y=RelativeAbundance,fill=Species),backgroundColor = "white")+geom_boxplot(color="black", fill="green", alpha=0.3)+ stat_summary(fun=mean, geom="point", shape=18, size=2, color="black")+
  ggtitle("Abundances of LAB species in coculture with minimial nutrients")+theme(plot.title = element_text(hjust = 0.5, size=12))+
   ylab("Relative Abundances")+ theme(legend.position="none")+theme(axis.text.x = element_blank(), axis.title.x=element_blank())+
  theme(axis.text.y = element_text(face = "bold", size = 10))+theme(axis.title.y = element_text(size = 10))+ theme(panel.background = element_rect(fill="white"),axis.line = element_line(colour = "black")) 

Commspecific <- ggplot(GXcommspecific_abd, aes(x=Species,y=RelativeAbundance,fill=Species),backgroundColor = "white")+geom_boxplot(color="black", fill="red", alpha=0.3)+ stat_summary(fun=mean, geom="point", shape=18, size=2, color="black")+
  ggtitle("Abundances of LAB species in coculture with community-specific nutrients")+theme(plot.title = element_text(hjust = 0.5, size=12))+theme(axis.text.x = element_blank(), axis.title.x=element_blank())+
  ylab("Relative Abundances")+ theme(legend.position="none")+theme(axis.text.y = element_text(face = "bold", size = 10))+theme(axis.title.y = element_text(size = 10))+theme(axis.title.x = element_text(size = 10))+ theme(panel.background = element_rect(fill="white"),axis.line = element_line(colour = "black")) 

Excess <- ggplot(GXexcess_abd, aes(x=Species,y=RelativeAbundance,fill=Species),backgroundColor = "white")+geom_boxplot(color="black", fill="blue", alpha=0.3)+ stat_summary(fun=mean, geom="point", shape=18, size=2, color="black")+
  ggtitle("Abundances of LAB species in coculture with excess nutrients")+theme(plot.title = element_text(hjust = 0.5, size=12))+
  xlab("LAB species")+ ylab("Relative Abundances")+ theme(legend.position="none")+theme(axis.text.x = element_text(face = "bold.italic", angle=90, size = 12,hjust = 1))+
  theme(axis.text.y = element_text(face = "bold", size = 10))+theme(axis.title.y = element_text(size = 10))+theme(axis.title.x = element_text(size = 12))+ theme(panel.background = element_rect(fill="white"),axis.line = element_line(colour = "black")) 

figure <- ggarrange(Minimal,Commspecific,Excess,labels = c("A", "B", "C"),font.label = list(size = 10),ncol = 1, nrow = 3,align=c("v"),heights = c(2, 2, 6))
figure 

ggsave(figure,filename="abundance.png", dpi=900)
dev.off

