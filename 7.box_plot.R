setwd("")
diff<-read.csv("diff_sig_genus_fecal_all_fdr_abuntop25.csv",header = T)

diff<-diff[,-2]
diff2<-as.data.frame(t(diff))
colnames(diff2)<-diff$Genus
diff2<-diff2[-1,]
diff2$Group<-data6$Group

data_frame=melt(diff, id='Genus')
names(data_frame)[2]='id'
data_frame$Group<-"Health"
data_frame[301:600,4] <- "Disease"
data_frame$num<-log10(data_frame$value)
data_frame$num[is.infinite(data_frame$num)]=0


#倒下的箱型图
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library("ggsci")
library("ggplot2")
library("gridExtra")


p<-ggplot(data_frame,aes(x=Genus,y=log10(value)))+
  stat_boxplot(geom="errorbar",position = position_dodge(0.9),width=0.3,aes(color=Group))+
  geom_boxplot(aes(fill=Group),position = position_dodge(0.9),lwd=0.6)+#注意画图顺序，如果反过来会出现十字先stat再geom
  theme_test()+
  xlab("")+coord_flip()+labs(y="Log10(Relative abundance)")+
  theme(legend.position="right",
        panel.grid = element_line(color = 'gray', linetype = 0, size = 0.2), 
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 15,face = "bold"),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.border = element_rect(fill=NA,color="black", size=1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 1),
        legend.text=element_text(size=10),
        legend.title=element_blank())
p+scale_fill_d3()+scale_color_manual(values = c("#000000","#000000"))
