library(ggrepel)
library(ggplot2)
library(ggsci)
library(ggplot2)
library(gridExtra)
setwd("")


bubble<-read.csv("3.4.genus_bubble.csv",header = TRUE)
bubble<-read.csv("3.4.genus_bubble_2group.csv",header = TRUE)

bubble$Group <- factor(bubble$Group, levels = c('Health', 'Disease',"Bile","Stone"))
bubble$Group <- factor(bubble$Group, levels = c('Health', 'Disease'))

pdf("1.figureC.pdf",width = 3,height = 3) 
p1<-ggplot(data=bubble,aes(x=Group,y=Genus))+
  geom_point(aes(size=Size,fill=Group),shape=21,stroke = 0,alpha=1)+
  theme_test()+
  theme(legend.position="right",
        #panel.grid = element_line(color = 'gray', linetype = 0, size = 0.2), 
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 8,face = "bold"),
        axis.title.y = element_text(size = 8,face = "bold"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        panel.border = element_rect(fill=NA,color="black", size=1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 1),
        legend.text=element_text(size=9),
        legend.title=element_blank())+
  labs(size="Relative abundance (%)",x=" ",y=" ")+
  scale_size_continuous(range = c(1,3.8))+
  scale_fill_simpsons(alpha = 1)
p1
dev.off()  