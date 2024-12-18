setwd("")
library(tidyr)
library(ggplot2)
library(ggsci)
library(xlsx)
library(dplyr)
library(picante)
library(ggpubr)
cli1<-read.xlsx("1.biochemical.xlsx",5,encoding = "UTF-8")
cli1$Group<-factor(cli1$Group,levels = c("LD","Cg","Rg","Fp"))

pdf("fig5c.pdf",
    width = 1.1,height = 1.4)

p<-ggplot(cli1, aes(x = Group, y = value, fill = factor(Group))) +
  geom_bar(position = position_dodge(0.9), stat = "summary", fun = "mean", color = NA) +
  geom_errorbar(position = position_dodge(0.9), stat = "summary", fun.data = mean_se, width = 0.2,aes(color = factor(Group))) +
  #geom_jitter(color=cli1$col,fill=NA,width =0.2,shape =cli1$Group, size = 0.65)+
  theme_classic() +
  xlab("")+labs(y="-△△CT")+
  theme(legend.position="none",
        panel.grid = element_line(color = 'gray', linetype = 0, size = 0.2), 
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 7,face = "bold"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        legend.text=element_text(size=9),
        legend.title=element_blank(),
        plot.margin =unit(c(0.05,0,0,0), "cm"))
p1 <- p + scale_fill_manual(values = c("#B2B9C1","#35395A","#E4B6B4","#C0C880"))
p2=p1+scale_y_continuous(breaks = c(0, 2, 4, 6,8,10,12),expand=c(0,0))+expand_limits(y=c(0,12))+
  scale_color_manual(values =  c("#B2B9C1","#35395A","#E4B6B4","#C0C880"))
p2

dev.off()
