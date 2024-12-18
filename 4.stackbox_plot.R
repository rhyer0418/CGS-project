setwd("")
group<-read.csv("1.group.csv",header = T)

library(reshape2)
library(ggplot2)
library(randomcoloR)
library(ggsci)
library(ggplot2)
library(gridExtra)
library(showtext)


data<- read.csv("3.3.phylum_relative_5.csv",header=T)

#所有样本的堆叠图
id=data$Taxonomy
data_g=data.frame(id)
data_g$group=data$group
phylum=data[2:6]
Taxonomy=colnames(phylum)
data_frame=data.frame(t(phylum), Taxonomy)
data_g$id<-as.character(data_g$id)
names(data_frame)[1:191]=data_g$id
data_frame<-data_frame[1:5,]
data_frame=melt(data_frame, id='Taxonomy')
names(data_frame)[2]='id'
data_frame=merge(data_frame, data_g, by='id')
data_frame$id <- as.character(data_frame$id)
data_frame$value<-as.numeric(data_frame$value)
a<-tapply(data_frame$value,data_frame$id,sum)

data_frame$Taxonomy=factor(data_frame$Taxonomy, levels=c('Firmicutes','Bacteroidota','Proteobacteria','Actinobacteriota',
                                                      "Others"))
p1=ggplot(data_frame, aes(x=reorder(id,-as.numeric(Taxonomy=='Firmicutes')*value/a[id]), fill=Taxonomy,value))+
  geom_col(position='stack') +
  labs(x='Samples', y='Relative Abundance (%)')+
  scale_y_continuous(expand=c(0, 0))+
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())+
  theme(legend.position='right',
        panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_blank(), 
        legend.title=element_blank())
p1_npg=p1+scale_fill_simpsons()
p1_npg




data<- read.csv("3.3.phylum_relative_5.csv",header=T)

#分组样本的堆叠图
id=data$Phylum
data_g=data.frame(id)
data_g$group=data$Phylum
phylum=data[2:6]
Taxonomy=colnames(phylum)
data_frame=data.frame(t(phylum), Taxonomy)
data_g$id<-as.character(data_g$id)
names(data_frame)[1:4]=data_g$id
data_frame<-data_frame[1:5,]
data_frame=melt(data_frame, id='Taxonomy')
names(data_frame)[2]='id'
data_frame=merge(data_frame, data_g, by='id')
data_frame$id <- as.character(data_frame$id)
data_frame$value<-as.numeric(data_frame$value)
a<-tapply(data_frame$value,data_frame$id,sum)

data_frame$Taxonomy=factor(data_frame$Taxonomy, levels=c('Firmicutes','Bacteroidota','Proteobacteria','Actinobacteriota',
                                                         "Others"))

pdf("4.phylum_4group.pdf",width = 3,height = 1)  ## 必须保存；pdf、png、jpg
#showtext_begin()  #必须有，制定字体

p1=ggplot(data_frame, aes(x=reorder(id,-as.numeric(Taxonomy=='Firmicutes')*value/a[id]), fill=Taxonomy,value))+
  geom_col(position='stack', width = 0.8) +
  labs(x='', y='Relative Abundance (%)')+
  coord_flip() + #倒置
  #scale_fill_manual(values = c("#b5182b","#f9cb45","#34a186","#DC143C","#00AFBB"))+
  theme_classic()+
  scale_y_continuous(expand=c(0, 0))+
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())+
  theme(legend.position="right",
        panel.grid = element_line(color = 'gray', linetype = 0, size = 0.2), 
        #panel.border = element_rect(fill=NA,colour = 'black',size=2),
        panel.background = element_blank(), 
        #text=element_text(family = "Arial"),
        axis.title.x = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(size = 10,face = "bold"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 1),
        legend.text=element_text(size=7),
        legend.title=element_text(size = 8),
        legend.key.size = unit(7, "points"),  # 设置图例色框
        plot.margin = margin(0, 0, 0, 0.1, "cm")  # 设置左边空白
        )
p2<-p1+scale_fill_simpsons(alpha = 0.9)
p2
#showtext_end();  ##必须有
dev.off()  ##必须有

#genus
data<- read.csv("3.4.genus_relative_4group_top10.csv",header=T)

id=data$Genus
data_g=data.frame(id)
data_g$group=data$Genus
phylum=data[2:26]
Taxonomy=colnames(phylum)
data_frame=data.frame(t(phylum), Taxonomy)
data_g$id<-as.character(data_g$id)
names(data_frame)[1:4]=data_g$id
data_frame<-data_frame[1:25,]
data_frame=melt(data_frame, id='Taxonomy')
names(data_frame)[2]='id'
data_frame=merge(data_frame, data_g, by='id')
data_frame$id <- as.character(data_frame$id)
data_frame$value<-as.numeric(data_frame$value)
a<-tapply(data_frame$value,data_frame$id,sum)

data_frame$Taxonomy=factor(data_frame$Taxonomy, levels=c("Bacteroides","Faecalibacterium","Bifidobacterium","Blautia",
                                                         "Escherichia.Shigella","Agathobacter","Others","Dialister","Subdoligranulum",
                                                         "Megamonas","Fusicatenibacter","Akkermansia","Parabacteroides","Enterococcus",
                                                         "Muribaculaceae","Mitochondria","Morganella","Acinetobacter","Pseudomonas",
                                                         "Lactococcus","Cutibacterium","Ralstonia","Comamonas","Vibrio","Castellaniella"))
library(RColorBrewer)
library(ggsci)
library('scales')
#show_col(colorRampPalette(pal_nejm('default',alpha = 0.8)(8))(8))#展示所选颜色 
#mycolor = colorRampPalette(pal_nejm('default',alpha = 0.8)(8))(8)#返回颜色值
#mycolor = colorRampPalette(pal_npg('nrc',alpha = 0.8)(10))(10)#返回颜色值
#mycolor = colorRampPalette(pal_locuszoom('default',alpha = 0.8)(7))(7)#返回颜色值


p1=ggplot(data_frame, aes(x=reorder(id,-as.numeric(Taxonomy=='Bacteroides')*value/a[id]), fill=Taxonomy,value))+
  geom_col(position='stack') +
  labs(x='', y='Relative Abundance (%)')+
  scale_fill_manual(values = c("#BC3C29","#0072B5","#E18727","#20854E","#7876B1","#6F99AD","#FFDC91","#EE4C97","#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85","#D43F3A","#EEA236","#5CB85C","#46B8DA","#357EBD","#9632B8","#B8B8B8"))+
  theme_classic()+
  scale_y_continuous(expand=c(0, 0))+
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())+
  theme(legend.position="right",
        panel.grid = element_line(color = 'gray', linetype = 0, size = 0.2), 
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 15,face = "bold"),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        #panel.border = element_rect(fill=NA,color="black", size=1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 1),
        legend.text=element_text(size=10),
        legend.title=element_blank())
p1




#带起始点折现变化的折线图
library(reshape2)
library(ggplot2)
library(randomcoloR)
library(tidyverse)
phylum<- read.csv("phlum.csv",header=T)
colnames(phylum)[1]<-"phylum"
# 计算连线起始点Y轴坐标，即累计丰度的值
link_dat <- phylum %>% 
  arrange(by=desc(phylum)) %>% 
  mutate(N=cumsum(N), C=cumsum(C))

# 数据格式转换，宽表格转换为ggplot2使用的长表格
df.long <- phylum %>% gather(group, abundance, -phylum)
df.long$group <- factor(df.long$group, levels = c('N', 'C'))

ggplot(df.long, aes(group, fill=phylum,abundance)) + 
  geom_col(position='stack',width = 0.5) + 
  scale_fill_manual(values = rev(distinctColorPalette(7))) +
  labs(x='Group', y='Relative Abundance (%)')+
  theme_classic()+
  #geom_segment(data=link_dat,aes(x=1.25, xend=1.75, y=N, yend=C))+
  theme(legend.position='bottom',
      panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
      panel.background = element_blank(), 
      legend.title=element_blank())

