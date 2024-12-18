setwd("")

library(dplyr)
library(ggplot2)
library(vegan)
library(ade4)
library(factoextra)
library(RColorBrewer)
library(ggsci)

#PCA
data<- read.csv("3.0.ASV_relative_fecal.csv",header=T)
data2<- data[1:119,2:1537]
#match_fecal
data3<-as.data.frame(t(data))
colnames(data3)<-data$ASV
data3<-data3[2:1543,]
data4<- read.table("asv_fecal_meta_match.txt",header=T)
vars<-data4$OUT_ID
data5<-data3[vars]
data6<-as.data.frame(t(data5))
data2<- data6[1:60,2:1536]
data2=as.data.frame(lapply(data2,as.numeric))

write.csv(data6, file= "5.ASV_relative_fecal_match.csv",row.names = T ,quote = F)


res.pca <- dudi.pca(data2,
                    scannf = FALSE,   # Hide scree plot
                    nf = 5           # Number of components kept in the results
)
fviz_eig(res.pca)

pca_eig <- (res.pca$eig)[1:2] / sum(res.pca$eig)
#提取样本点坐标（前两轴）
sample_site <- data.frame({res.pca$li})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCA1', 'PCA2')
sample_site$group<-data6$Group
sample_site$group<-data6$Group_age
sample_site$group<-data6$gender
#sample_site$group<-data$BMI_group_1
#sample_site$group<-data$BMI_group_2
#write.csv(sample_site,"sample_site.csv")

pca_plot <- ggplot(sample_site, aes(PCA1, PCA2,color=group)) +
  #scale_colour_manual(values = c("#00AFBB", "#E7B800"))+
  #scale_fill_manual(values = c("#00AFBB", "#E7B800"))+
  xlim(-8,-1)+ylim(-6,1)+labs(x="PC1 (12.9%)",y="PC2 (6.6%)")+
  theme_test()+#去掉背景框
  geom_point(size = 3)+ #可在这里修改点的透明度、大小
  stat_ellipse(aes(PCA1, PCA2,fill=group),level = 0.95,geom="polygon",linetype=0,alpha=0.2) +
  theme(legend.position=c(0.92,0.90),
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
pca_plot+scale_fill_d3()+scale_color_d3()

pca_plot <- ggplot(sample_site, aes(PCA1, PCA2,color=group)) +
  scale_colour_manual(values = c("#00AFBB", "#E7B800","#DC143C"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800","#DC143C"))+
  xlim(-5.5,-1.5)+ylim(-1.5,3)+labs(x="PC1 (8.8%)",y="PC2 (5.0%)")+
  theme_test()+#去掉背景框
  geom_point(size = 3)+ #可在这里修改点的透明度、大小
  stat_ellipse(aes(PCA1, PCA2,fill=group),level = 0.95,geom="polygon",linetype=0,alpha=0.2) +
  theme(legend.position=c(0.92,0.90),
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

otu.adonis <- adonis(data2 ~ Group_age, data = data6, permutations = 999)
otu.adonis <- adonis(data2 ~ Group, data = data6, permutations = 999)
otu.adonis <- adonis(data2 ~ gender, data = data6, permutations = 999)
otu.adonis <- adonis(data2 ~ BMI_group_1, data = data6, permutations = 999)

#PCA放法2：
#data<- read.table("asv_table_fecal_group.txt",header=T)
#data2<- data[1:119,2:1684]
#otu.pca <- prcomp(data2,scale= TRUE)
#get_eigenvalue(otu.pca)[1:3,]
#p <- fviz_pca_ind(otu.pca,
#     geom.ind = "point", # show points only ( not "text")
 #    col.ind = data$Group_age, # color by groups
  #   palette = c("#00AFBB", "#E7B800","#DC143C"),
   #  pointsize = 1, font.family = "Arial",
    # addEllipses = TRUE, # Concentration ellipses
     #legend.title = "Groups")

#pcoa
#data2

#物种数据 Hellinger 预转化（处理包含很多 0 值的群落物种数据时，推荐使用）
data_hel <- decostand(data2, method = 'hellinger')
bray <- vegdist(data_hel , method = 'bray')
pcoa <- cmdscale(bray, k = (nrow(data_hel ) - 1), eig = TRUE)
pcoa_eig<-data.frame(scores(pcoa)[,1:2])
pcoa_eig$sample<-rownames(pcoa_eig)
pcoa_eig$group<-data6$Group
pcoa_eig1<-round(100*pcoa$eig[1:2]/sum(pcoa$eig),2)

write.csv(pcoa_eig, file= "pcoa_matched_index_bray.csv",row.names = T ,quote = F)

# PERMANOVA分析
# 整体水平比较
adonis_result_dis = adonis(bray~Group, data6, permutations = 999)
adonis_result_dis

library("ggsci")

pdf("sup2.figureB.pdf",width = 1.5,height = 1.2) 
p<-ggplot(pcoa_eig, aes(Dim1, Dim2, color= group)) +
  xlim(-0.3,0.3)+ylim(-0.3,0.4)+
  labs(x = paste('PCoA1: ',pcoa_eig1[1], '%'), y = paste('PCoA2: ', pcoa_eig1[2],'%'),title="Bray-curtis PCoA Analysis")+
  geom_point(size=0.15)+
  stat_ellipse(aes(Dim1, Dim2,fill= group),level = 0.9,geom="polygon",linetype=0,alpha=0.2)+
  theme_test()+
  theme(legend.position="none",
        panel.grid = element_line(color = 'gray', linetype = 0, size = 0.2), 
        panel.background = element_blank(), plot.title = element_text(siz=3),
        axis.title.x = element_text(size = 7,face = "bold"),
        axis.title.y = element_text(size = 7,face = "bold"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        panel.border = element_rect(fill=NA,color="black", size=0.5),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        legend.text=element_text(size=6),
        legend.title=element_blank(),
        legend.key  = element_blank())
p+scale_fill_manual(values = c("#1F77B4","#739884"))+scale_color_manual(values = c("#1F77B4","#739884"))

dev.off()


#有坐标表时
pca<-read.csv("pca.csv",header = TRUE)
library(vegan)
library(gplots)
library(RColorBrewer)
library(plyr)
library(ggplot2)
pca$group <- factor(pca$group, levels = c('C', 'N'))
#scale_color_manual(values = brewer.pal(6,"Set2"))
ggplot(pca, aes(PC1, PC2, color= group,fill=group)) +
  theme_test()+#去掉背景框
  geom_vline(xintercept = 0,color="white",size = 0.4) +
  geom_hline(yintercept = 0,color="white",size = 0.4) +
  geom_point(aes(color = group),size = 3)+ #可在这里修改点的透明度、大小
  stat_ellipse(level = 0.95,geom = "polygon",alpha=1/5)+
  xlim(-0.5,0.6)+ylim(-0.5,0.5)+
  labs(x='PC1(21.03%)', y='PC2(13.36%)')+
  theme(legend.position=c(0.94,0.82),
        panel.grid = element_line(color = 'gray', linetype = 2, size = 0.2), 
        panel.background = element_blank(), 
        legend.title=element_blank())

nmds<-read.csv("nmds.csv",header = TRUE)
nmds$group <- factor(nmds$group, levels = c('N', 'C'))
ggplot(nmds, aes(MDS1,MDS2, color= group,fill=group)) +
  theme_test()+#去掉背景框
  geom_vline(xintercept = 0,color="white",size = 0.4) +
  geom_hline(yintercept = 0,color="white",size = 0.4) +
  geom_point(aes(color = group),size = 3)+ #可在这里修改点的透明度、大小
  stat_ellipse(level = 0.95,geom = "polygon",alpha=1/5)+
  xlim(-0.6,0.4)+ylim(-0.4,0.4)+
  labs(x='NMDS1', y='NMDS2')+
  theme(legend.position=c(0.95,0.90),
        panel.grid = element_line(color = 'gray', linetype = 2, size = 0.2), 
        panel.background = element_blank(), 
        legend.title=element_blank())

pcoa<-read.csv("pcoa2.csv",header = TRUE)
pcoa$group <- factor(pcoa$group, levels = c('N', 'C'))
ggplot(pcoa, aes(PC1, PC2, color= group,fill=group)) +
  theme_test()+#去掉背景框
  geom_vline(xintercept = 0,color="white",size = 0.4) +
  geom_hline(yintercept = 0,color="white",size = 0.4) +
  geom_point(aes(color = group),size = 3)+ #可在这里修改点的透明度、大小
  stat_ellipse(level = 0.95,geom = "polygon",alpha=1/5)+
  xlim(-0.6,-0.1)+ylim(-0.8,0.6)+
  labs(x='PC1(31.68%)', y='PC2(11.98%)')+
  theme(legend.position=c(0.94,0.82),
        panel.grid = element_line(color = 'gray', linetype = 2, size = 0.2), 
        panel.background = element_blank(), 
        legend.title=element_blank())