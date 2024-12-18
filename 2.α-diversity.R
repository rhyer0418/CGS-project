library(vegan)	#用于计算 Shannon 熵指数、Simpson 指数、Chao1 指数、ACE 指数等，同时用于抽样
library(picante)	#用于计算 PD_whole_tree，若不计算它就无需加载。事实上，picante 包加载时默认同时加载 vegan
library(ggplot2)	#用于 ggplot2 作图
library(doBy)	#用于分组统计
library(ggalt)
library(RColorBrewer)
setwd("")

#1.稀释曲线
#定义函数
alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
  if (method == 'richness') result <- rowSums(x > 0)	#丰富度指数
  else if (method == 'chao1') result <- estimateR(x)[3, ]	#Chao1 指数
  else if (method == 'ace') result <- estimateR(x)[5, ]	#ACE 指数
  else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)	#Shannon 指数
  else if (method == 'simpson') result <- diversity(x, index = 'simpson')	#Gini-Simpson 指数
  else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)	#Pielou 均匀度
  else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)	#goods_coverage
  else if (method == 'pd' & !is.null(tree)) {	#PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- pd[ ,1]
    names(result) <- rownames(pd)
  }
  result
}

#根据抽样步长（step），统计每个稀释梯度下的 Alpha 多样性指数，结果返回至列表
alpha_curves <- function(x, step, method = 'richness', rare = NULL, tree = NULL, base = exp(1)) {
  x_nrow <- nrow(x)
  if (is.null(rare)) rare <- rowSums(x) else rare <- rep(rare, x_nrow)
  alpha_rare <- list()
  
  for (i in 1:x_nrow) {
    step_num <- seq(0, rare[i], step)
    if (max(step_num) < rare[i]) step_num <- c(step_num, rare[i])
    
    alpha_rare_i <- NULL
    for (step_num_n in step_num) alpha_rare_i <- c(alpha_rare_i, alpha_index(x = rrarefy(x[i, ], step_num_n), method = method, tree = tree, base = base))
    names(alpha_rare_i) <- step_num
    alpha_rare <- c(alpha_rare, list(alpha_rare_i))
  }
  
  names(alpha_rare) <- rownames(x)
  alpha_rare
}
#导入文件
otu<- read.csv("2.asv_absolute_clean1.csv",row.names = 1,header=T)
otu <- t(otu)
rownames(otu)<-otu_1$Group
otu_1<- read.csv("1.group.csv",row.names = 1,header=T)


richness_curves <- alpha_curves(otu, step = 2500, method = 'richness')

#获得 ggplot2 作图文件
plot_richness <- data.frame()
for (i in names(richness_curves)) {
  richness_curves_i <- (richness_curves[[i]])
  richness_curves_i <- data.frame(rare = names(richness_curves_i), alpha = richness_curves_i, sample = i, stringsAsFactors = FALSE)
  plot_richness <- rbind(plot_richness, richness_curves_i)
}

rownames(plot_richness) <- NULL
plot_richness$rare <- as.numeric(plot_richness$rare)
plot_richness$alpha <- as.numeric(plot_richness$alpha)
library(Rmisc)

#作图
#library(randomcoloR)
#palette <- randomColor(count = 60)  #随机生成60种颜色，其实里面有重复的
#palette <- distinctColorPalette(191)#差异明显的191种

plot_richness$sample <- factor(plot_richness$sample, levels = c('Health', 'Disease',"Bile","Stone"))
#无置信区间，为解决，参数未起作用？？？
p1<-ggplot() +
  stat_summary(data=plot_richness,
               aes(x=rare, y=alpha,color=sample,fill=sample),
               fun.data="mean_cl_boot",geom="ribbon",
               fun.args=list(conf.int=0.95),
               alpha=1,size=1.5) +
  labs(x = 'Number of sequences sampled', y = 'Rarefaction Measure:Richness', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'),
        axis.title.x = element_text(size = 15,face = "bold"),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.ticks.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 1),
        panel.border = element_rect(fill=NA,color="black", size=1),
        legend.text=element_text(size=10),
        legend.title=element_blank()) +
  geom_vline(xintercept = min(rowSums(otu)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 55000, 5000), labels = as.character(seq(0, 55000, 5000)))
p1+scale_color_nejm()+scale_fill_nejm()

#方法2：
#install.packages("iNEXT")
library(iNEXT)
otu_data<- read.csv("2.asv_absolute_clean1_group.csv",row.names = 1,header=T)
otu <- iNEXT(otu_data, q=0, datatype="abundance",endpoint=25000)
ggiNEXT(otu, type=1)

#2.α多样性指数
library(picante)
library(ggpubr)
alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[3, ]
  ACE <- est[5, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  result
}
otu<- read.csv("2.asv_absolute_clean1.csv",row.names = 1,header=T)
otu <- t(otu)
alpha_all <- alpha(otu, base = exp(1))
write.csv(alpha_all, 'alpha_16S_all.csv', quote = FALSE)
#导入数据
diff_alpha <- read.csv("alpha_16S_all.csv",row.names = 1,header=T)

# 分组添加levels，分组有时间顺序
diff_alpha$Group<- factor(diff_alpha$Group, levels = c('Health', 'Disease',"Bile","Stone"))


pdf("1.figureE.pdf",width = 1.2,height = 1.5) 
p<-ggplot(diff_alpha, aes(x=Group, y=ACE))+
  #stat_boxplot(geom="errorbar",position = position_dodge(0.9),width=0.2,aes(color=Group))+
  geom_violin(aes(fill = Group), position = position_dodge(0.9),alpha = 0.9 ,width = 1.0,linetype = "blank") +
  geom_boxplot(aes(fill=Group),position = position_dodge(0.9),lwd=0.1,width = 0.3,outlier.shape = NA)+
  #geom_jitter(aes(fill=Group),width =0.2,shape = 21,size=0.1)+
  #注意画图顺序，如果反过来会出现十字先stat再geom
  theme_classic()+
  xlab("")+labs(y="ACE index of ASV level")+
  theme(legend.position="none",
        panel.grid = element_line(color = 'gray', linetype =" blank", size = 0.1), 
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 7,face = "bold"),
        axis.title.y = element_text(size = 7,face = "bold"),
        axis.text.x = element_text(size = 6,angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        #panel.border = element_rect(fill=NA,color="black", size=1),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        legend.text=element_text(size=5),
        legend.title=element_blank())
p1<-p+scale_fill_manual(values = c("#A00002","#000000","#0850A1","#999999"))+
  scale_color_manual(values = c("#A00002","#000000","#0850A1","#999999"))+
  stat_compare_means(comparisons = list(c('Health', 'Disease'),
                                        c('Health', 'Bile'),
                                        c('Health', 'Stone'),
                                        c('Disease', 'Bile'),
                                        c('Disease', 'Stone'),
                                        c('Bile', 'Stone')),size=2)
p1
dev.off()  

pdf("1.figureD.pdf",width = 1.2,height = 1.5) 
p<-ggplot(diff_alpha, aes(x=Group, y=Chao1))+
  #stat_boxplot(geom="errorbar",position = position_dodge(0.9),width=0.2,aes(color=Group))+
  geom_violin(aes(fill = Group), position = position_dodge(0.9),alpha = 0.9 ,width = 1.0,linetype = "blank") +
  geom_boxplot(aes(fill=Group),position = position_dodge(0.9),lwd=0.1,width = 0.3,outlier.shape = NA)+
  #geom_jitter(aes(fill=Group),width =0.2,shape = 21,size=0.1)+
  #注意画图顺序，如果反过来会出现十字先stat再geom
  theme_classic()+
  xlab("")+labs(y="Chao1 index of ASV level")+
  theme(legend.position="none",
        panel.grid = element_line(color = 'gray', linetype =" blank", size = 0.1), 
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 7,face = "bold"),
        axis.title.y = element_text(size = 7,face = "bold"),
        axis.text.x = element_text(size = 6,angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        #panel.border = element_rect(fill=NA,color="black", size=1),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        legend.text=element_text(size=5),
        legend.title=element_blank())
p1<-p+scale_fill_manual(values = c("#A00002","#000000","#0850A1","#999999"))+
  scale_color_manual(values = c("#A00002","#000000","#0850A1","#999999"))+
  stat_compare_means(comparisons = list(c('Health', 'Disease'),
                                        c('Health', 'Bile'),
                                        c('Health', 'Stone'),
                                        c('Disease', 'Bile'),
                                        c('Disease', 'Stone'),
                                        c('Bile', 'Stone')),size=2)
p1
dev.off()  
#粪便配对两组间
diff_alpha <- read.csv("alpha_16S_all_t.csv",row.names = 1,header=T)
data4<- read.table("asv_fecal_meta_match.txt",header=T)
vars<-data4$OUT_ID
diff_alpha<-diff_alpha[vars]
diff_alpha<-as.data.frame(t(diff_alpha))
diff_alpha$Group<- factor(diff_alpha$Group, levels = c('Health', 'Disease'))
diff_alpha$Chao1<-as.numeric(diff_alpha$Chao1)
diff_alpha$Shannon<-as.numeric(diff_alpha$Shannon)
diff_alpha$Simpson<-as.numeric(diff_alpha$Simpson)
diff_alpha$ACE<-as.numeric(diff_alpha$ACE)
diff_alpha$Pielou<-as.numeric(diff_alpha$Pielou)
#write.csv(diff_alpha,"alpha_16S_match_240129.csv")
#先用各组箱线图查看数据分布
p<-ggplot(diff_alpha, aes(x=Group, y=ACE))+
  stat_boxplot(geom="errorbar",position = position_dodge(0.9),width=0.2)+
  geom_boxplot(aes(fill=Group),position = position_dodge(0.9),lwd=0.6)+#注意画图顺序，如果反过来会出现十字先stat再geom
  geom_jitter(aes(fill=Group),width =0.2,shape = 21,size=2)+
  scale_color_manual(values=c("black","black"))+
  theme_classic()+
  xlab("")+labs(y="ACE index of ASV level")+
  theme(legend.position = "none",
        panel.grid = element_line(color = 'gray', linetype = 0, size = 0.2), 
        panel.border = element_rect(fill=NA,colour = 'black',linewidth =2),
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 16,face = "bold"),
        axis.title.y = element_text(size = 16,face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.ticks.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 1))
p+scale_fill_d3()+
  scale_color_manual(values = c("#000000","#000000"))+
  stat_compare_means(comparisons = list(c('Health', 'Disease')),method ="wilcox.test",size=5)

pdf("sup2.figureA.pdf",width = 1,height = 1.2) 
p<-ggplot(diff_alpha, aes(x=Group, y=ACE))+
  geom_violin(aes(fill = Group), position = position_dodge(0.9),alpha = 0.9 ,width = 1.0,linetype = "blank") +
  geom_boxplot(aes(fill=Group),position = position_dodge(0.9),lwd=0.1,width = 0.3,outlier.shape = NA)+
  theme_classic()+
  xlab("")+labs(y="ACE index of ASV level")+
  theme(legend.position="none",
        panel.grid = element_line(color = 'gray', linetype =" blank", size = 0.1), 
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 8,face = "bold"),
        axis.title.y = element_text(size = 8,face = "bold"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        legend.text=element_text(size=5),
        legend.title=element_blank())
p1<-p+scale_fill_manual(values = c("#1F77B4","#739884"))+
  scale_color_manual(values = c("#1F77B4","#739884"))+
  stat_compare_means(comparisons = list(c('Health', 'Disease')),method ="wilcox.test",size=2)
p1
dev.off()  


