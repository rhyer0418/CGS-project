setwd("")
#install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("DESeq2")
library(DESeq2)

#准备文件
#OTU丰度表格，otu.txt；(每一列为一个样本，每一行为一种OTU，交叉区域为每种OTU在各样本中的丰度,DESeq2计算时只能识别整数，不识别小数，所以请不要使用相对丰度表格)
#样本分组信息表格，即metadata.txt。(行为样本名称，第二列为各样本对应的分组信息)
otu<- read.csv("3.4.genus_abs_match_clean.csv",row.names=1,header=T)
group<- read.csv("1.group_match.csv",row.names=1,header=T)
group<-group[,1:2]

dds <- DESeqDataSetFromMatrix(countData=otu, colData=group, design=~Group) #构建 DESeqDataSet 对象  
dds <- DESeq(dds) #差异分析
suppressMessages(dds)

res <- results(dds, contrast=c('Group', 'Health', 'Disease'))#提取分析结果
res = res[order(res$pvalue),]
res #查看结果
summary(res)  #简要统计结果
table(res$padj<0.05) #查看fdr校正后的P<0.05的个数

#提取差异OTU并进行注释
#获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异OTU

diff_OTU_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
#或
# diff_OTU_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
dim(diff_OTU_deseq2)
head(diff_OTU_deseq2)
write.csv(diff_OTU_deseq2,file= "DEG_yes_vs_no_match.csv")

#合并结果
resdata <-merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
write.table(resdata,file= "DESeq2_Group_volvano.txt",sep="\t",quote=F,row.names=F)


#绘图 
library(ggplot2)
library(ggrepel)
resdata$color <- ifelse(resdata$padj<0.05 & abs(resdata$log2FoldChange)>= 1,ifelse(resdata$log2FoldChange > 1,'red','blue'),'gray')
color <- c(red = "#81A290",gray = "#C5C9CC",blue = "#3A85BC")

pdf("sup2.figureC.pdf",width = 2,height = 2.5) 

p <- ggplot(resdata, aes(log2FoldChange, -log10(padj), col = color)) +
  geom_point(size=1.2) +
  theme_test() +
  scale_color_manual(values = color) +
  geom_text_repel(
    data = resdata[resdata$padj<0.05&abs(resdata$log2FoldChange)>1,],
    aes(label = Row.names),
    size = 1,
    segment.color = "black", show.legend = FALSE )+
  labs(x="log2(FoldChange)",y="-log10(q-value)") +
  scale_y_continuous(breaks=seq(0,12,1))+
  scale_x_continuous(breaks=seq(-12,12,2))+
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.3) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.3) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 8,face = "bold"),
        axis.title.y = element_text(size = 8,face = "bold"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        panel.border = element_rect(fill=NA,color="black", size=1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 1),
        legend.text=element_text(size=5),
        legend.title=element_blank())
p
dev.off()



