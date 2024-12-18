#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")

library(limma)
library(reshape2)
library(ggpubr)
geneFile="gene_sup fig 4.txt"         
cliFile="1.group_m.txt"             

setwd("")

genus<-read.csv("kegg_ko.csv",header = T)
group<- read.csv("1.group_match.csv",row.names=1,header=T)
gene=read.table(geneFile, header=T, sep="\t", check.names=F)
name<-rownames(group)
genus2<-genus[name]
genus2$Species<-genus$KO

#write.csv(genus2,"gene_KO_match.csv",row.names = FALSE)

##选择相关基因
genus1<-as.data.frame(t(genus))
colnames(genus1)<-genus1[1,]
genus1<-genus1[-1,]
sameGene=intersect(as.vector(gene[,1]), colnames(genus1))
row.names(genus2)<-genus2$Species
data=genus2[sameGene,]

#全部基因
#data=genus2
#rownames(data)=data$Species
#supplementary figure 4

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
expdata=data
id=substr(colnames(expdata),1,12)
sameSample=intersect(id, row.names(cli))
expdata<-expdata[,-61]
colnames(expdata)= substr(colnames(expdata), 1, 12)
expdata<-as.data.frame(t(expdata))
expdata=expdata[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data2=cbind(expdata, cli)


data3=melt(data2, id.vars=c("Group"))
colnames(data3)=c("Type", "Gene", "Expression")
data3$Expression<-10000*data3$Expression

pdf(file="Supplementary Figure 4.pdf", width=9, height=7)

p=ggboxplot(data3, x="Gene", y="Expression", color = "Type", 
	     ylab="Relative abundance (per 104 reads)",
	     xlab="",
	     legend.title="Group",
	     palette = c("#0072B5","#749884"),
	     width=0.7)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

print(p1)
dev.off()




