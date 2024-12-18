setwd("")
data<- read.csv("5.heatmap_genus_4group_top20.csv",row.names = 1,header=T)

library(ggcorrplot)
library(ggthemes)
library(pheatmap)
library(psych)
library(gridExtra)
library(grid)

drows = dist(data, method = "euclidean")
dcols = dist(t(data), method = "euclidean")

annotation_col = read.csv("1.group_name.csv",row.names=1,header = T)
annotation_col$Group<-as.factor(annotation_col$Group)

ann_colors = list(
  Group =c("Health"="#FED854","Disease"="#5A88FA","Bile"="#D7B78E","Stone"="#FE825C"))
head(ann_colors)

pdf("1.figureF.pdf",width = 4.5,height = 3) 
heatmap<-pheatmap(data,cluster_cols=T,cluster_row=F,scale = "column",clustering_method="ward.D",
                  annotation_col = annotation_col,annotation_colors = ann_colors,
         clustering_distance_rows = drows, clustering_distance_cols = dcols,
         cellwidth=1,cellheight=6.5,fontsize_number=1,number_color="white",
         fontsize_row=8,fontsize_col=2.5,
         color=colorRampPalette(colors = c("#0750A1","white","#386092"))(100),
         border="white" ,
         legend=FALSE
  )
dev.off()  

pdf("1.figureF1.pdf",width = 4.5,height = 3) 
heatmap<-pheatmap(data,cluster_cols=T,cluster_row=F,scale = "column",clustering_method="ward.D",
                  annotation_col = annotation_col,annotation_colors = ann_colors,
                  clustering_distance_rows = drows, clustering_distance_cols = dcols,
                  cellwidth=1,cellheight=6.5,fontsize_number=1,number_color="white",
                  fontsize_row=8,fontsize_col=2.5,
                  color=colorRampPalette(colors = c("#0750A1","white","#386092"))(100),
                  border="white" ,
                  #legend=FALSE
)
dev.off()  