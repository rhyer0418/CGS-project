setwd("")
df <- read.csv("very_easy_input1.csv")
head(df)

geneSpecial <- read.table("gene_special.txt", header = T)
geneCol <- geneSpecial$Type
names(geneCol) <- geneSpecial$Gene
geneCol
#加载函数：`gather_graph_node()` 和 `gather_graph_edge()`
nodes <- gather_graph_node(df, index = c("term", "genes"), value = "logFC", root="all")
edges <- gather_graph_edge(df, index = c("term", "genes"), root = "all")
nodes <- nodes %>% mutate_at(c("node.level","node.branch"),as.character)
head(nodes, 10)
head(edges, 10)

# 把要突出显示的基因类型信息加到nodes里
nodes$color <- "normal"
nodes[nodes$node.short_name %in% geneSpecial$Gene,]$color <- geneCol[nodes[nodes$node.short_name %in% geneSpecial$Gene,]$node.short_name]
nodes[nodes$node.short_name %in% geneSpecial$Gene,]
nodes$color <- factor(nodes$color, levels = unique(nodes$color))

library("tidygraph")
library(ggraph)
# 有了节点和边的数据，使用 `tbl_graph()` 便可以得到一个图。
graph <- tbl_graph(nodes, edges)

# 用 `ggraph` 出图
gc <- ggraph(graph, layout = 'dendrogram', circular = TRUE) + 
  # 使用 filter 参数去掉 root（前面设置为"all"）节点及与其相连的边
  geom_edge_diagonal(aes(color = node1.node.branch,
                         filter=node1.node.level!="all"), 
                     alpha = 1,edge_width=1) + 
  geom_node_point(aes(size = node.size, 
                      color = node.branch,
                      filter=node.level!="all"), alpha = 0.9) + 
  scale_size(range = c(0.5,8)) + #做均一化处理，让点的大小介于range之间
  theme(legend.position = "none") + #不画图例
  
  # 点和边的配色
  # 如果要改变边的配色，需要同时给边和点变色，否则会对应不上
  #scale_edge_color_brewer(palette = "Set1") + #用?scale_color_brewer查看更多配色方案
  #scale_color_brewer(palette = "Set1") +
  scale_edge_color_manual(values = c("#FF7F00","#E18727FF","#7E6148B2","#EE4C97FF",
                                     "#FFDC91FF","#F39B7FB2","#BC3C29FF","#20854EFF","#7876B1FF",
                                     "#80796BFF","#984EA3","#0072B5FF","#FB8072","#61C3ED","#BEBADA",  
                                     #"#D95F02","#80B1D3","#BEBADA","#8DD3C7","#E41A1C",
                                     "#a06fda","#dd6bbb","#334c22","#5b83db","#6F99ADFF")) + #自定义颜色
  scale_color_manual(values = c("#FF7F00","#E18727FF","#7E6148B2","#EE4C97FF",
                                "#FFDC91FF","#F39B7FB2","#BC3C29FF","#20854EFF","#7876B1FF",
                                "#80796BFF","#984EA3","#0072B5FF","#FB8072","#61C3ED","#BEBADA",  
                                #"#D95F02","#80B1D3","#BEBADA","#8DD3C7","#E41A1C",
                                "#a06fda","#dd6bbb","#334c22","#5b83db","#6F99ADFF")) + #自定义颜色
  
  # 添加周围注释文字，此处是基因名gene
  geom_node_text(
    aes(
      x = 1.05 * x, #控制字跟点的距离
      y = 1.05 * y, #控制字跟点的距离
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf
    ),
    color="black", #统一为黑色字
    size = 1.5, hjust = 'outward') + theme(legend.position = "none")+
  
  # 添加内环文字，此处是通路名term
  #geom_node_text(
   # aes(label=node.short_name,
    #    filter = !leaf & (node.level != "all")
    #),
    #color="black", #统一为黑色字
    #fontface="bold",
    #size=3,
    #family="sans"
  #) + 
  # 添加内环文字
  theme(legend.position = "bottom",legend.title = element_blank(),legend.spacing = unit(0.1, "lines"),
        legend.key = element_blank(),legend.text = element_text(size = 1),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        legend.background = element_blank()) + # 设置图例位置为底部
  theme(panel.background = element_rect(fill = NA)) +
  guides(color = guide_legend(override.aes = list(size = 0.1), keyheight = unit(0.1, "cm"), keywidth = unit(0.1, "cm")))
  #coord_cartesian(xlim=c(-1.3,1.3),ylim = c(-1.3,1.3)) #扩大坐标系
gc

ggsave("ccgraph_color_new.pdf", width = 3, height =4)
ggsave("ccgraph_color_new1.pdf", width = 6, height =6)
#---在上面的图形中，线条的颜色由 
#` geom_edge_diagonal(aes(color = node1.node.branch)) ` 指定。 
#`node1.node.branch` 指的是出发点（`node1`）的 `node.branch` 属性。
 # 如果要改变线条颜色，可以修改 nodes 表，添加一个属性（如 `color`） ，
   #然后在 `geom_edge_diagonal()` 中将其映射到 `color` 上即可。
### 按例文配色，然后后期加背景色
nodes<-write_csv(nodes,"nodes.csv")
nodes<-read.csv("nodes.csv",header = T)
nodes$color<-factor(nodes$color,levels = unique(nodes$color))

graph <- tbl_graph(nodes, edges)
gc1 <- ggraph(graph, layout = 'dendrogram', circular = TRUE) + 
  #画连线
  geom_edge_diagonal(aes(color = node2.color,
                         filter=node1.node.level!="all"), 
                     alpha = 1, #透明度
                     edge_width=2) + #连线的粗细
  scale_edge_color_manual(values = c("#61C3ED","#F39B7FB2","#EE4C97FF","#7E6148B2","#E18727FF",
                                     "#FFDC91FF","#0072B5FF","#BC3C29FF","#20854EFF","#7876B1FF","#6F99ADFF")) + #自定义颜色
  
  #画点
  geom_node_point(aes(size = node.size,
                      filter=node.level!="all"), 
                  #alpha = 1/3,
                  color = "#61C3ED") + #统一为淡蓝色
  scale_size(range = c(0.5,40)) + #做均一化处理，让点的大小介于range之间
  theme(legend.position = "right") + #不画图例
  
  # 添加周围注释文字，此处是基因名gene
  geom_node_text(
    aes(
      x = 1.05 * x, #控制字跟点的距离
      y = 1.05 * y, #控制字跟点的距离
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf
    ),
    color="black", #统一为黑色字
    size = 2, hjust = 'outward') +
  
  # 添加内环文字，此处是通路名term
  geom_node_text(
    aes(label=node.short_name,
        filter = !leaf & (node.level != "all")
    ),
    color="black", #统一为黑色字
    fontface="bold",
    size=4,
    family="sans"
  ) + 
  theme(panel.background = element_rect(fill = NA)) + #背景透明色
  coord_cartesian(xlim=c(-1.3,1.3),ylim = c(-1.3,1.3)) #扩大坐标系

gc1

#保存
ggsave("ccgraph1.pdf",width = 17,height = 14)
#or
library(export)
graph2ppt(gc1, file="ccgraph.pptx",width = 14,height = 14)


  