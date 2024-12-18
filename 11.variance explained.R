setwd("")
library(tidyverse)
diff<-read.csv("diff_sig_species_match_unadjusted_lefse_boxplot_roc_s_test.csv",header = T)
diff2<-as.data.frame(t(diff))
colnames(diff2)<-diff$Species
diff2<-diff2[-1,]
diff2 <- mutate_all(diff2, as.numeric)
diff2$Group<-group$Group
diff2$Group <- as.factor(diff2$Group)
diff2$Group <- ifelse(diff2$Group == "Health", 0, 1)
rsquared_p_df <- data.frame(
  Microbe = colnames(diff2)[1:(ncol(diff2)-24)], # 使用原始数据框的列名作为菌种名称
  R_squared = numeric(length = (ncol(diff2)-24)),
  P_value = numeric(length = (ncol(diff2)-24))
)

# 循环计算 R squared 和 p 值，并存储到数据框中
for (i in 1:(ncol(diff2)-24)) {
  lm_model <- lm(TC ~ diff2[,i], data = diff2, family = binomial)
  rsquared_p_df$R_squared[i] <- summary(lm_model)$r.squared
  rsquared_p_df$P_value[i] <- summary(lm_model)$coefficients[, 4] # 提取第四列，即 p 值这个数值有问题
}

# 计算 univariate linear regression 的 R squared 值  ，只能是数值变量
var<-c('TC','TG','TBIL','HDLC','LDLC','ALT','AST')
rsquared_p_list <- list()
# 循环计算 R squared 和 p 值，并存储到列表中
for (var_name in var) {  # 循环遍历 var 中的每个因变量
  temp_df <- data.frame(Microbe = character(),
                        Var_Name = character(),
                        R_squared = numeric(),
                        P_value = numeric(),
                        stringsAsFactors = FALSE)
  for (i in 1:(ncol(diff2)-24)) {
    lm_model <- lm(as.formula(paste0(var_name, "~",  "`",colnames(diff2)[i],"`+NEUT")), data = diff2, family = binomial)
    temp_df <- rbind(temp_df, data.frame(
      Microbe = colnames(diff2)[i],
      Var_Name = var_name,
      R_squared = summary(lm_model)$r.squared,
      P_value = summary(lm_model)$coefficients[, 4]
    ))
  }
  rsquared_p_list[[var_name]] <- temp_df[-1, ]  # 去除第一行空白行
}
rsquared_p_df <- do.call(rbind, rsquared_p_list)
significant_microbes <- rsquared_p_df  %>%
  filter(0.0001< P_value & P_value< 0.05) # 这里根据需要调整显著性水

write.csv(rsquared_p_df,"species_test_lm_variance explained.csv")
significant_microbes<-read.csv("species_test_lm_variance explained.csv",header = T)
# 假设显著性水平为 0.05
significant_microbes <- significant_microbes %>%
  filter(P_value < 0.05) # 这里根据需要调整显著性水平
significant_microbes$Var_Name<-factor(significant_microbes$Var_Name,levels=c('HDLC','LDLC','ALT','AST','TC','TG','TBIL'))

#pdf("2.figure2E.pdf",width = 10,height =5) 
ggplot(significant_microbes, aes(x = 100*R_squared, y = reorder(Microbe, R_squared))) +
  geom_bar(stat = "identity", aes(fill = R_squared), width = 0.5) +
  facet_wrap(~ Var_Name, ncol = 7) +
  scale_fill_gradient(low = "#4AACC7", high = "#0072B5",guide = "none") + # 设置渐变色范围
  labs(x = "Variance explained (%)", y = "Microbe") +
  theme_classic()+
  theme(legend.key = element_blank(),
        panel.grid = element_line(color = 'gray', linetype = 0, size = 0.2), 
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 10,face = "bold"),
        #axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.ticks.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 1),
        axis.title.y=element_blank()) + # 调整y轴标签字体大小
  ggtitle("")
dev.off()

#ROC分类
library(pROC)
# 初始化一个空的数据框
auc_df <- data.frame(
  Microbe = character(),
  AUC = numeric(),
  stringsAsFactors = FALSE
)
# 对每个微生物菌种进行逻辑回归建模，并计算AUC值
for (i in 1:(ncol(diff2) - 1)) {
  glm_model <- glm(Group ~ diff2[,i], data = diff2, family = binomial)
  roc_obj <- roc(diff2$Group, fitted(glm_model))
  
  # 将菌种名称和AUC值添加到数据框中
  auc_df <- rbind(auc_df, data.frame(Microbe = colnames(diff2)[i], AUC = auc(roc_obj)))
}
print(auc_df)



setwd("/")
library(tidyverse)
diff<-read.csv("diff_sig_species_all_unadjusted_fdr_test.csv",header = T)
diff2<-as.data.frame(t(diff))
colnames(diff2)<-diff$Species
diff2<-diff2[-1,]
data<-diff2
data <- mutate_all(data, as.numeric)
data$Group<-diff2$Group
data$Group <- as.factor(data$Group)
data$Group <- ifelse(data$Group == "Health", 0, 1)

var<-c('TC','TG','TBIL','HDLC','LDLC','ALT','AST')
rsquared_p_list <- list()
# 循环计算 R squared 和 p 值，并存储到列表中
for (var_name in var) {  # 循环遍历 var 中的每个因变量
  temp_df <- data.frame(Microbe = character(),
                        Var_Name = character(),
                        R_squared = numeric(),
                        P_value = numeric(),
                        stringsAsFactors = FALSE)
  for (i in 1:(ncol(data)-23)) {
    #lm_model <- lm(as.formula(paste0(var_name, "~",  "`",colnames(data)[i],"`+NEUT")), data = data, family = binomial)
    lm_model <- lm(as.formula(paste0(var_name, "~",  "`",colnames(data)[i],"`")), data = data, family = binomial)
    temp_df <- rbind(temp_df, data.frame(
      Microbe = colnames(data)[i],
      Var_Name = var_name,
      R_squared = summary(lm_model)$r.squared,
      P_value = summary(lm_model)$coefficients[, 4]
    ))
  }
  rsquared_p_list[[var_name]] <- temp_df[-1, ]  # 去除第一行空白行
}
rsquared_p_df <- do.call(rbind, rsquared_p_list)
significant_microbes <- rsquared_p_df  %>%
  filter(P_value< 0.05) # 这里根据需要调整显著性水

#write.csv(significant_microbes,"species_test_lm_variance explained_all.csv")
significant_microbes<-read.csv("species_test_lm_variance explained_all.csv",header = T)
significant_microbes$Var_Name<-factor(significant_microbes$Var_Name,levels=c('TC','TG','ALT','AST','HDLC','LDLC','TBIL'))

pdf("2.figure2E.pdf",width =9,height =4.5) 
ggplot(significant_microbes, aes(x = 100*R_squared, y = reorder(Microbe, R_squared))) +
  geom_bar(stat = "identity", aes(fill = R_squared), width = 0.5) +
  facet_wrap(~ Var_Name, ncol = 7) +
  scale_fill_gradient(low = "#4AACC7", high = "#0072B5",guide = "none") + # 设置渐变色范围
  labs(x = "Variance explained (%)", y = "Microbe") +
  theme_classic()+
  theme(legend.key = element_blank(),
        panel.grid = element_line(color = 'gray', linetype = 0, size = 0.2), 
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 10,face = "bold"),
        #axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 12.5),
        axis.ticks.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 1),
        axis.title.y=element_blank()) + # 调整y轴标签字体大小
  ggtitle("")
dev.off()

