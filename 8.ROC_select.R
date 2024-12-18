setwd("")

#结果合并需要的包
library(plyr)
library(pROC)
library(ROCR) 

diff<-read.csv("diff_sig_genus_match_unadjusted_roc.csv",header = T,row.names = 1)

#roc=roc(diff$group,diff$Catenibacterium)

outTab <- data.frame()
for(i in colnames(diff[,2:ncol(diff)])){
  roc=roc(diff$group,diff[,i]) 
  if(roc$auc>0.7){
    outTab=rbind(outTab,cbind(genus=i,AUC=roc$auc)) 
    }
  }
#write.csv(outTab,"roc_16s_0.7_12_match.csv",row.names = FALSE)


diff$group<-as.factor(diff$group)
fit1=glm(data=diff,group~Terrisporobacter+Erysipelatoclostridium+	
          Lachnospiraceae_UCG.010+	
          Ruminococcus_gnavus_group+
        Lachnospiraceae_UCG.001+
          Clostridium_innocuum_group+
        Senegalimassilia+	
          Lachnospiraceae_UCG.004+
        #UCG.010+
        Roseburia+Sutterella+	Megamonas,
         family = binomial()) 

prob1 <- predict(fit1, diff, type = "response")
rocboj1<-roc(diff$group,prob1,levels=c("Disease","Health"))
roc_result <- coords(rocboj1, "best")
auc1<-round(auc(diff$group,prob1),4)

plot(rocboj1,
     legacy.axes = TRUE,
    main="ROC曲线最佳阈值点",
    thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
     print.thres="best") # 在roc曲线上显示最佳阈值点

p=ggroc(rocboj1,color="#00A1D5FF", linetype=1,size=1,alpha=1,legacy.axes = T) +
  geom_abline(intercept = 0,slope = 1,color="grey",size=1,linetype=2)+
  labs(x="False Positive Rate (1 - Specificity)", y="True Positive Rate (Sensitivity)") +
  annotate("text",x=.75,y=.25,label=paste("AUC=",auc1),size=5,family="serif")+
  coord_cartesian(xlim = c(0,1),ylim = c(0,1))+
  theme_bw() +
  theme(panel.grid = element_line(color = 'gray', linetype = 0, size = 0.2), 
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 15,face = "bold"),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.border = element_rect(fill=NA,color="black", size=2),
        axis.ticks.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 1),
        legend.text=element_text(size=10),
        legend.title=element_blank())
p1=p+scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00),expand=c(0,0))+expand_limits(y=c(0,1))
p1

setwd("")
diff<-read.csv("diff_sig_species_match_unadjusted_lefse_boxplot_roc_select.csv",header = T,row.names = 1)
outTab <- data.frame()
for(i in colnames(diff[,2:ncol(diff)])){
  roc=roc(diff$group,diff[,i]) 
  if(roc$auc>0.68){
    outTab=rbind(outTab,cbind(genus=i,AUC=roc$auc)) 
  }
}

diff$group<-as.factor(diff$group)
fit2=glm(data=diff,group~	
          Thermodesulfovibrio.thiophilus+
        Leucothrix.mucor+
        Clostridium.perfringens+Gordonibacter.massiliensis+Gordonibacter.sp..An230+	
          Lachnospira.multipara+Desulfitobacterium.hafniense+Mycobacterium.avium,family = binomial()) 

fit2=glm(data=diff,group~	
           Faecalibacterium.prausnitzii,family = binomial()) 
fit2=glm(data=diff,group~	
           X.Ruminococcus..gnavus,family = binomial()) 
fit2=glm(data=diff,group~	
           X.Clostridium..glycyrrhizinilyticum,family = binomial()) 


prob2 <- predict(fit2, diff, type = "response")
rocboj2<-roc(diff$group,prob2)
roc_result <- coords(rocboj2, "best")
auc2<-round(auc(diff$group,prob2),4)

pdf("Fp_ROC.pdf",
    width = 2.1,height = 2)
p=ggroc(rocboj2,color="#0072B5", linetype=1,size=1,alpha=1,legacy.axes = T) +
  geom_abline(intercept = 0,slope = 1,color="grey",size=1,linetype=2)+
  labs(x="False Positive Rate (1 - Specificity)", y="True Positive Rate (Sensitivity)") +
  annotate("text",x=.75,y=.25,label=paste("AUC=",auc2),size=2.5,family="serif")+
  coord_cartesian(xlim = c(0,1),ylim = c(0,1))+
  theme_bw() +
  theme(panel.grid = element_line(color = 'gray', linetype = 0, size = 0.2), 
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 8.5,face = "bold"),
        axis.title.y = element_text(size =8.5,face = "bold"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        panel.border = element_rect(fill=NA,color="black", size=1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 1),
        legend.text=element_text(size=10),
        legend.title=element_blank())
p1=p+scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00),expand=c(0,0))+expand_limits(y=c(0,1))
p1
dev.off()

pdf("Rg_ROC.pdf",
    width = 2.1,height = 2)
p=ggroc(rocboj2,color="#E38A04", linetype=1,size=1,alpha=1,legacy.axes = T) +
  geom_abline(intercept = 0,slope = 1,color="grey",size=1,linetype=2)+
  labs(x="False Positive Rate (1 - Specificity)", y="True Positive Rate (Sensitivity)") +
  annotate("text",x=.75,y=.25,label=paste("AUC=",auc2),size=2.5,family="serif")+
  coord_cartesian(xlim = c(0,1),ylim = c(0,1))+
  theme_bw() +
  theme(panel.grid = element_line(color = 'gray', linetype = 0, size = 0.2), 
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 8.5,face = "bold"),
        axis.title.y = element_text(size =8.5,face = "bold"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        panel.border = element_rect(fill=NA,color="black", size=1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 1),
        legend.text=element_text(size=10),
        legend.title=element_blank())
p1=p+scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00),expand=c(0,0))+expand_limits(y=c(0,1))
p1
dev.off()

pdf("Cg_ROC.pdf",
    width = 2.1,height = 2)
p=ggroc(rocboj2,color="#208548", linetype=1,size=1,alpha=1,legacy.axes = T) +
  geom_abline(intercept = 0,slope = 1,color="grey",size=1,linetype=2)+
  labs(x="False Positive Rate (1 - Specificity)", y="True Positive Rate (Sensitivity)") +
  annotate("text",x=.75,y=.25,label=paste("AUC=",auc2),size=2.5,family="serif")+
  coord_cartesian(xlim = c(0,1),ylim = c(0,1))+
  theme_bw() +
  theme(panel.grid = element_line(color = 'gray', linetype = 0, size = 0.2), 
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 8.5,face = "bold"),
        axis.title.y = element_text(size =8.5,face = "bold"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        panel.border = element_rect(fill=NA,color="black", size=1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 1),
        legend.text=element_text(size=10),
        legend.title=element_blank())
p1=p+scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00),expand=c(0,0))+expand_limits(y=c(0,1))
p1
dev.off()


plot(rocboj2,
   legacy.axes = TRUE,
   main="ROC曲线最佳阈值点",
   thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
   print.thres="best") # 在roc曲线上显示最佳阈值点


setwd("")
diff<-read.csv("roc_metabo_match.csv",header = T,row.names = 1,check.names=FALSE)
diff2<-read.csv("diff_metabo2.csv",header = T,row.names = 1,check.names=FALSE)


outTab <- data.frame()
for(i in colnames(diff[,2:ncol(diff)])){
  roc=roc(diff$group,diff[,i]) 
  if(roc$auc>0.68){
    outTab=rbind(outTab,cbind(genus=i,AUC=roc$auc)) 
  }
}
#write.csv(outTab,"roc_metabolites_0.68_33_match.csv",row.names = FALSE)

diff$group<-as.factor(diff$group)
fit3=glm(data=diff,group~	
          `1-Methylhistidine`+
          #`Hypoglycin A`+
           `Ureidosuccinic acid`+
           `L-Asparagine`+
          `N-Acetylphenylalanine`+
           #`2-(Formylamino)Benzoic Acid`+
           #Norfenefrine+	
          `4-Aminobenzoic acid`+
           `N-Glycolylneuraminic acid`+
           `4-Coumaric acid`+
          `3-Coumaric acid`+
          `3-Hydroxyvaleric acid`+
           Cholecalciferol+
           `delta-Tocopherol`+
           Neopterin+
           #Biopterin+
         `Nicotinic acid`
          #`Uridine 5'-diphospho-D-glucose`
           ,family = binomial()) 

prob3 <- predict(fit3, diff, type = "response")
rocboj3<-roc(diff$group,prob3,levels=c("Disease","Health"))
roc_result <- coords(rocboj3, "best")
auc3<-round(auc(diff$group,prob3),4)
#levels=c("Disease","Health"),

pdf("3.fig3E.pdf",
    width =2.3,height = 2.3)

p=ggroc(list(Genus=rocboj1,Species=rocboj2,Metabolites=rocboj3)
  , linetype=1,size=1,alpha=1,legacy.axes = T) +
  geom_abline(intercept = 0,slope = 1,color="grey",size=1,linetype=2)+
  labs(x="False Positive Rate (1 - Specificity)", y="True Positive Rate (Sensitivity)") +
  annotate("text",x=.75,y=.25,label=paste("AUC=",auc1),size=5,family="serif")+
  annotate("text",x=.75,y=.35,label=paste("AUC=",auc2),size=5,family="serif")+
  annotate("text",x=.75,y=.45,label=paste("AUC=",auc3),size=5,family="serif")+
  coord_cartesian(xlim = c(0,1),ylim = c(0,1))+
  theme_bw() +
  theme(panel.grid = element_line(color = 'gray', linetype = 0, size = 0.2), 
        panel.background = element_blank(), 
        axis.title.x = element_text(size = 9,face = "bold"),
        axis.title.y = element_text(size = 9,face = "bold"),
        axis.text.x = element_text(size = 8.5),
        axis.text.y = element_text(size = 8.5),
        panel.border = element_rect(fill=NA,color="black", size=2),
        axis.ticks.x=element_line(color="black",size=0.5,lineend = 1),
        axis.ticks.y=element_line(color="black",size=0.5,lineend = 1),
        legend.text=element_text(size=6),legend.position ="NA",
        legend.title=element_blank())+
  scale_color_manual(values = c("#334A7D",  "#E38A06","#B82630"))
p1=p+scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00),expand=c(0,0))+expand_limits(y=c(0,1))
p1

dev.off()

roc=roc(diff$group,diff$`DL-α-Methoxyphenylacetic acid`)

plot(rocboj3,
     legacy.axes = TRUE,
     main="ROC曲线最佳阈值点",
     thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
     print.thres="best") # 在roc曲线上显示最佳阈值点


