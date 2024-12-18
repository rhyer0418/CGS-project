setwd("workspace")
data<-read.csv("1.asv_absolute.csv",header = T)
d2<-data[,2:192]
d2[d2>0]<-1

health<-d2[,1:53]
health<-cbind(data,health)
disease<-d2[,54:119]
disease<-cbind(data,disease)
bile<-d2[,120:166]
bile<-cbind(data,bile)
stone<-d2[,167:191]
stone<-cbind(data,stone)

health$sum<-rowSums(health[,194:246],na.rm = T)
disease$sum<-rowSums(disease[,194:259],na.rm = T)
bile$sum<-rowSums(bile[,194:240],na.rm = T)
stone$sum<-rowSums(stone[,194:218],na.rm = T)

health<- subset(health, sum>5.3)
disease<- subset(disease, sum>6.6)
bile<-subset(bile, sum>4.7)
stone<-subset(stone, sum>2.5)

a=health$锘縊TU_ID
b=disease$锘縊TU_ID
c=bile$锘縊TU_ID
d=stone$锘縊TU_ID
a<-as.data.frame(a)
b<-as.data.frame(b)
c<-as.data.frame(c)
d<-as.data.frame(d)
#write.csv(a, file= "a.csv",row.names = T ,quote = F)
#write.csv(b, file= "b.csv",row.names = T ,quote = F)
#write.csv(c, file= "c.csv",row.names = T ,quote = F)
#write.csv(d, file= "d.csv",row.names = T ,quote = F)
colnames(a)<-"ID"
colnames(b)<-"ID"
colnames(c)<-"ID"
colnames(d)<-"ID"

e=rbind(a,b,c,d)
e=unique(e)
vars=e$ID

d1<-as.data.frame(t(data))
d1$sample_ID<-rownames(d1)
colnames(d1)<-data$锘縊TU_ID
d1<-d1[2:193,]

clean<-d1[vars]
write.csv(clean,"2.asv_absolute_clean.csv")

d4<-as.data.frame(t(clean))
d4<-d4[,1:191]
d4<-as.data.frame(t(d4))
d4=as.data.frame(lapply(d4,as.numeric))
ASV_relative <- d4 / rowSums(d4)*100

ASV_relative<-as.data.frame(t(ASV_relative))
colnames(ASV_relative)<-colnames(d2)
write.csv(ASV_relative, file= "3.ASV_relative.csv",row.names = T ,quote = F)


#分成门纲目科属种
d5<-as.data.frame(t(data))
d5$sample_ID<-rownames(d5)
colnames(d5)<-data$锘縊TU_ID
tax<-d5[vars]
tax<-tax[2:193,]
tax<-as.data.frame(t(tax))
library(tidyverse)
tax<-separate(data = tax, col = Taxonomy, into = c("kingdom","phylum", "class","order","family","genus","species"), sep = ";")
tax$kingdom<-as.factor(tax$kingdom)
tax$phylum<-as.factor(tax$phylum)
tax$genus<-as.factor(tax$genus)

tax1<-tax[,1:191]
taxonomy=as.data.frame(lapply(tax1,as.numeric))
taxonomy$phylum<-tax$phylum
#taxonomy$class<-tax$class
#taxonomy$order<-tax$order
#taxonomy$family<-tax$family
taxonomy$genus<-tax$genus
#taxonomy$species<-tax$species

phylum<-aggregate(taxonomy[,1:191],list(taxonomy$phylum),sum)
genus<-aggregate(taxonomy[,1:191],list(taxonomy$genus),sum)

#write.csv(phylum, file= "3.1.phylum_abun.csv",row.names = T ,quote = F)
#write.csv(genus, file= "3.2.genus_abun.csv",row.names = T ,quote = F)

genus<-read.csv("3.2.genus_abun.csv",header = T)
phylum<-read.csv("3.1.phylum_abun.csv",header = T)
p1<-phylum[,2:192]
p1<-as.data.frame(t(p1))
p1=as.data.frame(lapply(p1,as.numeric))
phylum_relative <- p1 / rowSums(p1)*100

phylum_relative<-as.data.frame(t(phylum_relative))
colnames(phylum_relative)<-colnames(phylum[,2:192])
rownames(phylum_relative)<-phylum$phylum
write.csv(phylum_relative, file= "3.3.phylum_relative.csv",row.names = T ,quote = F)
#----------------
p1<-genus[,2:192]
p1<-as.data.frame(t(p1))
p1=as.data.frame(lapply(p1,as.numeric))
genus_relative <- p1 / rowSums(p1)*100

genus_relative<-as.data.frame(t(genus_relative))
colnames(genus_relative)<-colnames(genus[,2:192])
rownames(genus_relative)<-genus$Genus
write.csv(genus_relative, file= "3.4.genus_relative.csv",row.names = T ,quote = F)

ph_m<-phylum_relative[,1:119]
ph_m<-ph_m[vars]

ge_m<-genus_relative[,1:119]
ge_m<-ge_m[vars]

write.csv(ge_m, file= "3.4.genus_rel_match.csv",row.names = T ,quote = F)
write.csv(ph_m, file= "3.4.phylum_rel_match.csv",row.names = T ,quote = F)

g_t<-as.data.frame(t(genus_relative))


#write.csv(genus, file= "3.4.genus_abs_match.csv",row.names = T ,quote = F)
otu<- read.csv("3.4.genus_abs_match.csv",header=T)
d2<-otu[,2:61]
d2[d2>0]<-1

health<-d2[,1:30]
health<-cbind(otu,health)
disease<-d2[,31:60]
disease<-cbind(otu,disease)

health$sum<-rowSums(health[,62:91],na.rm = T)
disease$sum<-rowSums(disease[,62:91],na.rm = T)

health<- subset(health, sum> 3 )
disease<- subset(disease, sum > 3)

a=health$Genus
b=disease$Genus
a<-as.data.frame(a)
b<-as.data.frame(b)
#write.csv(a, file= "a.csv",row.names = T ,quote = F)
#write.csv(b, file= "b.csv",row.names = T ,quote = F)
#write.csv(c, file= "c.csv",row.names = T ,quote = F)
#write.csv(d, file= "d.csv",row.names = T ,quote = F)
colnames(a)<-"Genus"
colnames(b)<-"Genus"

e=rbind(a,b)
e=unique(e)
vars=e$Genus

d1<-as.data.frame(t(otu))
d1$sample_ID<-rownames(d1)
colnames(d1)<-otu$Genus
d1<-d1[2:61,]

clean<-d1[vars]
#write.csv(clean, file= "3.4.genus_abs_match_clean.csv",row.names = T ,quote = F)
phylum<-read.csv("3.4.genus_abs_match_clean.csv",header = T)
p1<-phylum[,2:61]
p1<-as.data.frame(t(p1))
p1=as.data.frame(lapply(p1,as.numeric))
phylum_relative <- p1 / rowSums(p1)*100

phylum_relative<-as.data.frame(t(phylum_relative))
colnames(phylum_relative)<-colnames(phylum[,2:61])
rownames(phylum_relative)<-phylum$Genus
write.csv(phylum_relative, file= "3.4.genus_rel_match_clean.csv",row.names = T ,quote = F)


#lefse数据表清理
#http://www.ehbio.com/Cloud_Platform/front/#/analysis?page=b%27MzY%3D%27
#https://cloud.tencent.com/developer/article/1936404
lefse<-read.csv("2.asv_absolute_clean_match.csv",header = T)
taxonomy<-lefse[,8:67]
taxonomy=as.data.frame(lapply(taxonomy,as.numeric))
taxonomy$phylum<-lefse$phylum
taxonomy$phylum<-as.factor(taxonomy$phylum)
class<-aggregate(taxonomy[,1:60],list(taxonomy$phylum),sum)

p1<-class[,2:61]
p1<-as.data.frame(t(p1))
p1=as.data.frame(lapply(p1,as.numeric))
ASV_relative_match <- p1 / rowSums(p1)*100

ASV_relative_match<-as.data.frame(t(ASV_relative_match))
colnames(ASV_relative_match)<-colnames(lefse[,8:67])
rownames(ASV_relative_match)<-class$Group.1
write.csv(ASV_relative_match, file= "3.4.lefse_p_match.csv",row.names = T ,quote = F)


