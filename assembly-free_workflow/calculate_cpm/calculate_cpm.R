library(ggplot2)
library(ggpubr)
library(tidyr)

data<-read.csv("matrix_ResFinder.txt", sep="\t", header=TRUE, stringsAsFactors = TRUE, check.names = FALSE)
reads<-read.table("viromeqc_all.txt", sep=" ", header=TRUE)

datagenes<-data[,c(1:4)]
data<-data[,-c(1:4)]
#data<-data[,order(colnames(data))]
#reads<-reads[order(reads$Sample),]

reads$Sample[1747:nrow(reads)]=gsub("_1","",reads$Sample[1747:nrow(reads)])

table(colnames(data)==gsub("_R1","",reads$Sample))
colnames(data)[559]
reads$Sample[559]

#colnames(data)=gsub(".1","",colnames(data))

fam<-datagenes[,1]
#data<-as.data.frame(data[,-c(1,2)])
families<-aggregate(data, by=list(fam), FUN=sum)
rownames(families)<-families$Group.1
families<-families[,-1]

data2<-families
for (i in 1:nrow(families)) { 
  for (j in 1:ncol(families)) { data2[i,j]<-families[i,j]*1000000*reads$Bacterial_Markers.alignment.rate[j]/reads$Reads[j]}
}

genefam<-datagenes[,1]
#data<-as.data.frame(data[,-c(1,2)])
genefamilies<-aggregate(data, by=list(genefam), FUN=sum)
rownames(genefamilies)<-genefamilies$Group.1
genefamilies<-genefamilies[,-1]

data3<-genefamilies
for (i in 1:nrow(genefamilies)) { 
  for (j in 1:ncol(genefamilies)) { data3[i,j]<-genefamilies[i,j]*1000000*reads$Bacterial_Markers.alignment.rate[j]/reads$Reads[j]}
}

write.table(data2, "data_ARGs_cpms_family.csv")
write.table(data3, "data_ARGs_cpms_genefamily.csv")
