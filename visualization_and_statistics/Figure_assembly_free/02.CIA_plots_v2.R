
setwd("D:/00.COSAS/02.RESISTOME_MASTER/Figuras_Final/plots_CIA")
#Check extra column on the first 2 lines (genes and family)
library(ggplot2)
library(ggpubr)
library(tidyr)

data<-read.csv("ARGs_cpms_genelevel.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE, check.names=F)
meta<-read.csv("MASTER_metadata_JNJ-2023-03-10.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE)

all(rownames(data)==rownames(meta))

data<-data[-which(meta$Surface %in% c("Negative_control_industry", "Negative_control_laboratory", "Food_contact_defect", "Positive_control_10e6_CFU-mL",
  "Intermediate_ripening","Positive_control_10e4_CFU-mL", "Positive_control_10e2_CFU-mL", "Operator", "-", "Pump")),]
meta<-meta[-which(meta$Surface %in% c("Negative_control_industry", "Negative_control_laboratory", "Food_contact_defect", "Positive_control_10e6_CFU-mL",
	"Intermediate_ripening","Positive_control_10e4_CFU-mL", "Positive_control_10e2_CFU-mL", "Operator", "-", "Pump")),]

meta$Surface[which(meta$Surface=="Before_ripening")]<-meta$Room[which(meta$Surface=="Before_ripening")]
meta$Surface[which(meta$Surface=="Raw_materials")]<-"Raw_material"

listado<-read.csv("ARGs_fam_gene_hit.txt", sep="\t", header = TRUE, check.names=F)

data<-data[,-ncol(data)]
dataa<-aggregate(t(data), list(listado$Gene), FUN=sum)
rownames(dataa)<-dataa$Group.1
data<-t(dataa[-1])

selected<-read.csv("list_cia.txt", sep="\t", header = TRUE, check.names=F)

data2<-data[,which(colnames(data) %in% selected$gene)]

data3<-cbind(meta,data2)
#data3<-data3a[,-c(21,40)]	# to remove qac and heat
data3$Sum<-rowSums(data3[,(ncol(meta)+1):ncol(data3)])

data3$Surface<-as.factor(data3$Surface)
data3$Industry_type<-as.factor(data3$Industry_type)

data3$Surface<-factor(data3$Surface, levels=c("Raw_material","Final_product","Food_contact","Non_food_contact"))
data3$Industry_type<-factor(data3$Industry_type, levels=c("Meat","Cheese_Dairy","Processed_fish","Vegetables"))

#mycolors<-c("#795548","coral3","#EF5350","#FF8A65","#FFCA28","#29B6F6","#2E7D32")
mycolors<-c("#FF8A65","#FFCA28","#29B6F6","#2E7D32")
mycolors2<-c("seagreen2" , "seagreen4", "gold","coral2")
  



stat.test1<-compare_means(Sum ~ Industry_type, data = data3, method = "wilcox")#,group.by="Industry_type")  
stat.test.signif1<-stat.test1[stat.test1$p.signif!="ns",]%>%
  mutate(y.position = c(4.5))

aggregate(data3$Sum, by=list(data2$Industry_type), FUN=mean)
aggregate(data3$Sum, by=list(data2$Industry_type), FUN=sd)
#data3$Sum[which(data3$Sum==0)]<-0.00000001

(mod.a1<-ggplot(data3, aes(x=Industry_type, y=Sum, colour=Industry_type)) + 
    geom_boxplot(outlier.shape = NA, size = 0.8)+
    geom_jitter(aes(fill=Surface),colour = 'black', width = 0.15, height = 0, pch = 21, size = 1.0) +
    scale_colour_manual(values=mycolors)+
    scale_fill_manual(values=mycolors2)+ 
#    stat_summary(fun=mean,shape=1,fill='black',geom='point')+
    #facet_grid(~Industry_type, scales="free", space="free")+
    scale_y_continuous(trans='log10') +
    stat_pvalue_manual(stat.test.signif1, 
                       tip.length=0.01, bracket.size=0.6,# step.group.by=c("Industry_type"),
                       step.increase=c(0.05), label="p.signif") +
    theme_bw()+
    theme(legend.position = "right",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold')) +
    xlab(' ') +
    ylab("ARGs (counts per million reads)"))

stat.test2<-compare_means(Sum ~ Surface, data = data3, method = "wilcox" ,group.by="Industry_type")  
stat.test.signif2<-stat.test2[stat.test2$p.signif!="ns",]%>%
  mutate(y.position = c(3.5))

(mod.a2<-ggplot(data3, aes(x=Surface, y=Sum, colour=Surface)) + 
    geom_boxplot(outlier.shape = NA, size = 0.8)+
    geom_jitter(aes(fill=Surface),colour = 'black', width = 0.15, height = 0, pch = 21, size = 1.0) +
    scale_colour_manual(values=mycolors2)+
    scale_fill_manual(values=mycolors2)+ 
    facet_grid(~Industry_type, scales="free", space="free")+
    scale_y_continuous(trans='log10') +
    stat_pvalue_manual(stat.test.signif2, 
                       tip.length=0.005, bracket.size=0.5, step.group.by=c("Industry_type"),
                       step.increase=c(0.05), label="p.signif") +
    theme_bw()+
    theme(legend.position = "right",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold')) +
    xlab(' ') +
    ylab("ARGs (CPM per bacterial markers)"))


stat.test3<-compare_means(Sum ~ Industry_type, data = data3, method = "wilcox" ,group.by="Surface")  
stat.test.signif3<-stat.test3[stat.test3$p.signif!="ns",]%>%
  mutate(y.position = c(3.5))

(mod.a3<-ggplot(data3, aes(x=Industry_type, y=Sum, colour=Industry_type)) + 
    geom_boxplot(outlier.shape = NA, size = 0.8)+
    geom_jitter(aes(fill=Industry_type),colour = 'black', width = 0.15, height = 0, pch = 21, size = 1.0) +
    scale_colour_manual(values=mycolors)+
    scale_fill_manual(values=mycolors)+ 
    facet_grid(~Surface, scales="free", space="free")+
    scale_y_continuous(trans='log10') +
    stat_pvalue_manual(stat.test.signif3, 
                       tip.length=0.005, bracket.size=0.5, step.group.by=c("Surface"),
                       step.increase=c(0.05), label="p.signif") +
    theme_bw()+
    theme(legend.position = "right",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold')) +
    xlab(' ') +
    ylab("ARGs (CPM per bacterial markers)"))

pdf("resistome_reads_CIA_v2.pdf",width = 10, height = 9)
ggarrange(mod.a2, mod.a3, nrow=2, ncol=1)
dev.off()
pdf("resistome_reads_CIA_v1.pdf",width = 5, height = 4)
mod.a1
dev.off()


############################################################## VAAAAAAMOS

cia_data<-data3[,20:(ncol(data3)-1)]
cia_list<-listado[which(listado$Gene %in% colnames(cia_data)),]
cia_list<-cia_list[!duplicated(cia_list[,'Gene']),]


#ORDENAR CIALIST Y CIADATA PARA HACER EL BARPLOT

fam1<-cia_data
families2<-fam1[,order(-colMeans(fam1))]
#families2<-families2[,-c(1,5)] # to remove Heat and Quaternary Amonium Compounds
data4<-cbind(data3[,1:ncol(meta)],families2[,1:20],rowSums(families2[,21:ncol(families2)]))
colnames(data4)[ncol(data4)]<-"Other antibiotic families"

#data4<-data4[,-c(2:5)]
#data4$Sum<-rowSums(data3)
gathered<-as.data.frame(gather(data4, key="Family", value="CPMs",-colnames(data4)[1:ncol(meta)]))
levels(gathered$Family)<-factor(c(colnames(data4)[(ncol(meta)+1):ncol(data4)]))
gathered2<-aggregate(gathered$CPMs, by=list(gathered$Industry_type, gathered$Surface,gathered$Family), FUN=mean )
colnames(gathered2)<-c("Industry_type","Surface","Family","CPMs")
levels(gathered2$Family)<-levels(gathered$Family)
library(RColorBrewer)
colorines<-c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),"grey")

pdf("barplot_CIA_argfamilies.pdf", height = 5, width = 8)
ggplot(aes(y = CPMs, x = Surface), data = gathered2) +
  geom_bar(aes(fill = factor(Family, levels=levels(Family))),
  stat = "summary", fun = "mean") +
 # stat = "summary", fun = "mean", position="fill") +
  scale_y_continuous() +
  facet_grid(~Industry_type, scales = "free", space = "free") +
  labs(y= "ARGs (counts per million reads)")+#, fill = "Family") + #fill makes reference to the legend
  theme_bw()+
  labs(fill = "Antibiotic family") + # To change legend title
  scale_fill_manual(values = colorines) +
#  guide_legend(title="Antibiotic family") +
  theme(#axis.title.x=element_blank(),
    axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.ticks.x=element_blank(),
    legend.position="bottom"
#    legend.text = element_text(face = "italic")
    )
dev.off()

pdf("barplot_CIA_argfamilies_fill.pdf", height = 5, width = 8)
ggplot(aes(y = CPMs, x = Surface), data = gathered2) +
  geom_bar(aes(fill = factor(Family, levels=levels(Family))),
  stat = "summary", fun = "mean", position="fill") +
  scale_y_continuous() +
  facet_grid(~Industry_type, scales = "free", space = "free") +
  labs(y= "ARGs (counts per million reads)")+#, fill = "Family") + #fill makes reference to the legend
  theme_bw()+
  labs(fill = "Antibiotic family") + # To change legend title
  scale_fill_manual(values = colorines) +
  theme(#axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
    axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x=element_blank(),
        legend.position="bottom"
  )
dev.off()

