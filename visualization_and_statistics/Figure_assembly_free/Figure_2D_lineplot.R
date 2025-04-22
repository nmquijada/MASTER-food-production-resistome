
setwd("D:/00.COSAS/02.RESISTOME_MASTER/Figuras_Final/Figure_1")

library(ggpubr)
library(vegan)
library(tidyr)
library(ggplot2)

data<-read.csv("ARGs_cpms_genelevel.txt", sep="\t", header = TRUE, check.names=F)
meta<-read.csv("MASTER_metadata_JNJ-2023-03-10.txt", sep="\t", header = TRUE, row.names=1, stringsAsFactors=FALSE)
listado<-read.csv("ARGs_genefam.txt", sep="\t", header = TRUE, check.names=F)

data<-data[,-ncol(data)]
data<-data[,order(colnames(data))]
listado<-listado[order(listado$Hit),]
all(colnames(data)==listado$Hit)

data<-data[-which(meta$Surface %in% c("Negative_control_industry", "Negative_control_laboratory", "Food_contact_defect", "Positive_control_10e6_CFU-mL",
	"Positive_control_10e4_CFU-mL", "Positive_control_10e2_CFU-mL", "Operator", "-", "Pump","Intermediate_ripening")),]
meta<-meta[-which(meta$Surface %in% c("Negative_control_industry", "Negative_control_laboratory", "Food_contact_defect", "Positive_control_10e6_CFU-mL",
	"Positive_control_10e4_CFU-mL", "Positive_control_10e2_CFU-mL", "Operator", "-", "Pump","Intermediate_ripening")),]
# we have removed previously "intermediate ripening" and with these 2 lines we change before_ripening to raw_material on meat and to final_product on cheese
meta$Surface[which(meta$Surface=="Before_ripening")]<-meta$Room[which(meta$Surface=="Before_ripening")]
meta$Surface[which(meta$Surface=="Raw_materials")]<-"Raw_material"


#meta<-meta[which(rowSums(data)!=0),]
#data<-data[which(rowSums(data)!=0),]
dataa<-aggregate(t(data), list(listado$Genefam), FUN=sum)
rownames(dataa)<-dataa$Group.1
data<-as.data.frame(t(dataa[-1]))

data<-data[,order(colMeans(data),decreasing=TRUE)]
#data<-data[,-1]	# to remove Sums column
data2<-cbind(meta,data[,which(colMeans(data)>0.5)])

datata<-data[,which(colMeans(data)>0.5)]
datata<-datata[,order(colnames(datata))]
good_order<- match(colnames(datata),listado$Genefam)
selected<-listado[good_order,]
selected<-droplevels(selected)
selected<-selected[order(selected$Genefam),]

data2$Surface<-factor(data2$Surface,
	levels=c("Raw_material", "Final_product","Food_contact","Non_food_contact"))
data2$Industry_type<-factor(data2$Industry_type,
	levels=c(#"Aged_beef","Cured_meats","Fermented_sausages","Raw_meat",
	  "Meat","Cheese_Dairy","Processed_fish","Vegetables"))

#mycolors<-c("#795548","coral3","#EF5350","#FF8A65","#FFCA28","#29B6F6","#2E7D32")
mycolors<-c("Meat"="#FF8A65","Cheese_Dairy"="#FFCA28","Processed_fish"="#29B6F6","Vegetables"="#2E7D32")
#mycolors2<-c("Raw_materials"='grey50', "Before_ripening"="seagreen2", "Intermediate_ripening"="seagreen3", 
#             "Final_product"="seagreen4", "Food_contact"="gold","Non_food_contact"="coral2")
mycolors2<-c("Raw_material"="seagreen2",  "Final_product"="seagreen4", "Food_contact"="gold","Non_food_contact"="coral2")

library(RColorBrewer)
#mycolors3<-brewer.pal(12, "Paired")
display.brewer.pal(12,"Paired")
mycolors3<-c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F",
	"#FF7F00","#CAB2D6","#6A3D9A","#B15928")

###########   STATISTICAL DIFFERENCES    ##############################

mean2<-aggregate(data2[,20:ncol(data2)], by=list(data2$Industry_type), FUN=mean)

names<-mean2$Group.1
mean <- as.data.frame(t(mean2[,-1])) 
rownames(mean)<-colnames(data2)[20:ncol(data2)]
colnames(mean)<-names
#mean<-mean[,order(colnames(mean))]

krusk<- lapply(20:ncol(data2), function(x) kruskal.test(data2[[x]]~ as.factor(data2$Industry_type))) 
kruskout <- sapply(krusk, function(x) {p <- x$p.value})
wilk <- lapply(20:ncol(data2), function(x) pairwise.wilcox.test(data2[[x]], as.factor(data2$Industry_type), p.adjust.method = "BH")) #loop to create dataframe with the wilcoxon results
wilkout <- sapply(wilk, function(x) {p <- x$p.value})
wilkout<-wilkout[-c(4,7,8),]


gatheredmean<-as.data.frame(gather(mean2, key="Family", value="CPMs",-c("Group.1")))
pdf("lineplot_Figure2D_mean.pdf", height=4, width=3.8) 
ggplot(gatheredmean,aes(x=Family,y=CPMs,col=Group.1,group=Group.1))+geom_point()+coord_flip()+#geom_line()+
  scale_colour_manual(values = mycolors) +
  scale_y_continuous(trans='sqrt') +
  theme_bw()+
  theme(axis.text.y = element_text(face="italic"), axis.text.x = element_text(size=10, angle=270),
	        strip.background = element_blank(), strip.text = element_text(face = 'bold')) +
  xlab('Antimicrobial Resistance Gene') + ylab("Count Per Million reads")
dev.off()

write.table(wilk[[1]]$p.value, "wilcon_Host.txt", sep="\t",quote=FALSE)

kruskwilk <- cbind(rownames(mean),mean,kruskout,t(wilkout))
head<-c("ARG","Meat","Cheese/dairy","Processed fish","Vegetables","Kruskal","Wilcox_Meat-Cheese/dairy","Wilcox_Meat-Fish","Wilcox_Meat-Veget","Wilcox_Cheese-Fish","Wilcox_Cheese-Veget","Wilcox_Fish-Veget")
write.table(rbind(head,kruskwilk), "ARG_stat_Industry_type.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

######## ALL SAMPLES
#(mod.a3<-ggplot(data3, aes(x=Industry_type, y=Sum, colour=Industry_type)) + 
#    geom_boxplot(outlier.shape = NA, size = 0.8)+
#.    geom_jitter(aes(fill=Industry_type),colour = 'black', width = 0.15, height = 0, pch = 21, size = 1.0) +
data_boxplot<-data2[,c(8,10,20:ncol(data2))]
gathered<-as.data.frame(gather(data_boxplot, key="Family", value="CPMs",-c("Industry_type","Surface")))    
levels(gathered$Family)<-unique(gathered$Family)

data3<-aggregate(data2[,(ncol(meta)+1):ncol(data2)], list(data2$Industry_type,data2$Surface), FUN=mean)
#rownames(data3)<-data3$Group.1
colnames(data3)[1:2]<-c("Industry_type","Surface")
gathered2<-as.data.frame(gather(data3, key="Family", value="CPMs",-c("Industry_type","Surface")))
#gathered<-gathered[-which(gathered$Surface %in% c("Before_ripening","Intermediate_ripening")),]

datatest<-cbind(data2$Industry_type, data2$Surface, data2[,(ncol(meta)+1):ncol(data2)])
colnames(datatest)[1:2]<-c("Industry_type","Surface")
gathertest<-as.data.frame(gather(datatest,key="Family", value="CPMs",-c("Industry_type","Surface")))
stat.test1<-compare_means(CPMs ~ Surface, data = gathertest, method = "wilcox",group.by=c("Industry_type","Family"))  
stat.test.signif1<-stat.test1[stat.test1$p.signif!="ns",]%>%
  mutate(y.position = c(3.1))

pdf("lineplot_Figure2D.pdf", height=4, width=7) #width=8)
ggplot(gathered,aes(x=Family,y=CPMs,col=Surface,group=Surface))+geom_point()+coord_flip()+#geom_line()+
facet_grid(~ Industry_type, scales = "free") +
  scale_colour_manual(values = mycolors2) +
#  scale_y_continuous() +
  scale_y_continuous(trans='sqrt') +
#  stat_pvalue_manual(stat.test.signif1, label= "p.signif", linetype="blank", hjust=2,
 #                      tip.length=0, bracket.size=0.6, step.group.by=c("Industry_type")) +
  theme_bw()+
  #theme(axis.text.y = element_text(colour = selected$Family, face="italic"),
  theme(axis.text.y = element_text(face="italic"), axis.text.x = element_text(size=10, angle=270),
	        strip.background = element_blank(), strip.text = element_text(face = 'bold')) +
  xlab('Antimicrobial Resistance Gene') + ylab("Count Per Million reads")
dev.off()


######## BY SURFACES

data_cheese<-data[which(meta$Industry_type=="Cheese_Dairy"),]
data_fish<-data[which(meta$Industry_type=="Processed_fish"),]
data_meat<-data[which(meta$Industry_type=="Meat"),]
data_vegetables<-data[which(meta$Industry_type=="Vegetables"),]

meta_cheese<-meta[which(meta$Industry_type=="Cheese_Dairy"),]
meta_fish<-meta[which(meta$Industry_type=="Processed_fish"),]
meta_meat<-meta[which(meta$Industry_type=="Meat"),]
meta_vegetables<-meta[which(meta$Industry_type=="Vegetables"),]


#### Meat
data_meat<-data_meat[,order(colMeans(data_meat),decreasing=TRUE)]
data_meat2<-cbind(meta_meat,data_meat[,which(colMeans(data_meat)>0.1)])

datata<-data_meat[,which(colMeans(data_meat)>0.1)]
datata<-datata[,order(colnames(datata))]
good_order<- match(colnames(datata),listado$Genefam)
selected<-listado[good_order,]
selected<-droplevels(selected)
selected_meat<-selected[order(selected$Genefam),]

data3<-aggregate(data_meat2[,(ncol(meta_meat)+1):ncol(data_meat2)], list(data_meat2$Surface), FUN=mean)
#rownames(data3)<-data3$Group.1
colnames(data3)[1]<-"Surface"
data3$Surface<-factor(data3$Surface,
	levels=c("Raw_material", "Final_product","Food_contact","Non_food_contact"))
gathered_meat<-as.data.frame(gather(data3, key="Family", value="CPMs",-c("Surface")))

(line_meat<-ggplot(gathered_meat,aes(x=Family,y=CPMs,col=Surface,group=Surface))+geom_line()+geom_point()+coord_flip()+
#facet_grid(~ Industry_type, scales = "free", space = "free") +
#  scale_fill_manual(values = mycolors3) +
  scale_colour_manual(values = mycolors2) +
#  scale_x_continuous() +
  scale_y_continuous(trans='sqrt') +
  theme_bw()+
  #theme(axis.text.y = element_text(face="italic",colour = selected_meat$Family),
  theme(axis.text.y = element_text(),
	# axis.text.x = element_text(angle=90),
       # axis.ticks.x = element_blank(),
        strip.background = element_blank(),
#        strip.text = element_text(face = 'bold')) +
        strip.text = element_text()) +
  xlab('Antimicrobial Resistance Gene') +
  ylab("Count Per Million reads"))


#### Cheese
data_cheese<-data_cheese[,order(colMeans(data_cheese),decreasing=TRUE)]
data_cheese2<-cbind(meta_cheese,data_cheese[,which(colMeans(data_cheese)>0.1)])

datata<-data_cheese[,which(colMeans(data_cheese)>0.1)]
datata<-datata[,order(colnames(datata))]
good_order<- match(colnames(datata),listado$Genefam)
selected<-listado[good_order,]
selected<-droplevels(selected)
selected_cheese<-selected[order(selected$Genefam),]

data3<-aggregate(data_cheese2[,(ncol(meta_cheese)+1):ncol(data_cheese2)], list(data_cheese2$Surface), FUN=mean)
#rownames(data3)<-data3$Group.1
colnames(data3)[1]<-"Surface"
data3$Surface<-factor(data3$Surface,
                      levels=c("Raw_material","Final_product","Food_contact","Non_food_contact"))
gathered_cheese<-as.data.frame(gather(data3, key="Family", value="CPMs",-c("Surface")))

(line_cheese<-ggplot(gathered_cheese,aes(x=Family,y=CPMs,col=Surface,group=Surface))+geom_line()+geom_point()+coord_flip()+
    #facet_grid(~ Industry_type, scales = "free", space = "free") +
    #  scale_fill_manual(values = mycolors3) +
    scale_colour_manual(values = mycolors2) +
  scale_y_continuous(trans='sqrt') +
    #  scale_x_continuous() +
    theme_bw()+
    theme(axis.text.y = element_text(),#colour = selected_cheese$Family, face="italic"),
          # axis.text.x = element_text(angle=90),
          # axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold')) +
    xlab('Antimicrobial Resistance Gene') +
    ylab("Count Per Million reads"))

#### Fish
data_fish<-data_fish[,order(colMeans(data_fish),decreasing=TRUE)]
data_fish2<-cbind(meta_fish,data_fish[,which(colMeans(data_fish)>0.1)])

datata<-data_fish[,which(colMeans(data_fish)>0.1)]
datata<-datata[,order(colnames(datata))]
good_order<- match(colnames(datata),listado$Genefam)
selected<-listado[good_order,]
selected<-droplevels(selected)
selected_fish<-selected[order(selected$Genefam),]

data3<-aggregate(data_fish2[,(ncol(meta_fish)+1):ncol(data_fish2)], list(data_fish2$Surface), FUN=mean)
#rownames(data3)<-data3$Group.1
colnames(data3)[1]<-"Surface"
data3$Surface<-factor(data3$Surface,
                      levels=c("Raw_material","Final_product","Food_contact","Non_food_contact"))
gathered_fish<-as.data.frame(gather(data3, key="Family", value="CPMs",-c("Surface")))

(line_fish<-ggplot(gathered_fish,aes(x=Family,y=CPMs,col=Surface,group=Surface))+geom_line()+geom_point()+coord_flip()+
    #facet_grid(~ Industry_type, scales = "free", space = "free") +
    #  scale_fill_manual(values = mycolors3) +
    scale_colour_manual(values = mycolors2) +
  scale_y_continuous(trans='sqrt') +
    #  scale_x_continuous() +
    theme_bw()+
    theme(axis.text.y = element_text(), #colour = selected_fish$Family, face="italic"),
          # axis.text.x = element_text(angle=90),
          # axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold')) +
    xlab('Antimicrobial Resistance Gene') +
    ylab("Count Per Million reads"))


#### Vegetables
data_vegetables<-data_vegetables[,order(colMeans(data_vegetables),decreasing=TRUE)]
data_vegetables2<-cbind(meta_vegetables,data_vegetables[,which(colMeans(data_vegetables)>0.1)])

datata<-data_vegetables[,which(colMeans(data_vegetables)>0.1)]
datata<-datata[,order(colnames(datata))]
good_order<- match(colnames(datata),listado$Genefam)
selected<-listado[good_order,]
selected<-droplevels(selected)
selected_vegetables<-selected[order(selected$Genefam),]

data3<-aggregate(data_vegetables2[,(ncol(meta_vegetables)+1):ncol(data_vegetables2)], list(data_vegetables2$Surface), FUN=mean)
#rownames(data3)<-data3$Group.1
colnames(data3)[1]<-"Surface"
data3$Surface<-factor(data3$Surface,
                      levels=c("Raw_material","Final_product","Food_contact","Non_food_contact"))
gathered_vegetables<-as.data.frame(gather(data3, key="Family", value="CPMs",-c("Surface")))

(line_vegetables<-ggplot(gathered_vegetables,aes(x=Family,y=CPMs,col=Surface,group=Surface))+geom_line()+geom_point()+coord_flip()+
    #facet_grid(~ Industry_type, scales = "free", space = "free") +
    #  scale_fill_manual(values = mycolors3) +
    scale_colour_manual(values = mycolors2) +
  scale_y_continuous(trans='sqrt') +
    #  scale_x_continuous() +
    theme_bw()+
    theme(axis.text.y = element_text(), #colour = selected_vegetables$Family, face="italic"),
          # axis.text.x = element_text(angle=90),
          # axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold')) +
    xlab('Antimicrobial Resistance Gene') +
    ylab("Count Per Million reads"))


pdf("Lineplot_industry_individual_sqrt.pdf", height=10, width=15, onefile=FALSE)
ggarrange(line_meat, line_cheese, line_fish, line_vegetables, 
          labels=c("Meat","Dairy_Cheese", "Fish", "Vegetables"),
          common.legend = TRUE,
          legend="right",
          ncol=4#,nrow=2
)
dev.off()


