
setwd("E:/00.COSAS/02.RESISTOME_MASTER/000.Revision_Nat_Microb/Figure_1")

library(vegan)
library(devtools)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
#setwd("D:/0000.COSAS/02.RESISTOME_MASTER/")


#data<-read.csv("ARGs_cpms_genelevel.txt", sep="\t", header = TRUE, check.names=F)
data<-read.csv("CPM_Family_Group_Elena.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors = FALSE, )
meta<-read.csv("MASTER_metadata_JNJ-2023-03-10.txt", sep="\t", header = TRUE, row.names=1)
#listado<-read.csv("ARGs_genefam.txt", sep="\t", header = TRUE, check.names=F)

#data<-data[,-ncol(data)]
data<-data[,order(colMeans(data),decreasing = TRUE)]
#listado<-listado[order(listado$Hit),]
all(rownames(data)==rownames(meta))

data<-data[-which(meta$Surface %in% c("Negative_control_industry", "Negative_control_laboratory", "Food_contact_defect", "Positive_control_10e6_CFU-mL",
	"Positive_control_10e4_CFU-mL", "Positive_control_10e2_CFU-mL", "Operator", "-", "Pump","Intermediate_ripening")),]
meta<-meta[-which(meta$Surface %in% c("Negative_control_industry", "Negative_control_laboratory", "Food_contact_defect", "Positive_control_10e6_CFU-mL",
	"Positive_control_10e4_CFU-mL", "Positive_control_10e2_CFU-mL", "Operator", "-", "Pump","Intermediate_ripening")),]
# we have removed previously "intermediate ripening" and with these 2 lines we change before_ripening to raw_material on meat and to final_product on cheese
meta$Surface[which(meta$Surface=="Before_ripening")]<-meta$Room[which(meta$Surface=="Before_ripening")]
meta$Surface[which(meta$Surface=="Raw_materials")]<-"Raw_material"

meta<-meta[which(rowSums(data)!=0),]
data<-data[which(rowSums(data)!=0),]

data3<-cbind(meta$Surface,meta$Industry_type,meta$Room,data)
#data3<-data3a[,-c(21,40)]	# to remove qac and heat
#data3$Sum<-rowSums(data3[,(ncol(meta)+1):ncol(data3)])
colnames(data3)[1:3]<-c("Surface","Industry_type","Room")
data3$Surface<-as.factor(data3$Surface)
data3$Industry_type<-as.factor(data3$Industry_type)
data3$Room<-as.factor(data3$Room)
#families2<-fam1[,order(-colMeans(fam1))]
#families2<-families2[,-c(1,5)] # to remove Heat and Quaternary Amonium Compounds

data3<-data3[which(data3$Room %in% c("Raw_material","Processing","Drying","Cold","Smoking",
  "Ripening","Packing","Final_product")),]

#data4<-data4[,-c(2:5)]
#data4$Sum<-rowSums(data3)
gathered<-as.data.frame(gather(data3, key="Family", value="CPMs",-colnames(data3)[1:3]))
levels(gathered$Family)<-factor(c(colnames(data3)[4:ncol(data3)]))

gathered$Room<-factor(gathered$Room, levels=c("Raw_material","Processing","Drying","Cold","Smoking",
                                                      "Ripening","Packing","Final_product"))
gathered$Surface<-factor(gathered$Surface, levels=c("Raw_material","Food_contact","Non_food_contact","Final_product"))
#"Operator"            
#"Delivery"              
#"Water_ice"
#"-"
library(RColorBrewer)
colorines<-c(brewer.pal(12,"Paired"),brewer.pal(6,"Dark2"))


pdf("barplot_ROOM.pdf", height = 6, width = 12)
ggplot(aes(y = CPMs, x = Room), data = gathered) +
  geom_bar(aes(fill = factor(Family, levels=levels(Family))),
           #         stat = "identity") +
           stat = "summary", fun = "mean") +
  scale_y_continuous() +
  facet_grid(~Industry_type+Surface, scales = "free", space = "free") +
  labs(y= "AMRGs (counts per million reads)")+#, fill = "Family") + #fill makes reference to the legend
  theme_bw()+
  labs(fill = "Antibiotic family") + # To change legend title
  scale_fill_manual(values = colorines) +
  #  guide_legend(title="Antibiotic family") +
  theme(#axis.title.x=element_blank(),
    axis.text.x=element_text(angle=-90),
#    axis.ticks.x=element_blank(),
    legend.position="bottom"
    #    legend.text = element_text(face = "italic")
  )
dev.off()
