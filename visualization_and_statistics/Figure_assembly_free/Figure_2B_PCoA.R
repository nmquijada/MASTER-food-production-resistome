
setwd("D:/0000.COSAS/02.RESISTOME_MASTER/Figuras_Final/Figure_1")

library(vegan)
library(devtools)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
#setwd("D:/0000.COSAS/02.RESISTOME_MASTER/")


data<-read.csv("ARGs_cpms_genelevel.txt", sep="\t", header = TRUE, check.names=F)
meta<-read.csv("MASTER_metadata_JNJ-2023-03-10.txt", sep="\t", header = TRUE, row.names=1)
listado<-read.csv("ARGs_genefam.txt", sep="\t", header = TRUE, check.names=F)

data<-data[,-ncol(data)]
data<-data[,order(colnames(data))]
listado<-listado[order(listado$Hit),]
table(colnames(data)==listado$Hit)

data<-data[-which(meta$Surface %in% c("Negative_control_industry", "Negative_control_laboratory", "Food_contact_defect", "Positive_control_10e6_CFU-mL",
	"Positive_control_10e4_CFU-mL", "Positive_control_10e2_CFU-mL", "Operator", "-", "Pump","Intermediate_ripening")),]
meta<-meta[-which(meta$Surface %in% c("Negative_control_industry", "Negative_control_laboratory", "Food_contact_defect", "Positive_control_10e6_CFU-mL",
	"Positive_control_10e4_CFU-mL", "Positive_control_10e2_CFU-mL", "Operator", "-", "Pump","Intermediate_ripening")),]
# we have removed previously "intermediate ripening" and with these 2 lines we change before_ripening to raw_material on meat and to final_product on cheese
meta$Surface[which(meta$Surface=="Before_ripening")]<-meta$Room[which(meta$Surface=="Before_ripening")]
meta$Surface[which(meta$Surface=="Raw_materials")]<-"Raw_material"
meta$Industry_type[which(meta$Industry_type=="Cheese_Dairy")]<-"Cheese_dairy"

#meta<-meta[which(rowSums(data)!=0),]
#data<-data[which(rowSums(data)!=0),]
dataa<-aggregate(t(data), list(listado$Genefam), FUN=sum)
rownames(dataa)<-dataa$Group.1
data<-as.data.frame(t(dataa[-1]))

meta<-meta[which(rowSums(data)!=0),]
data<-data[which(rowSums(data)!=0),]


data_FC<-data[which(meta$Industry_type=="Meat"),]
data_NFC<-data[which(meta$Industry_type=="Cheese_dairy"),]
data_raw<-data[which(meta$Industry_type=="Processed_fish"),]
data_final<-data[which(meta$Industry_type=="Vegetables"),]
meta_FC<-meta[which(meta$Industry_type=="Meat"),]
meta_NFC<-meta[which(meta$Industry_type=="Cheese_dairy"),]
meta_raw<-meta[which(meta$Industry_type=="Processed_fish"),]
meta_final<-meta[which(meta$Industry_type=="Vegetables"),]

droplevels(meta_FC$Surface)
droplevels(meta_NFC$Surface)
droplevels(meta_raw$Surface)
droplevels(meta_final$Surface)

speciesPcoa <- cmdscale(vegdist(data_FC, method = 'hellinger'), k = 2, eig = TRUE)
speciesPoints <- as.data.frame(speciesPcoa$points)
eigs <- eigenvals(speciesPcoa)
eigs1 <- (eigs/sum(eigs))
speciesPoints1 <- setNames(data.frame(meta_FC, speciesPoints),
                          c(colnames(meta_FC), 'PC1', 'PC2'))

speciesPcoa <- cmdscale(vegdist(data_NFC, method = 'hellinger'), k = 2, eig = TRUE)
speciesPoints <- as.data.frame(speciesPcoa$points)
eigs <- eigenvals(speciesPcoa)
eigs2 <- (eigs/sum(eigs))
speciesPoints2 <- setNames(data.frame(meta_NFC, speciesPoints),
                          c(colnames(meta_NFC), 'PC1', 'PC2'))

speciesPcoa <- cmdscale(vegdist(data_raw, method = 'hellinger'), k = 2, eig = TRUE)
speciesPoints <- as.data.frame(speciesPcoa$points)
eigs <- eigenvals(speciesPcoa)
eigs3 <- (eigs/sum(eigs))
speciesPoints3 <- setNames(data.frame(meta_raw, speciesPoints),
                          c(colnames(meta_raw), 'PC1', 'PC2'))

speciesPcoa <- cmdscale(vegdist(data_final, method = 'hellinger'), k = 2, eig = TRUE)
speciesPoints <- as.data.frame(speciesPcoa$points)
eigs <- eigenvals(speciesPcoa)
eigs4 <- (eigs/sum(eigs))
speciesPoints4 <- setNames(data.frame(meta_final, speciesPoints),
                          c(colnames(meta_final), 'PC1', 'PC2'))

adon1<-adonis(vegdist(data_FC, method = 'hellinger') ~ meta_FC$Surface)
adon2<-adonis(vegdist(data_NFC, method = 'hellinger') ~ meta_NFC$Surface)
adon3<-adonis(vegdist(data_raw, method = 'hellinger') ~ meta_raw$Surface)
adon4<-adonis(vegdist(data_final, method = 'hellinger') ~ meta_final$Surface)
adon<-rbind(adon1$aov.tab,adon2$aov.tab,adon3$aov.tab,adon4$aov.tab)
write.table(adon, "Adonis_table_industry_hellinger.txt", sep="\t", quote=FALSE)

mycolors<-c("Meat"="#FF8A65","Cheese_Dairy"="#FFCA28","Processed_fish"="#29B6F6","Vegetables"="#2E7D32")
mycolors<-c("Raw_material"="seagreen2" , "Final_product"="seagreen4", "Food_contact"="gold","Non_food_contact"="coral2")

(pcoa_FC<-ggplot(speciesPoints1, aes(x = PC1, y = PC2, fill = Surface)) +
  stat_ellipse(aes(fill = Surface), size = 1.1, show.legend = F, geom = 'polygon', alpha = 0.4) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_colour_manual(values=mycolors)+
  scale_fill_manual(values=mycolors)+ 
  geom_point(size = 1, pch = 21, colour = 'black') +
  xlab(paste0("PC 1 [", round(eigs1[1]*100, digits = 2), "%]")) +
  ylab(paste0("PC 2 [", round(eigs1[2]*100, digits = 2), "%]")) +
  theme_bw() +
    ggtitle("Meat") +
  theme(#legend.title = element_text(sample),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.border = element_rect(colour = 'black', linetype = 1, size = 1)) 
)

(pcoa_NFC<-ggplot(speciesPoints2, aes(x = PC1, y = PC2, fill = Surface)) +
  stat_ellipse(aes(fill = Surface), size = 1.1, show.legend = F, geom = 'polygon', alpha = 0.4) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_colour_manual(values=mycolors)+
  scale_fill_manual(values=mycolors)+ 
  geom_point(size = 1, pch = 21, colour = 'black') +
  xlab(paste0("PC 1 [", round(eigs2[1]*100, digits = 2), "%]")) +
  ylab(paste0("PC 2 [", round(eigs2[2]*100, digits = 2), "%]")) +
  theme_bw() +
    ggtitle("Cheese/dairy") +
  theme(#legend.title = element_text(sample),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.border = element_rect(colour = 'black', linetype = 1, size = 1)) 
)

(pcoa_raw<-ggplot(speciesPoints3, aes(x = PC1, y = PC2, fill = Surface)) +
  stat_ellipse(aes(fill = Surface), size = 1.1, show.legend = F, geom = 'polygon', alpha = 0.4) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_colour_manual(values=mycolors)+
  scale_fill_manual(values=mycolors)+ 
  geom_point(size = 1, pch = 21, colour = 'black') +
  xlab(paste0("PC 1 [", round(eigs3[1]*100, digits = 2), "%]")) +
  ylab(paste0("PC 2 [", round(eigs3[2]*100, digits = 2), "%]")) +
  theme_bw() +
    ggtitle("Processed fish") +
  theme(#legend.title = element_text(sample),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.border = element_rect(colour = 'black', linetype = 1, size = 1)) 
)

(pcoa_final<-ggplot(speciesPoints4, aes(x = PC1, y = PC2, fill = Surface)) +
  stat_ellipse(aes(fill = Surface), size = 1.1, show.legend = F, geom = 'polygon', alpha = 0.4) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_colour_manual(values=mycolors)+
  scale_fill_manual(values=mycolors)+ 
  geom_point(size = 1, pch = 21, colour = 'black') +
  xlab(paste0("PC 1 [", round(eigs4[1]*100, digits = 2), "%]")) +
  ylab(paste0("PC 2 [", round(eigs4[2]*100, digits = 2), "%]")) +
  theme_bw() +
    ggtitle("Vegetables") +
  theme(#legend.title = element_text(sample),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.border = element_rect(colour = 'black', linetype = 1, size = 1)) 
)

pdf("PCoA_surface_hellinger.pdf", height=7, width=8, onefile=FALSE)
ggarrange(pcoa_raw, pcoa_final, pcoa_FC, pcoa_NFC, 
          #labels=c("Meat", "Cheese_dairy", "Processed_fish","Vegetables"),
          common.legend = TRUE,
          legend="right",
          ncol=2,
          nrow=2
)
dev.off()


