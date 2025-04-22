library(ggplot2)
library(ggpubr)
library(tidyr)
library(vegan)
library(dplyr)

data<-read.csv("ARGs_cpms_genelevel.txt", sep="\t", header = TRUE, check.names=F)
datafam<-read.csv("CPM_Family_Group_Elena.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors = FALSE, )
meta<-read.csv("MASTER_metadata_JNJ-2023-03-10.txt", sep="\t", header = TRUE, row.names=1)
listado<-read.csv("ARGs_genefam.txt", sep="\t", header = TRUE, check.names=F)

all(colnames(data)[-ncol(data)]==listado$Hit)
all(rownames(meta)==rownames(datafam))

data<-data[-which(meta$Surface %in% c("Negative_control_industry", "Negative_control_laboratory", "Food_contact_defect", "Positive_control_10e6_CFU-mL",
                                      "Positive_control_10e4_CFU-mL", "Positive_control_10e2_CFU-mL", "Operator", "-", "Pump", "Intermediate_ripening")),]
datafam<-datafam[-which(meta$Surface %in% c("Negative_control_industry", "Negative_control_laboratory", "Food_contact_defect", "Positive_control_10e6_CFU-mL",
                                      "Positive_control_10e4_CFU-mL", "Positive_control_10e2_CFU-mL", "Operator", "-", "Pump", "Intermediate_ripening")),]

meta<-meta[-which(meta$Surface %in% c("Negative_control_industry", "Negative_control_laboratory", "Food_contact_defect", "Positive_control_10e6_CFU-mL",
                                      "Positive_control_10e4_CFU-mL", "Positive_control_10e2_CFU-mL", "Operator", "-", "Pump","Intermediate_ripening")),]

# we have removed previously "intermediate ripening" and with these 2 lines we change before_ripening to raw_material on meat and to final_product on cheese
meta$Surface[which(meta$Surface=="Before_ripening")]<-meta$Room[which(meta$Surface=="Before_ripening")]
meta$Surface[which(meta$Surface=="Raw_materials")]<-"Raw_material"

#meta<-meta[which(rowSums(data)!=0),]
#data<-data[which(rowSums(data)!=0),]
data<-data[-ncol(data)]
dataa<-aggregate(t(data), list(listado$Genefam), FUN=sum)
rownames(dataa)<-dataa$Group.1

data2<-cbind(meta,t(dataa[-1]))

data2$Industry_type<-factor(data2$Industry_type,levels=c("Meat","Cheese_Dairy","Processed_fish","Vegetables"))
                      
data2$Surface<-factor(data2$Surface,levels=c("Raw_material","Final_product","Food_contact","Non_food_contact"))
                            

data2$Sum<-rowSums(data2[,(ncol(meta)+1):ncol(data2)])

mycolors<-c("#FF8A65","#FFCA28","#29B6F6","#2E7D32")
#mycolors2<-c('grey50', "seagreen2", "seagreen3", "seagreen4", "gold","coral2")
mycolors2<-c("seagreen2" , "seagreen4", "gold","coral2")
  
#### TOTAL ARGs (Counts per Mreads)

stat.test1<-compare_means(Sum ~ Surface, data = data2, method = "wilcox")#,group.by="Surface")  
stat.test.signif1<-stat.test1[stat.test1$p.signif!="ns",]%>%
  mutate(y.position = c(3.1))
aggregate(data2$Sum, by=list(data2$Surface), FUN=mean)
aggregate(data2$Sum, by=list(data2$Surface), FUN=sd)

(plot_count<-ggplot(data2, aes(x=Surface, y=Sum, colour=Surface)) + 
    geom_boxplot(outlier.shape = NA, size = 0.8)+
    geom_jitter(aes(fill=Industry_type),colour = 'black', width = 0.25, height = 0, pch = 21, size = 1.5) +
    scale_colour_manual(values=mycolors2)+
    scale_fill_manual(values=mycolors)+ 
#    stat_summary(fun=mean,shape=1,fill='black',geom='point')+
    #facet_grid(~Surface, scales="free", space="free")+
    #stat_summary(fun.y=mean, geom="point", shape=23, size=3, color="black", fill="black") +
    scale_y_continuous(trans='log10') +
    stat_pvalue_manual(stat.test.signif1, 
                       tip.length=0.01, bracket.size=0.6,# step.group.by=c("Surface"),
                       #step.increase=c(0.025), label="p.signif") +
                       step.increase=c(0.05), label= "p  = {signif(p,2)}") +
    theme_bw()+
    theme(legend.position = "right",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold')) +
    xlab(' ') +
    ylab("ARGs (counts per million reads)"))

#### ARG alpha-diversity

data2$simpson<-diversity(data, index="simpson",MARGIN=1)
data2$specnumber<-specnumber(data, MARGIN=1)

aggregate(data2$specnumber, by=list(data2$Surface), FUN=mean)
aggregate(data2$specnumber, by=list(data2$Surface), FUN=sd)
aggregate(data2$simpson, by=list(data2$Surface), FUN=mean)
aggregate(data2$simpson, by=list(data2$Surface), FUN=sd)



stat.test <- compare_means(specnumber ~ Surface, data = data2,method = "wilcox.test") %>%  
  mutate(y.position = c(450))
stat.test.signif2<-stat.test[stat.test$p.signif!="ns",]
(plot_rich<-ggplot(data2, aes(x=Surface, y=specnumber, colour=Surface))+
    geom_boxplot(outlier.shape = NA, size = 1.0) +
    geom_jitter(aes(fill = Industry_type), colour = 'black', width = 0.25, pch = 21, size = 1.5) +
    scale_fill_manual(values = mycolors) +
    scale_colour_manual(values = mycolors2) +
    scale_y_continuous() +
    stat_pvalue_manual(stat.test.signif2, tip.length=0.01, bracket.size=0.6,
                       step.increase=c(0.05), label= "p  = {signif(p,2)}") +
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold')) +
    xlab('Industry type') +
    ylab("Richness"))

stat.test <- compare_means(simpson ~ Surface, data = data2,method = "wilcox.test") %>%  
  #stat.test <- compare_means(simpson ~ Time, data = data3,method = "kruskal.test") %>%  
  mutate(y.position = c(1.05))
stat.test.signif3<-stat.test[stat.test$p.signif!="ns",]
(plot_simp<-ggplot(data2, aes(x=Surface, y=simpson, colour=Surface))+
    geom_boxplot(outlier.shape = NA, size = 1.0) +
    geom_jitter(aes(fill = Industry_type), colour = 'black', width = 0.25, pch = 21, size = 1.5) +
    scale_fill_manual(values = mycolors) +
    scale_colour_manual(values = mycolors2) +
    scale_y_continuous() +
    stat_pvalue_manual(stat.test.signif3, tip.length=0.01, bracket.size=0.6,
                       step.increase=c(0.05), label= "p  = {signif(p,2)}") +
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold')) +
    ylim(0,1.2)+
    xlab('Industry type') +
    ylab("Simpson's index"))

pdf("Alphadiv.pdf", height=7, width=12, onefile=FALSE)
ggarrange(plot_count,plot_rich, plot_simp,
          labels=c("A","",""),
          common.legend = TRUE,
          legend="bottom",
          ncol=3,
          nrow=1
)
dev.off()

############################################################## VAAAAAAMOS
### ROOM_SURFACE

stat.test <- compare_means(specnumber ~ Surface + Room, data = data2,method = "wilcox.test", group.by = "Industry_type") %>%  
  mutate(y.position = c(450))
stat.test.signif2<-stat.test[stat.test$p.signif!="ns",]
(plot_rich<-ggplot(data2, aes(x=Surface, y=specnumber, colour=Surface))+
    geom_boxplot(outlier.shape = NA, size = 1.0) +
    geom_jitter(aes(fill = Industry_type), colour = 'black', width = 0.25, pch = 21, size = 1.5) +
    scale_fill_manual(values = mycolors) +
    scale_colour_manual(values = mycolors2) +
    scale_y_continuous() +
    stat_pvalue_manual(stat.test.signif2, tip.length=0.01, bracket.size=0.6,
                       step.increase=c(0.05), label= "p  = {signif(p,2)}") +
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold')) +
    xlab('Industry type') +
    ylab("Richness"))

stat.test <- compare_means(simpson ~ Surface, data = data2,method = "wilcox.test") %>%  
  #stat.test <- compare_means(simpson ~ Time, data = data3,method = "kruskal.test") %>%  
  mutate(y.position = c(1.05))
stat.test.signif3<-stat.test[stat.test$p.signif!="ns",]
(plot_simp<-ggplot(data2, aes(x=Surface, y=simpson, colour=Surface))+
    geom_boxplot(outlier.shape = NA, size = 1.0) +
    geom_jitter(aes(fill = Industry_type), colour = 'black', width = 0.25, pch = 21, size = 1.5) +
    scale_fill_manual(values = mycolors) +
    scale_colour_manual(values = mycolors2) +
    scale_y_continuous() +
    stat_pvalue_manual(stat.test.signif3, tip.length=0.01, bracket.size=0.6,
                       step.increase=c(0.05), label= "p  = {signif(p,2)}") +
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold')) +
    ylim(0,1.2)+
    xlab('Industry type') +
    ylab("Simpson's index"))

pdf("Alphadiv.pdf", height=7, width=12, onefile=FALSE)
ggarrange(plot_count,plot_rich, plot_simp,
          labels=c("A","",""),
          common.legend = TRUE,
          legend="bottom",
          ncol=3,
          nrow=1
)
dev.off()
############################################################## VAAAAAAMOS
data3<-cbind(meta,datafam)
data3$Sum<-rowSums(data3[,(ncol(meta)+1):ncol(data3)])

data3$Industry_type<-as.factor(data3$Industry_type)
data3$Surface<-as.factor(data3$Surface)

data3$Surface<-factor(data3$Surface, levels=c("Raw_material", "Final_product","Food_contact","Non_food_contact"))
data3$Industry_type<-factor(data3$Industry_type, levels=c("Meat","Cheese_Dairy","Processed_fish","Vegetables"))

fam1<-data3[,(ncol(meta)+1):(ncol(data3)-1)]
families2<-fam1[,order(-colMeans(fam1))]
data4<-cbind(data3[,1:ncol(meta)],families2[,1:10],rowSums(families2[,11:ncol(families2)]))
colnames(data4)[ncol(data4)]<-"Other antibiotic families"
gathered<-as.data.frame(gather(data4, key="Family", value="CPMs",-colnames(data4)[1:ncol(meta)]))
#gathered<-gathered[-which(gathered$Industry_type %in% c("Before_ripening","Intermediate_ripening")),]
levels(gathered$Family)<-factor(c(colnames(data4)[(ncol(meta)+1):ncol(data4)]))

families<-families2[which(rowSums(families2)!=0),]
data4B<-data3[,1:ncol(meta)]
data4B<-data4B[which(rowSums(families2)!=0),]
familiesB<-families*100/rowSums(families)

aggregate(gathered$CPMs, by=list(gathered$Surface, gathered$Industry_type,gathered$Family), FUN=mean )
aggregate(data4$Tetracyclines, by=list(data4$Surface), FUN=sd)


library(RColorBrewer)
colorines<-c(brewer.pal(11,"Paired"))#,brewer.pal(8,"Dark2"),"grey")


(plot_fam<-ggplot(aes(y = CPMs, x = Surface), data = gathered) +
    geom_bar(aes(fill = factor(Family, levels=levels(Family))), 
             stat = "summary", fun = "mean") +
    scale_y_continuous() +
    facet_grid(~Industry_type, scales = "free", space = "free") +
    labs(y= "ARGs (counts per million reads)")+#, fill = "Family") + #fill makes reference to the legend
    theme_bw()+
    labs(fill = "Antibiotic family") + # To change legend title
    scale_fill_manual(values = colorines) +
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
      legend.position="bottom"
    )
)

pdf("barplot.pdf", height = 7, width = 8)
plot_fam
dev.off()


(plot_fam2<-ggplot(aes(y = CPMs, x = Surface), data = gathered) +
    geom_bar(aes(fill = factor(Family, levels=levels(Family))),
             stat = "summary", fun = "mean", position="fill") +
    scale_y_continuous() +
    facet_grid(~Industry_type, scales = "free", space = "free") +
    labs(y= "ARGs (counts per million reads)")+#, fill = "Family") + #fill makes reference to the legend
    theme_bw()+
    labs(fill = "Antibiotic family") + # To change legend title
    scale_fill_manual(values = colorines) +
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="bottom"
    )
)

pdf("barplot_fill.pdf", height = 7, width = 8)
plot_fam2
dev.off()

