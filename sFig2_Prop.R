###################
## 1. Plot Prop of cells
###################

library(Seurat)
library(ggplot2)
library(dplyr)

pbmc <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
Idents(pbmc) <- pbmc$celltype.abbr
DimPlot(pbmc,label=T)
CT <- c('CSC','C','SSC','SC','IC','PVC','EC')
Color <- c("#A6CEE3","#1F78B4","#CAB2D6" ,"#33A02C","#E31A1C" ,"#FF7F00","#6A3D9A")
names(Color) <- CT

### 1.1 Number
Number_df <- pbmc@meta.data[,'celltype.abbr']
Number_df <- data.frame(table(Number_df))
colnames(Number_df) <- c('CellType','cellNumber')
# p1 <- ggplot(Number_df,aes(x=CellType,y=cellNumber,fill=CellType))+
#   geom_bar(stat="identity")+
#   scale_fill_manual(values = Color)+
#   theme_classic()+
#   theme(plot.title = element_text(hjust = 0.5))+
#   theme(axis.text = element_text(face = 'bold',size = 18,colour = 'black'),
#         axis.title = element_text(face = 'bold',size = 18,colour = 'black'))+
#   labs(x='',y='Number of cells')
# p1


library(ggplot2)
library(ggpubr)
library("ggsci")
#画下面
p0 <- ggplot(Number_df,aes(x=CellType,y=cellNumber,fill=CellType)) + 
  geom_bar(stat='identity',width=0.8,position=position_dodge(0.6)) +
  labs(x=NULL,y=NULL,fill=NULL)+    #可自定义标签名字
  coord_cartesian(ylim = c(0,5000))+#设置下面一半的值域
  scale_fill_manual(values = Color)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(face = 'bold',size = 20,colour = 'black'),
        axis.title = element_text(face = 'bold',size = 20,colour = 'black'))+
  labs(x='',y='Number of cells')

#画上面
p1 <- ggplot(Number_df,aes(x=CellType,y=cellNumber,fill=CellType)) + 
  geom_bar(stat='identity',width=0.8,position=position_dodge(0.6)) +
  labs(x=NULL,y=NULL,fill=NULL) +   #不要标签
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +     #去掉X轴和X轴的文字
  coord_cartesian(ylim = c(10000,40000)) +  #设置上面一半的值域
  scale_y_continuous(breaks = c(10000,40000,10000))+#以5为单位划分Y轴
  theme_classic()+
  theme(axis.text.y = element_text(face = 'bold',size = 20,colour = 'black'),
        text = element_text(face = 'bold',size = 20,colour = 'black'),
        axis.line.y = element_line(size=0.8, colour = "black"),
        legend.position="none")+
  theme(axis.line.x = element_line(size=0, colour = "white"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ 
  scale_fill_manual(values = Color)+
  labs(x='',y='')

#拼起来
p <- ggarrange(p1,p0,heights=c(1/5, 4/5),ncol = 1, nrow = 2,common.legend = TRUE,legend="right",align = "v")
ggsave('JingMA_NEW/res/Harmony/ALL/FIG/sFig_cellNumber.pdf',p,width = 10,height = 10)


### 1.2 all samples ratio in each cluster
phylog_df <- pbmc@meta.data[,c('batch',"celltype.abbr")]
phylog_df <- table(phylog_df$batch,phylog_df[,"celltype.abbr"])
phylog_df <- data.frame(phylog_df)
colnames(phylog_df) <- c('SampleID','CellType','Freq')
phylog_df$CellType <- factor(phylog_df$CellType,levels = rev(CT))

p2 <- ggplot(phylog_df,aes(x=SampleID,y=Freq,fill=CellType))+
  geom_col(position = "fill", width = 0.8)+
  coord_flip()+
  theme_classic()+
  theme(axis.text = element_text(face = 'bold',size = 20,colour = 'black'),
        axis.title = element_text(face = 'bold',size = 20,colour = 'black'),
        axis.line = element_line(size=1, colour = "black"),
        legend.title = element_text(size=15,face = 'bold',colour = 'black'),
        legend.text = element_text(size=15,face = 'bold',colour = 'black'))+
  labs(x='',y='Cell proportion')+theme(legend.position="right")+
  scale_fill_discrete(guide = guide_legend(reverse=TRUE))+
  scale_fill_manual(values = rev(Color))
p2
ggsave('JingMA_NEW/res/Harmony/ALL/FIG/sFig_cellProportion.pdf',p2,width = 10,height = 10)


### 1.3 Ratio
type_df <- pbmc@meta.data[,c("celltype.abbr",'type')]
type_df <- table(type_df$type,type_df[,"celltype.abbr"])
type_df <- data.frame(type_df)
colnames(type_df) <- c('SampleType','CellType','Freq')
p3 <- ggplot(type_df,aes(x=CellType,y=Freq,fill=SampleType))+
  geom_col(position = "fill", width = 0.6)+
  coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(face = 'bold',size = 18))+
  labs(x='',y='Proportion of cells')+theme(legend.position="right")
p3


### 1.4 all clusters ratio in each sample
phylog_df <- pbmc@meta.data[,c('batch',"celltype.abbr")]
phylog_df <- table(phylog_df$batch,phylog_df[,"celltype.abbr"])
phylog_df <- data.frame(phylog_df)
colnames(phylog_df) <- c('SampleID','CellType','Freq')
phylog_df$CellType <- factor(phylog_df$CellType)

p4 <- ggplot(phylog_df,aes(x=CellType,y=Freq,fill=SampleID))+
  geom_col(position = "fill", width = 0.6)+
  coord_flip()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(face = 'bold',size = 18),
        axis.line = element_line(size=1, colour = "black"))+
  labs(x='',y='Proportion of cells')+theme(legend.position="right")
p4





library(ggpubr)
pdf('/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/Barplot_Prop.pdf',width = 24,height = 4)
ggarrange(p1, p2, p3, p4,ncol = 4, nrow = 1)
dev.off()


### 1.5 Fraction
library(reshape2)
df <- table(Idents(pbmc),pbmc$batch)
df <- as.data.frame(apply(df,2,function(x) x/sum(x)))
df <- df[1:4,]
df$Cluster <- rownames(df)
df <- melt(df)
df$Type <- 'Control'
df$Type[df$variable %in% c('M1','M2','M3')] <- 'Microtia'
head(df)

library(ggpubr)
p6 <- ggbarplot(df, x = "Cluster", y = "value", 
                add = c("mean_se"),
                color = "Type",fill = "Type", palette = "npg",
                position = position_dodge(0.8))
p6

p7 <- p6+
  #+stat_compare_means(aes(group = Type), label = "p.signif",method = 'kruskal.test') 
  theme(legend.position="right",axis.text = element_text(size = 20)) +
  theme(legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 20))+
  labs(x='',y='')
p7

pdf('/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/Barplot_Prop_2.pdf',width = 8,height = 6)
print(p6)
dev.off()

##################
### 1.6 Fig3A: 只关注 CSC lineage
##################
library(Seurat)
library(ggplot2)
library(dplyr)

pbmc_C <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
pbmc_C$batch= factor(pbmc_C$batch,levels = c('C6','C4','C2','C3','C1','C5'))
print(levels(pbmc_C$batch))

phylog_p8 <- pbmc_C@meta.data[,c('batch',"celltype")]
phylog_p8 <- table(phylog_p8$batch,phylog_p8[,"celltype"])
phylog_p8 <- data.frame(phylog_p8)
colnames(phylog_p8) <- c('SampleID','CellType','Freq')
levels(phylog_p8$CellType) <- c('CSC','C0','C1','C2')

Color <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
p8 <- ggplot(phylog_p8,aes(x=SampleID,y=Freq,fill=CellType))+
  geom_col(position = "fill", width = 0.7)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size = 8,colour = 'black'))+
  theme(axis.text = element_text(size = 8,colour = "black"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.title.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5))+
  labs(x='',y='% of chondrocytes',title = 'Per sample')+
  scale_fill_manual(values = (Color))
p8


phylog_p9 <- pbmc_C@meta.data[,c('Phase',"celltype")]
phylog_p9 <- table(phylog_p9$Phase,phylog_p9[,"celltype"])
phylog_p9 <- data.frame(phylog_p9)
colnames(phylog_p9) <- c('Phase','CellType','Freq')
phylog_p9$CellType <- factor(phylog_p9$CellType)
phylog_p9$Phase <- factor(phylog_p9$Phase,levels = c('Children','Adults'))


Color <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
p9 <- ggplot(phylog_p9,aes(x=Phase,y=Freq,fill=CellType))+
  geom_col(position = "fill", width = 0.9)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size = 8,colour = 'black'))+
  theme(axis.text = element_text(size = 8,colour = "black"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.title = element_text(size=8),legend.title = element_text(size = 6,color = 'black'),
        legend.text = element_text(size = 6,color = 'black'),legend.key.size = unit(0.3,'cm'),
        legend.key.width  = unit(0.3,'cm'),legend.key.height = unit(0.3,'cm'),)+
  labs(x='',y='',title = 'Mean')+theme(legend.position="right",plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = rev(Color))+
  guides(fill = guide_legend(reverse=TRUE))
p9

library(ggpubr)
p89 <- ggarrange(p8, p9,ncol = 2, nrow = 1)
ggsave('/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/Fig3A_Barplot_PropChond.pdf',p89,
       width = 10,height = 6,units = 'cm')

