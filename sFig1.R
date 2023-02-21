###############
## plot QC
###############
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
require(RColorBrewer)
library(gridExtra)
mycol=brewer.pal(3, "Set3")

## sFig1B: QC info of sample
QC.info <- readRDS('/local/yzj/JingMA_NEW/res/QC/ALL/RDS/QC_info.RDS')
before.info <- as.data.frame(QC.info$before)
after.info <- as.data.frame(QC.info$afer)
print(QC.info)
print(sum(before.info$cellNumber))
print(sum(after.info$cellNumber))

## sFig1C: vlnbox of nFeature/nCount/percent.mt of each batch after QC
pbmc <- readRDS('/local/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
print(mean(pbmc$nCount_RNA))
print(mean(pbmc$nFeature_RNA))

print(median(pbmc$nCount_RNA))
print(median(pbmc$nFeature_RNA))

plot_QC_bybatch<- function(feature){
  p <- VlnPlot(pbmc,features =feature ,group.by = 'batch',pt.size = 0)+
    theme(axis.text=element_text(size=8,colour="black"),axis.title=element_text(size = 8,colour="black"),
          legend.position = 'none',
          plot.margin =unit(c(-0.4,0,-0.2,0.1), "cm"),
          axis.line = element_line(colour = "black", size = 0.3))+
    labs(x="",y="Numbers",title = "")
  return(p)
}

plot.nFeature_RNA <- plot_QC_bybatch('nFeature_RNA')
plot.nCount_RNA <- plot_QC_bybatch('nCount_RNA')
plot.percent.mt <-plot_QC_bybatch('percent.mt')
plot.QC <- grid.arrange(plot.nFeature_RNA,plot.nCount_RNA,plot.percent.mt,nrow = 1, ncol = 3)
ggsave('/local/yzj/JingMA_NEW/res/QC/ALL/FIG/sFig1C.pdf',plot.QC,width = 13.5,height = 4,units = 'cm')


## sFig1D: vlnbox of nFeature/nCount/percent.mt of all cells after QC
df <- data.frame(cells=colnames(pbmc),Gene='Gene',UMI='UMI',MT='MT',
                 nFeature_RNA=pbmc$nFeature_RNA,nCount_RNA=pbmc$nCount_RNA,perMT=pbmc$percent.mt)

plot.gene <- ggplot(df, aes(x=Gene, y=nFeature_RNA,fill=Gene)) + 
  geom_violin(trim=FALSE,color='black',size=0.3) + 
  geom_boxplot(width=0.3,position=position_dodge(0.9),size=0.3)+
  scale_fill_manual(values = mycol[1])+ 
  theme_classic()+ 
  theme(axis.text.x=element_blank(),axis.ticks.x =element_blank(),
        plot.margin =unit(c(0.2,0,0.2,0.1), "cm"),
        panel.grid = element_blank(),
        axis.text=element_text(size=8,colour="black"),axis.title=element_text(size = 8,colour="black"),
        axis.ticks.y =element_line(colour="black",size=0.25),
        axis.line = element_line(colour = "black", size = 0.3),legend.position = 'none')+ 
  labs(x="All",y="Numbers")

plot.UMI <- ggplot(df, aes(x=Gene, y=nCount_RNA,fill=UMI)) + 
  geom_violin(trim=FALSE,color='black',size=0.3) + 
  geom_boxplot(width=0.3,position=position_dodge(0.9),size=0.3)+
  scale_fill_manual(values = mycol[2])+ 
  theme_classic()+ 
  theme(axis.text.x=element_blank(),axis.ticks.x =element_blank(),
        plot.margin = unit(c(0.2,0,0.2,0.1), "cm"),
        panel.grid = element_blank(),
        axis.text=element_text(size=8,colour="black"),axis.title=element_text(size = 8,colour="black"),
        axis.ticks.y =element_line(colour="black",size=0.25),
        axis.line = element_line(colour = "black", size = 0.3),legend.position = 'none')+ 
  labs(x="All",y="Numbers")

plot.MT <- ggplot(df, aes(x=MT, y=perMT,fill=MT)) + 
  geom_violin(trim=FALSE,color='black',size=0.3) + 
  geom_boxplot(width=0.3,position=position_dodge(0.9),size=0.3)+
  scale_fill_manual(values = mycol[3])+ 
  theme_classic()+ 
  theme(axis.text.x=element_blank(),axis.ticks.x =element_blank(),
        plot.margin = unit(c(0.2,0,0.2,0.1), "cm"),
        panel.grid = element_blank(),
        axis.text=element_text(size=8,colour="black"),axis.title=element_text(size = 8,colour="black"),
        axis.ticks.y =element_line(colour="black",size=0.25),
        axis.line = element_line(colour = "black", size = 0.3),
        legend.position = 'none')+ 
  labs(x="All",y="Numbers")

library(gridExtra)
plot.all <- grid.arrange(plot.gene,plot.UMI, plot.MT, nrow = 1, ncol = 3)
ggsave('/local/yzj/JingMA_NEW/res/QC/ALL/FIG/sFig1D.pdf',plot.gene,width = 6,height = 4,units = 'cm')

## sFig1e: UMAP of  harmony by batch 
pbmc <- readRDS('/local/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')

plot1E <- DimPlot(pbmc, group.by='batch',label=F,pt.size = 0.02)+
  theme(axis.text = element_text(size=24,colour = 'black'),
        axis.title  = element_text(size=24,colour = 'black'),
        panel.background=element_rect(fill='transparent', color='black',size = 1),
        legend.key=element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size=24,colour = 'black'))
ggsave("/local/yzj/JingMA_NEW/res/Harmony/ALL/FIG/sFig1E.pdf",plot1E,width = 7,height = 6)

## sFig1f: UMAP of harmony by type 
plot1F <- DimPlot(pbmc, group.by = 'type',pt.size = 0.02)+
  theme(axis.text = element_text(size=24,colour = 'black'),
        axis.title = element_text(size=24,colour = 'black'),
        panel.background=element_rect(fill='transparent', color='black',size = 1),
        legend.key=element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size=24,colour = 'black'))
ggsave("/local/yzj/JingMA_NEW/res/Harmony/ALL/FIG/sFig1F.pdf",plot1F,width = 7,height = 5.5)


## sFig1g. UMAP by cluster resol=1.2
plot1G <- DimPlot(pbmc, group.by='RNA_snn_res.1.2',label=T,pt.size = 0.01,label.size = 8)+
  theme(axis.text = element_text(size=24,colour = 'black'),
        axis.title = element_text(size=24,colour = 'black'),
        panel.background=element_rect(fill='transparent', color='black',size = 2),
        legend.key=element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size=24,colour = 'black'))
ggsave("/local/yzj/JingMA_NEW/res/Harmony/ALL/FIG/sFig1G.pdf",plot1G,
       width = 7,height = 6)


## sFigh. Feature plot of mk genes
marker.genes <- rev(c('CDH5','CLDN5','PDGFRB','ACTA2','PTPRC','HLA-DRA','LUM','CYTL1','ELN','COL2A1','HES1','EGR1','VCAN','COL1A1','ACAN','COL9A2'))

Idents(pbmc) <- pbmc$celltype
subpbmc <- subset(pbmc,downsample=10000)

plot.lst <- list()
for(i in 1:length(marker.genes)){
  gene <- marker.genes[i]
  plot.lst[[i]] <- FeaturePlot(subpbmc,gene,pt.size = 0.1,order = F)+ 
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.1,0.1,0.1,0), "cm"),
          axis.text = element_blank(),axis.title = element_blank(), plot.title = element_text(hjust = 0.5),
          axis.ticks = element_blank(),legend.position = 'none')+labs(title=gene)
}

library(cowplot)
p <-  plot_grid(plot.lst[[1]],plot.lst[[2]],plot.lst[[3]],plot.lst[[4]],
                plot.lst[[5]],plot.lst[[6]],plot.lst[[7]],plot.lst[[8]],
                plot.lst[[9]],plot.lst[[10]],plot.lst[[11]],plot.lst[[12]],
                plot.lst[[13]],plot.lst[[14]],plot.lst[[15]],plot.lst[[16]],ncol = 4)
ggsave('/local/yzj/JingMA_NEW/res/Harmony/ALL/FIG/sFig1_FeaturePlot.pdf',p,width = 21,height = 21,units = 'cm')

## sFig1I/J. cell cycle
library(Seurat)
library(ggplot2)
pbmc <- readRDS('/local/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
Idents(pbmc) <- pbmc$celltype

pbmc <- CellCycleScoring(pbmc,g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
data=as.data.frame(pbmc@reductions$harmony@cell.embeddings[,c(2,1)])
data$CellType <- pbmc$celltype
data$S.Score <- pbmc$S.Score
data$G2M.Score <- pbmc$G2M.Score
data$Score <- (data$S.Score+data$G2M.Score)
colors <-  c("#EE9572","#B2DF8A" )
names(colors) <- c('CSC', 'SSC')

data.CSC <- data[data$CellType %in% c('CSC','SSC'),]
plot.density <- ggplot(data.CSC,aes(x=Score,color=CellType,fill=CellType))+
  geom_density(size=0.3, alpha=.4)+  
  theme_classic()+
  scale_color_manual(breaks = names(colors), values = colors) + 
  scale_color_discrete(breaks = c('CSC', 'SSC'),
                       labels = c('CSPC', 'SSC')) + 
  scale_fill_discrete(breaks = c('CSC', 'SSC'),
                      labels = c('CSPC', 'SSC')) + 
  theme(panel.background=element_rect(fill='transparent', color='black',size = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.1),plot.margin = unit(c(0.1,0.1,0.1,0.1),'cm')) + 
  theme(axis.title = element_text(size = 6,colour = 'black'),
        legend.position = 'none',
        axis.line = element_line(arrow = arrow(length = unit(0.2, 'cm'))),
        axis.text = element_text(size = 6,colour = 'black'), 
        axis.ticks = element_line(size = 0.3))+
  labs(x='Cell cyle score',y='Density')
plot.density
ggsave('JingMA_NEW/res/Harmony/ALL/FIG/sFig1_density.pdf',plot.density,width = 6,height = 3,units = 'cm')


pbmc <- readRDS('/local/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
Idents(pbmc) <- pbmc$celltype

get_cc <- function(pbmc,CT){
  Idents(pbmc) <- pbmc$celltype
  subpbmc <- subset(pbmc,idents=CT)
  subpbmc <- CellCycleScoring(subpbmc,g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
  data=as.data.frame(subpbmc@reductions$harmony@cell.embeddings[,c(2,1)])
  data$CellType <- subpbmc$celltype
  data$S.Score <- subpbmc$S.Score
  data$G2M.Score <- subpbmc$G2M.Score
  data$Score <- (data$S.Score+data$G2M.Score)/2
  return(data)
}

cc_CSC <- get_cc(pbmc,'CSC')
cc_SSC <- get_cc(pbmc,'SSC')
cc_C <- get_cc(pbmc,'C')
cc_SC <- get_cc(pbmc,'SC')

colors <-  c("red","blue","green","yellow")
names(colors) <- c('CSC', 'SSC','C','SC')

plot.df <- rbind(rbind(rbind(cc_CSC,cc_SSC),cc_C),cc_SC)
plot.df$CT <- 'non-SC'
plot.df$CT[plot.df$ce]
plot.density <- ggplot(plot.df,aes(x=Score,color=CellType,fill=CellType))+
  geom_density(size=0.3, alpha=.4)+  
  theme_classic()+
  scale_color_manual(breaks = names(colors), values = colors) +
  theme(panel.background=element_rect(fill='transparent', color='black',size = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.1),plot.margin = unit(c(0.1,0.1,0.1,0.1),'cm')) + 
  theme(axis.title = element_text(size = 8,colour = 'black'),
        axis.line = element_line(arrow = arrow(length = unit(0.2, 'cm'))),
        axis.text = element_text(size = 8,colour = 'black'), 
        axis.ticks = element_line(size = 0.3))+
  labs(x='Cell cyle score',y='Density')
plot.density

plot.cc <- ggplot(data,aes(S.Score,G2M.Score))+
  geom_point(size=0.5,colour="#EE9572")+xlim(c(-0.1,0.2))+ylim(c(-0.1,0.2))+
  theme(panel.background=element_rect(fill='transparent', color='black',size = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  theme_classic()+
  theme(axis.title = element_text(size = 6,colour = 'black'),
        legend.position = 'none',
        axis.line = element_line(arrow = arrow(length = unit(0.2, 'cm'))),
        axis.text = element_text(size = 6,colour = 'black'), 
        axis.ticks = element_line(size = 0.3))+
  labs(x='S.Score',y='G2M.Score')
plot.cc
ggsave('JingMA_NEW/res/Harmony/ALL/FIG/sFig1_cc.pdf',plot.cc,width = 6,height = 6,units = 'cm')

# MKI67, TOP2A
pbmc <- readRDS('/local/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
marker.genes <- c('MKI67','TOP2A','CDK1','MCM2')
Idents(pbmc) <- pbmc$celltype
subpbmc <- subset(pbmc,downsample=10000)

plot.lst <- list()
for(i in 1:length(marker.genes)){
  gene <- marker.genes[i]
  plot.lst[[i]] <- FeaturePlot(subpbmc,gene,pt.size = 0.001,order = F)+ 
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.1,0.1,0.1,0), "cm"),
          axis.text = element_blank(),axis.title = element_blank(), plot.title = element_text(hjust = 0.5,size=24),
          axis.ticks = element_blank(),legend.position = 'none')+labs(title=gene)
}
library(cowplot)
p <-  plot_grid(plot.lst[[1]],plot.lst[[2]],plot.lst[[3]],plot.lst[[4]],ncol = 2)
ggsave('/local/yzj/JingMA_NEW/res/Harmony/ALL/FIG/sFig1_proliferMK.pdf',p,width = 5,height = 6)

