library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratWrappers)
library(harmony)

pbmc <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')

###################
## 1. clustering
###################
pbmc_harmony <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/PBMC_harmony.RDS')
dimN=50
resol=1.2
pbmc <- FindNeighbors(pbmc, dims = 1:dimN,reduction = 'harmony',features = VariableFeatures(object = pbmc))
pbmc <- FindClusters(pbmc,resolution = resol,verbose = 0,algorithm = 2)

pdf(paste("/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/UMAP_resol",resol,'.pdf',sep=''),width = 7,height = 6)
DimPlot(pbmc, group.by='RNA_snn_res.1.2',label=T,pt.size = 1,label.size = 5)+
  theme(axis.text = element_text(size=15),
        panel.background=element_rect(fill='transparent', color='black',size = 1.5),
        legend.key=element_rect(fill='transparent', color='transparent'))
dev.off()

###################
## 2. Naming
###################
Idents(pbmc) <- pbmc$celltype
require('RColorBrewer')
CT <- c('CSC','C','SSC','SC','IC','PVC','EC')
Color <- c("#A6CEE3","#1F78B4","#CAB2D6" ,"#33A02C","#E31A1C" ,"#FF7F00","#6A3D9A")

names(Color) <- CT
pbmc$color <- Color[pbmc$celltype]
pbmc$celltype <- factor(pbmc$celltype,levels = CT)
Idents(pbmc) <- pbmc$celltype
saveRDS(pbmc,'/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')


###################
## 3.合并细胞亚型，后面分析产生的
###################
pbmc <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
pbmc_chond <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
pbmc_stroma <- readRDS('/home/disk/drizzle/wgk/data/AllSample_2_merge/stroma_lineage/seurat_celltype.Rdata')

pbmc_chond$celltype <- as.character(pbmc_chond$celltype)
pbmc_chond$celltype[pbmc_chond$celltype=='TC'] <- 'C0'
pbmc_chond$celltype <- factor(pbmc_chond$celltype,levels = c('CSC','C0','C1','C2'))
DimPlot(pbmc_chond,group.by = 'celltype')

C0_cells <- colnames(pbmc_chond)[pbmc_chond$celltype=='C0']
C1_cells <- colnames(pbmc_chond)[pbmc_chond$celltype=='C1']
C2_cells <- colnames(pbmc_chond)[pbmc_chond$celltype=='C2']

SSC_cells <- colnames(pbmc_stroma)[pbmc_stroma$celltype=='SSC']
SC1_cells <- colnames(pbmc_stroma)[pbmc_stroma$celltype=='SC1']
SC2_cells <- colnames(pbmc_stroma)[pbmc_stroma$celltype=='SC2']

setdiff(colnames(pbmc)[pbmc$celltype %in% c('SSC','SC')],colnames(pbmc_stroma))

unique(pbmc$celltype)
pbmc$celltype[colnames(pbmc) %in% C0_cells] <- 'C0'
pbmc$celltype[colnames(pbmc) %in% C1_cells] <- 'C1'
pbmc$celltype[colnames(pbmc) %in% C2_cells] <- 'C2'
pbmc$celltype[colnames(pbmc) %in% SSC_cells] <- 'SSC'
pbmc$celltype[colnames(pbmc) %in% SC1_cells] <- 'SC1'
pbmc$celltype[colnames(pbmc) %in% SC2_cells] <- 'SC2'

pbmc$celltype <- factor(pbmc$celltype,levels = c('CSC','C0','C1','C2','SSC','SC1','SC2','IC','PVC','EC'))
DimPlot(pbmc,group.by = 'celltype')

saveRDS(pbmc,'/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
saveRDS(pbmc_chond,'/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
#####


################
## Fig1B
################
CT <- c('CSC','C','SSC','SC','IC','PVC','EC')
Color <- c("#A6CEE3","#1F78B4","#CAB2D6" ,"#33A02C","#E31A1C" ,"#FF7F00","#6A3D9A")
names(Color) <- CT

library(ggplot2)
p <- DimPlot(pbmc, group.by='celltype',label=F,pt.size = 0.1,
        cols = Color)+
  theme(axis.text = element_text(size=8,colour = 'black'),
        axis.title = element_text(size=8,colour = 'black'),
        panel.background=element_rect(fill='transparent', color='black',size = 1),
        legend.key=element_blank(),legend.key.size  = unit(0.5,'cm'),
        legend.key.height = unit(0.5,'cm'),legend.key.width  = unit(0.5,'cm'),
        legend.text = element_text(size=8,colour = 'black'),
        plot.margin = unit(c(0,0,-0.1,-0.1),'cm'))
ggsave("/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/Fig1B_UMAP_CellType.pdf",p,width = 10,height = 8,units = 'cm')

pdf("/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/UMAP_CellType_noLable.pdf",width = 5,height = 4)
DimPlot(pbmc, group.by='celltype', label=F,pt.size = 1,cols = Color)+
  theme(axis.text = element_text(size=20),
        panel.background=element_rect(fill='transparent', color='black',size = 1),
        legend.key=element_rect(fill='transparent', color='transparent'))
dev.off()
########

########
pdf("/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/UMAP_Harmony_batch.pdf",width = 12,height = 12)
DimPlot(pbmc, group.by='batch', label=F,split.by = 'batch',pt.size = 1,ncol = 3)+
  theme(axis.text = element_text(size=20),
        panel.background=element_rect(fill='transparent', color='black',size = 1.5),
        legend.key=element_rect(fill='transparent', color='transparent'))
dev.off()

pdf("/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/UMAP_Harmony_type_1.pdf",width = 12,height = 6)
DimPlot(pbmc, group.by='batch', label=F,split.by = 'type',pt.size = 0.2)+
  theme(axis.text = element_text(size=20),
        panel.background=element_rect(fill='transparent', color='black',size = 1.5),
        legend.key=element_rect(fill='transparent', color='transparent'))
dev.off()

pdf("/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/UMAP_Celltype_Type.pdf",width = 12,height = 6)
DimPlot(pbmc, group.by = 'celltype',split.by = 'type',pt.size = 0.3,label.size = 5,
        cols = Color)+
  theme(axis.text = element_text(size=20),
        panel.background=element_rect(fill='transparent', color='black',size = 1.5),
        legend.key=element_rect(fill='transparent', color='transparent'))
dev.off()
########

########
## celltype表明general的细胞类型，celltype表明详细的chond的细胞类型
########
pbmc <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
pbmc_chond <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')

Idents(pbmc) <- pbmc$celltype
subpbmc <- subset(pbmc,idents = c('CSC','SSC','SC','IC','PVC','EC'))


Cell2CT <- c(as.character(subpbmc$celltype),as.character(pbmc_chond$celltype))
names(Cell2CT) <- c(colnames(subpbmc),colnames(pbmc_chond))
unique(pbmc$celltype)
pbmc$celltype <- Cell2CT[colnames(pbmc)]
DimPlot(pbmc,group.by = 'celltype')
saveRDS(pbmc,'/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')

########
## Fig1C violin plot
########
library(ggplot2)
library(Seurat)

pbmc <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
marker.genes <- c('CDH5','CLDN5','PDGFRB','ACTA2','PTPRC','HLA-DRA','COL1A1','LUM','VCAN','ACAN','COL9A2','CYTL1','ELN','COL2A1','EGR1','HES1')
df.gene <- data.frame(stringsAsFactors = F)
for (gene in marker.genes) {
  df.sub <- data.frame(expvalue = pbmc@assays$RNA@data[gene,],
                       gene = rep(gene, ncol(pbmc@assays$RNA@data)),
                       celltype = pbmc$celltype)
  df.gene <- rbind(df.gene, df.sub)
}
df.plot <- df.gene
df.plot$gene <- factor(df.gene$gene, levels = marker.genes)
df.plot$celltype <- factor(df.gene$celltype, 
                           levels = c('CSC', 'C', 'SSC', 'SC','IC','PVC','EC'))
color.cell <- c("#A6CEE3" ,"#1F78B4","#CAB2D6","#33A02C","#E31A1C","#FF7F00" ,"#6A3D9A")
plot.vln <- 
  ggplot(data = df.plot, aes(x = gene, y = expvalue, color = celltype, fill = celltype)) + 
  geom_violin(trim = T, scale = 'width') + 
  scale_color_manual(labels = c('CSC','C', 'SSC','SC','IC','PVC','EC'),values = color.cell) + 
  scale_fill_manual(labels = c('CSC','C', 'SSC','SC','IC','PVC','EC'),values = color.cell) + 
  facet_grid( ~ celltype) + 
  theme_classic() + coord_flip() +
  stat_summary(fun= mean, geom = "point",shape = 23, size = 2, color = "black") + 
  labs(x = 'Gene', y = 'Expression Level') + 
  theme(axis.text.y = element_text(size = 8, color = "black", face = 'italic'), 
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8,color = "black"), 
    strip.text.x = element_text(size = 8, color = "black"), legend.position = 'none')
plot.vln
ggsave("/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/Fig1C_Vln.pdf",plot.vln,
       height = 8, width = 16, units = 'cm')

########
## Fig1E heatmap plot
########
library(ggplot2)
library(Seurat)

pbmc <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
Idents(pbmc) <- pbmc$celltype
MK <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype.RDS')

CT <- rev(c( 'EC', 'PVC','IC', 'SC','SSC', 'C', 'CSC'))
color.cell <- c("#A6CEE3" ,"#1F78B4","#CAB2D6","#33A02C","#E31A1C","#FF7F00" ,"#6A3D9A")
sel.genes <- c()
for (cell in CT) {
  sub.markers <- MK[[cell]]
  sub.markers$diff.pct <- sub.markers$pct.1 - sub.markers$pct.2
  sub.markers <- sub.markers[sub.markers$p_val_adj < 0.01 & sub.markers$avg_logFC > 0.5 & 
                               sub.markers$diff.pct > 0.1, ]
  sub.markers <- sub.markers[order(sub.markers$avg_logFC, decreasing = T),]
  sel.genes <- c(sel.genes, rownames(sub.markers)[1:min(nrow(sub.markers), 50)])
}
sel.genes <- unique(sel.genes)

tmp <- AverageExpression(pbmc, return.seurat = TRUE)
tmp@active.ident <- factor(tmp@active.ident, levels = (CT))

plot.heatmap <- 
  DoHeatmap(tmp, features = sel.genes,draw.lines = FALSE, 
            label = T, group.colors = color.cell, size = 3,
            slot = 'scale.data', disp.min = -2, disp.max = 2) +
  guides(color = F) + 
  labs(fill = 'Scaled Expression') + 
  theme(axis.text = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.position = 'bottom', legend.direction = 'horizontal',
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.2,"cm") ,
        plot.margin = unit(c(0.5,0,0,0), "cm"))+
  scale_fill_gradientn(colors = c("navy", "white", "firebrick3"), 
                       breaks = c(-1, 0, 1)) 

ggsave('JingMA_NEW/res/Harmony/ALL/FIG/Fig1E_DoHeatmap_markers.pdf',plot.heatmap,
       height = 10, width = 4, units = 'cm')




########
## Fig1E GO enrichemnt bar plot
########
library(xlsx)
# plot
list.go.BP <- list()
CT=c('CSC', 'C', 'SSC', 'SC', 'IC', 'PVC', 'EC')
for(ct in CT){
  list.go.BP[[ct]] <- read.xlsx('JingMA_NEW/res/compCT/All/FC1.5_Qvalue0.05/GO_BP_all.xlsx',sheetName = ct) 
}

GO.BP <- c('extracellular matrix organization',
           'response to BMP',
           'chondrocyte differentiation',
           'cellular response to transforming growth factor beta stimulus',
           'cell cycle arrest',
           'response to metal ion',
           'extracellular structure organization',
           'positive regulation of cell-substrate adhesion',
           'extracellular matrix organization',
           'connective tissue development',
           'cellular response to transforming growth factor beta stimulus',
           'cell cycle arrest',
           'BMP signaling pathway',
           'extracellular matrix organization',
           'collagen fibril organization',
           'regulation of smooth muscle cell proliferation',
           'regulation of inflammatory response',
           'regulation of vasculature development',
           'neutrophil activation involved in immune response',
           'leukocyte migration',
           'antigen processing and presentation',
           'muscle system process',
           'muscle cell proliferation',
           'endothelium development',
           'positive regulation of angiogenesis'
)
GO.cell <- c('CSC', 'CSC', 'CSC', 'CSC', 'CSC',  
             'C', 'C', 'C',
             'SSC', 'SSC', 'SSC', 'SSC', 'SSC',
             'SC', 'SC', 'SC', 'SC','SC',
             'IC', 'IC', 'IC','PVC', 'PVC', 'EC', 'EC')

sort.cells <- CT
color.cell <- c("#A6CEE3","#1F78B4","#CAB2D6" ,"#33A02C","#E31A1C" ,"#FF7F00","#6A3D9A")

df.plot <- data.frame(stringsAsFactors = F)
for (i in 1:length(GO.BP)) {
  print(i)
  GO.term <- GO.BP[i]
  cell <- GO.cell[i]
  sub.BP <- list.go.BP[[cell]]
  df.plot[i, 'GO'] <- GO.term
  df.plot[i, 'Cell'] <- cell
  df.plot[i, 'pvalue'] <- sub.BP[
    sub.BP$Description == GO.term, 'pvalue']
  df.plot[i, 'name'] <- paste0(cell, '_', GO.term)
}
df.plot$log10Pval <- -log10(df.plot$pvalue)
df.plot$Cell <- factor(df.plot$Cell, levels = sort.cells)
df.plot$name <- factor(df.plot$name, levels = rev(paste(GO.cell, GO.BP, sep = '_')))
df.plot$log10Pval[df.plot$log10Pval > 20] <- 20

plot.bar <- 
  ggplot(data = df.plot, aes(x = name, y = log10Pval, color = Cell, fill = Cell)) + 
  geom_bar(stat = 'identity',width=0.8) + 
  theme_classic() + coord_flip() +
  scale_x_discrete(breaks = paste(GO.cell, GO.BP, sep = '_'), 
                   labels = GO.BP,
                   position = "top") + 
  scale_color_manual(values = color.cell) +
  scale_fill_manual(values = color.cell) +
  labs(x = 'Enriched GO terms', y = expression(paste("-log"[10], "(", italic("P"), "-value)"))) +
  theme(axis.title = element_text(size = 6, color = 'black'), 
        axis.text.y = element_text(size = 6, color = 'black'), 
        axis.text.x = element_text(size = 6,  color = 'black'),
        legend.position = 'none')

ggsave(filename = 'JingMA_NEW/res/Harmony/ALL/FIG/Fig1E_heatmap_bar.pdf',plot.bar,height = 8.5, width = 10,units = 'cm')
