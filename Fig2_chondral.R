library(Seurat)
library(harmony)
library(ggplot2)
require("RColorBrewer")

pbmc <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
pbmc_chond <- subset(pbmc, subset = celltype %in% c('CSC', 'C'))

pbmc_chond <- NormalizeData(pbmc_chond)
pbmc_chond <- FindVariableFeatures(pbmc_chond, nfeatures = 3000)
pbmc_chond <- ScaleData(pbmc_chond, split.by = "batch")
pbmc_chond <- RunPCA(pbmc_chond, verbose = F, npcs = 100)
pbmc_chond <- RunHarmony(pbmc_chond, "batch", reduction.save = "harmony2")
pbmc_chond <- RunUMAP(pbmc_chond, reduction = "harmony2", 
                       dims = 1:80, n_neighbors = 50, min.dist = 1)

DimPlot(pbmc_chond, group.by = "celltype", label = T, reduction = 'umap', pt.size = 1)

pbmc_chond <- FindNeighbors(pbmc_chond, reduction = "harmony", dims = 1:80)
pbmc_chond <- FindClusters(pbmc_chond, resolution = 0.3)

saveRDS(pbmc_chond,'JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')

########
## Fig2A. UMAP
########
pbmc_chond <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
pbmc_chond$celltype <- as.character(pbmc_chond$celltype)
pbmc_chond$celltype <- factor(pbmc_chond$celltype, 
                               levels = c('CSC', 'C0', 'C1', 'C2'))
colors <-  c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
names(colors) <- c('CSC', 'C0', 'C1', 'C2')

data=as.data.frame(pbmc_chond@reductions$umap@cell.embeddings)
data$CellType <- pbmc_chond$celltype
plot.umap <- ggplot(data = data,aes(x=UMAP_1,y=UMAP_2,color=CellType))+
  geom_point(size=0.1)+
  theme_bw()+
  theme(axis.title = element_text(size = 7,colour = 'black'),
        panel.background=element_rect(fill='transparent', color='black',size = 0.3),
        legend.key=element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size = 7,colour = 'black'),
        legend.title = element_text(size = 7,colour = 'black'),
        legend.key.size = unit(0.1, "cm"),legend.key.width = unit(0.1,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank()) + 
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(values = colors)
ggsave( 'JingMA_NEW/res/Harmony/ALL/FIG/Fig2A_chon_lineage_umap.pdf',plot.umap,
        height = 4, width = 6, units = 'cm')




########
## Fig2B. 2 PCA
########
pc.stdev <- pbmc_chond@reductions$harmony2@stdev[1:10]
pc.stdev[1]/sum(pc.stdev)
pc.stdev[2]/sum(pc.stdev)

data=as.data.frame(pbmc_chond@reductions$harmony@cell.embeddings[,c(2,1)])
data$CellType <- pbmc_chond$celltype
plot.lineage <- ggplot(data = data,aes(x=harmony_2,y=harmony_1,color=CellType))+
  geom_point(size=0.1)+
  labs(x = 'PC2 (15%)', y = 'PC1 (17%)') + 
  theme(panel.background=element_rect(fill='transparent', color='black',size = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_color_manual(values = colors)+
  theme_classic()+
  theme(axis.title = element_text(size = 7,colour = 'black'),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.position = 'none',
        axis.line = element_line(arrow = arrow(length = unit(0.2, 'cm'))))
ggsave('JingMA_NEW/res/Harmony/ALL/FIG/Fig2B_chon_lineage_harmony.pdf',plot.lineage,
       height = 4, width = 6, units = 'cm')


########
## Fig2C: pick genes
########
# PC2 <- pbmc_chond.bak@reductions$pca@cell.embeddings[, 'PC_2']
sample.cells <- sample(colnames(pbmc_chond),10000)
sample.pbmc_chond <- subset(pbmc_chond,cells = sample.cells)
Harmony2 <- sample.pbmc_chond@reductions$harmony@cell.embeddings[, 'harmony_2']
mat.gene <- sample.pbmc_chond@assays$RNA@data
df.pc.gene <- data.frame(t(as.matrix(mat.gene)))
df.pc.gene$Harmony2 <- Harmony2
df.pc.gene$celltype <- sample.pbmc_chond$celltype
df.pc.gene$status <- sample.pbmc_chond$type

genes <- c('EGR1', 'HES1','FRZB', 'CTGF','ACAN', 'COL9A2',  'MMP3',  'IL8')
list.plot <- list()
for (gene in genes) {
  df.plot <- df.pc.gene[, c('Harmony2', 'celltype', gene)]
  names(df.plot) <- c('Harmony2', 'celltype', 'gene')
  p.gene <- 
    ggplot(data = df.plot, aes(x = Harmony2, y = gene)) + 
    geom_point(aes(color = celltype), size = 0.1) + 
    scale_color_manual(labels = c('CSC', 'C0', 'C1', 'C2'),values = colors) + 
    xlim(-30, 10) + 
    geom_smooth(color = '#696969',size = 0.5) + 
    labs(x = '', y = '') + theme_bw()+
    theme(panel.background=element_rect(color='black',size =0.4,linetype ="solid"),panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),plot.margin = unit(c(0.1,0.1,-0.5,-0.5), "cm"),
          legend.position = 'none') +
    annotate('text', label = gene, x = 10, y = max(df.plot$gene), 
             hjust = 1, vjust = 1, size = 3)
  list.plot[[gene]] <- p.gene
}

library(gridExtra)
p.merge <- marrangeGrob(list.plot, ncol = 4, nrow = 2, top = NULL)
ggsave(plot = p.merge, path = 'JingMA_NEW/res/Harmony/ALL/FIG/', 
       filename = 'Fig2C_chon_genes.pdf',
       height = 4, width = 12, units = 'cm')


grid.arrange(plot3, plot2, plot4, plot1, nrow = 2, ncol = 2)

########
## Fig2D: Heatmap
########
# cell marker
clusters <- unique(pbmc_chond$celltype)
list.marker.all <- list()
for (cluster in clusters) {
  sub.markers <- FindMarkers(pbmc_chond, ident.1 = cluster, group.by = 'celltype',
                             logfc.threshold = 0.3, min.diff.pct = 0.05, only.pos = T)
  list.marker.all[[cluster]] <- sub.markers
}
file.marker <- paste0(path.lineage, 'marker_genes.Rdata')
saveRDS(list.marker.all, file.marker)
list.marker.all <- readRDS(file.marker)

# heatmap
list.marker.all <- readRDS('/home/disk/drizzle/wgk/data/AllSample_2_merge/chon_lineage/marker_genes.Rdata')
names(list.marker.all)[4] <- 'C0'
sort.cells <- c('CSC', 'C0', 'C1', 'C2')

sel.genes <- c()
for (cell in sort.cells) {
  sub.markers <- list.marker.all[[cell]]
  sub.markers$diff.pct <- sub.markers$pct.1 - sub.markers$pct.2
  sub.markers <- sub.markers[sub.markers$p_val_adj < 0.01 & sub.markers$avg_logFC > 0.5 & 
                               sub.markers$diff.pct > 0.05, ]
  sub.markers <- sub.markers[order(sub.markers$avg_logFC, decreasing = T),]
  sel.genes <- c(sel.genes, rownames(sub.markers)[1:min(nrow(sub.markers), 50)])
}
sel.genes <- intersect(unique(sel.genes), rownames(pbmc_chond@assays$RNA@scale.data))

pbmc_chond.sample <- subset(pbmc_chond, 
                             cells = sample(rownames(pbmc_chond@meta.data), 5000))
plot.heatmap <- 
  DoHeatmap(pbmc_chond.sample, features = sel.genes,
            draw.lines = T, lines.width = 20,
            group.by = 'celltype',
            label = T, group.colors = colors, size = 3,
            slot = 'scale.data', disp.min = -1.5, disp.max = 1.5) +
  guides(color = F) + 
  labs(fill = 'Scaled Expression') + 
  theme(axis.text = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.2,"cm") ,
        legend.position = 'bottom', legend.direction = 'horizontal',
        plot.margin = unit(c(0.5,0,0,0), "cm"))+
  scale_fill_gradientn(colors = c("navy", "white", "firebrick3"), 
                       breaks = c(-1, 0, 1)) 
ggsave(filename = '/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/Fig2D_DoHeatmap_markers_sample.pdf',plot.heatmap,
       height = 9, width = 9, units = 'cm')


########
## Fig2D: GO terms
########
# GO
fc.cutoff <- 0.3
pct.cutoff <- 0.05
clusters <- unique(pbmc_chond$celltype)
list.marker.all <- list()
for (cluster in clusters) {
  sub.markers <- FindMarkers(pbmc_chond, ident.1 = cluster, group.by = 'celltype',
                             logfc.threshold = fc.cutoff, min.diff.pct = pct.cutoff, only.pos = T)
  sub.markers <- sub.markers[sub.markers$p_val_adj < 0.01,]
  list.marker.all[[cluster]] <- sub.markers
}

path.cutoff <- paste0(path.lineage, 'cutoff_', fc.cutoff, '_', pct.cutoff, '/')
if (!file.exists(path.cutoff)) {
  dir.create(path.cutoff)
}
file.marker.all <- paste0(path.cutoff, 'marker_go.Rdata')
saveRDS(list.marker.all, file = file.marker.all)
list.marker.all <- readRDS(file.marker.all)

# enrich GO
library(dplyr)
df.gene_id <- read.delim('/home/disk/drizzle/wgk/ncbi/Homo_sapiens.gene_info')
df.gene_id <- df.gene_id %>% distinct(Symbol, .keep_all = T)
df.symbol2geneid <- data.frame(Gene_id = df.gene_id$GeneID, row.names = df.gene_id$Symbol)
df.geneid2symbol <- data.frame(Symbol = df.gene_id$Symbol, row.names = df.gene_id$GeneID)
library(clusterProfiler)
library(org.Hs.eg.db)
all.genes <- rownames(pbmc_chond@assays$RNA@data)
use.genes <- intersect(keys(org.Hs.eg.db, keytype = "SYMBOL"), all.genes)
list.go.BP <- list()
list.go.MF <- list()
list.go.CC <- list()
list.kegg <- list()
for (cell in clusters) {
  sub.markers <- list.marker.all[[as.character(cell)]]
  genes.input <- intersect(row.names(sub.markers), use.genes)
  egmt <- enrichGO(gene = genes.input, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                   universe = use.genes, pvalueCutoff = 0.1, ont = 'BP')
  # res.egmt <- simplify(egmt)@result
  res.egmt <- egmt@result
  list.go.BP[[as.character(cell)]] <- res.egmt
  # egmt <- enrichGO(gene = genes.input, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
  #                  universe = use.genes, pvalueCutoff = 0.5, ont = 'MF')
  # # res.egmt <- simplify(egmt)@result
  # res.egmt <- egmt@result
  # list.go.MF[[as.character(cell)]] <- res.egmt
  # egmt <- enrichGO(gene = genes.input, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
  #                  universe = use.genes, pvalueCutoff = 0.5, ont = 'CC')
  # # res.egmt <- simplify(egmt)@result
  # res.egmt <- egmt@result
  # list.go.CC[[as.character(cell)]] <- res.egmt
  # res.kegg <- enrichKEGG(gene = df.symbol2geneid[genes.input, 'Gene_id'], 
  #                        universe = as.character(df.symbol2geneid[use.genes, 'Gene_id']), 
  #                        pvalueCutoff = 0.5, keyType = 'kegg')
  # list.kegg[[as.character(cell)]] <- res.kegg@result
}

file.go.BP <- paste0(path.cutoff, 'GO_BP.Rdata')
saveRDS(list.go.BP, file = file.go.BP)

## 画图
list.go.BP <- readRDS('/home/disk/drizzle/wgk/data/AllSample_2_merge/chon_lineage/cutoff_0.3_0.05/GO_BP.Rdata')

list.sel.GO <- list()
list.sel.GO$CSC <- c('cartilage development','extracellular matrix organization', 'chondrocyte differentiation',
                     'BMP signaling pathway', 'cell cycle arrest', 'Wnt signaling pathway')
list.sel.GO$TC <- c('cellular response to zinc ion', 'extracellular matrix organization',
                    'cellular response to copper ion', 'cartilage development',
                    'chondrocyte differentiation')
list.sel.GO$C1 <- c('extracellular matrix organization', 'nitric oxide biosynthetic process',
                    'cellular response to interferon-gamma', 'lymphocyte migration',
                    'cellular response to tumor necrosis factor', 
                    'regulation of cell shape')
list.sel.GO$C2 <- c('interleukin-1-mediated signaling pathway', 'neutrophil migration', 
                    'response to oxidative stress', 'extracellular matrix organization',
                    'nitric oxide biosynthetic process', 'extracellular matrix disassembly')
# barplot
sort.cells <- c('CSC', 'TC', 'C1', 'C2')
colors <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
i = 0
for (cell in sort.cells) {
  i = i + 1
  sub.go <- list.go.BP[[cell]]
  sel.go.term <- list.sel.GO[[cell]]
  sel.go <- sub.go[sub.go$Description %in% sel.go.term, 
                   c('Description', 'pvalue')]
  sel.go$log10Pval <- -log10(sel.go$pvalue)
  sel.go$Description <- factor(sel.go$Description, levels = rev(sel.go.term))
  p <- ggplot(sel.go, aes(x = Description, y = log10Pval)) + 
    geom_bar(stat = 'identity', color = colors[i], fill = colors[i],width = 0.8) + 
    theme_classic() + coord_flip() +
    labs(x='',y = expression(paste("-log"[10], "(", italic("P"), "-value)"))) +
    theme(axis.title = element_text(size = 8, colour = 'black'), 
          axis.text.y = element_text(size = 8, colour = 'black'), 
          axis.text.x = element_text(size = 8, colour = 'black'))
  ggsave(paste('/home/yzj/JingMA_NEW/res/Harmony/ALL/FIG/Fig2D_',cell,'_bar.pdf',sep=''),p,
         height = 4, width = 8, units = 'cm')
}



#############
## 补充材料
##############
library(Seurat)
library(ggplot2)

pbmc_chond <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')

pdf('JingMA_NEW/res/Harmony/ALL/FIG/sFig_chon_lineage_umap.pdf',height = 3, width = 4)
DimPlot(pbmc_chond, group.by = "RNA_snn_res.0.3", label = T, reduction = 'umap', pt.size = 0.3)
dev.off()

pdf('JingMA_NEW/res/Harmony/ALL/FIG/sFig_chon_lineage_type.pdf',height = 3.5, width = 5)
DimPlot(pbmc_chond, group.by = "type", label = F, reduction = 'umap', pt.size = 0.01)
dev.off()


########
## sFig: pick genes
########
colors <-  c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
names(colors) <- c('CSC', 'C0', 'C1', 'C2')

pbmc_chond <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
Idents(pbmc_chond) <- pbmc_chond$celltype
sample.cells <- sample(colnames(pbmc_chond),10000)
sample.pbmc_chond <- subset(pbmc_chond,cells = sample.cells)
Harmony2 <- sample.pbmc_chond@reductions$harmony@cell.embeddings[, 'harmony_2']
mat.gene <- sample.pbmc_chond@assays$RNA@data
df.pc.gene <- data.frame(t(as.matrix(mat.gene)))
df.pc.gene$Harmony2 <- Harmony2
df.pc.gene$celltype <- sample.pbmc_chond$celltype
df.pc.gene$status <- sample.pbmc_chond$type

genes <- c('ING1', 'CCNL1', 'VIT', 'SCRG1', 'ELN', 'CYTL1', 'SPARC', 'COL2A1', 'COL9A3',
           'MT1X','ICAM1','IL6','MMP1','CXCL3','CXCL2')
list.plot <- list()
for (gene in genes) {
  df.plot <- df.pc.gene[, c('Harmony2', 'celltype', gene)]
  names(df.plot) <- c('Harmony2', 'celltype', 'gene')
  p.gene <- 
    ggplot(data = df.plot, aes(x = Harmony2, y = gene)) + 
    geom_point(aes(color = celltype), size = 0.1) + 
    scale_color_manual(labels = c('CSC', 'C0', 'C1', 'C2'),values = colors) + 
    xlim(-30, 10) + 
    geom_smooth(color = '#696969',size = 0.5) + 
    labs(x = '', y = '') + theme_bw()+
    theme(panel.background=element_rect(color='black',size =0.4,linetype ="solid"),
          axis.text = element_blank(),panel.grid = element_blank(),
          axis.ticks = element_blank(),plot.margin = unit(c(0.1,0.1,-0.5,-0.5), "cm"),
          legend.position = 'none') +
    annotate('text', label = gene, x = 10, y = max(df.plot$gene), 
             hjust = 1, vjust = 1, size = 3)
  list.plot[[gene]] <- p.gene
}

library(gridExtra)
p.merge <- marrangeGrob(list.plot, ncol = 5, nrow = 3, top = NULL)
ggsave(plot = p.merge, path = 'JingMA_NEW/res/Harmony/ALL/FIG/', 
       filename = 'sFig3C_chon_genes.pdf',
       height = 7, width = 18, units = 'cm')







