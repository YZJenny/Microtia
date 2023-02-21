.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(ggplot2)
library(reshape2)

path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
path.lineage <- paste0(path.data, 'chon_lineage/')
file.chon <- paste0(path.lineage, 'seurat_celltype.Rdata')
seurat.chon <- readRDS(file.chon)
# M1 M2 M3
seurat.child <- subset(seurat.chon, subset = batch %in% c('C4', 'C6', 'M1', 'M2', 'M3'))

####################
## Fig4C GO term
####################
fc.cutoff <- 0.4
path.M123 <- '/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/'
path.cutoff <- paste0(path.M123, 'cutoff_', fc.cutoff, '/')
file.marker.go <- paste0(path.cutoff, 'marker_go.Rdata')
list.marker.go <- readRDS(file.marker.go)

file.go.BP <- paste0(path.cutoff, 'GO_BP.Rdata')
list.go.BP <- readRDS(file.go.BP)
file.go.MF <- paste0(path.cutoff, 'GO_MF.Rdata')
list.go.MF <- readRDS(file.go.MF)
# file.go.CC <- paste0(path.cutoff, 'GO_CC.Rdata')
# list.go.CC <- readRDS(file.go.CC)
# file.kegg <- paste0(path.cutoff, 'kegg.Rdata')
# list.kegg <- readRDS(file.kegg)

file.go.BP <- paste0(path.cutoff, 'GO_BP_all.Rdata')
list.go.BP <- readRDS(file.go.BP)
file.go.MF <- paste0(path.cutoff, 'GO_MF_all.Rdata')
list.go.MF <- readRDS(file.go.MF)
file.go.CC <- paste0(path.cutoff, 'GO_CC_all.Rdata')
list.go.CC <- readRDS(file.go.CC)

# select GO
# df.GO <- data.frame(stringsAsFactors = F)
# Chondral stem cell
# GO.BP.CSC.M <- c('response to oxidative stress', 'response to unfolded protein', 
#                  'response to tumor necrosis factor', 'RNA splicing',
#                  'RNA localization', 'positive regulation of defense response')
# GO.BP.CSC.N <- c('ribosome biogenesis', 'response to copper ion', 'oxidative phosphorylation',
#                  'extracellular matrix organization', 'skeletal system development',
#                  'cell aggregation', 'cellular zinc ion homeostasis')
sel.GO.BP.N <- c('protein localization to endoplasmic reticulum', 
                 'translational initiation', 
                 # 'ribosome biogenesis', 
                 'oxidative phosphorylation', 'electron transport chain', 
                 'cellular respiration',
                 # 'mitochondrion organization', 'cellular oxidant detoxification', 
                 'cartilage development', 'chondrocyte differentiation',
                 # 'cartilage condensation',
                 'cellular response to zinc ion', 'cellular response to copper ion',
                 'cell aggregation', 'extracellular matrix organization')
sel.GO.BP.M <- c('response to oxidative stress', 
                 'cell death in response to oxidative stress',
                 # 'nitric oxide biosynthetic process',
                 'I-kappaB kinase/NF-kappaB signaling',
                 'p38MAPK cascade', 
                 # 'ERK1 and ERK2 cascade',
                 'intrinsic apoptotic signaling pathway', 
                 'extrinsic apoptotic signaling pathway',
                 'response to unfolded protein', 
                 'regulation of RNA stability',
                 # 'positive regulation of defense response', 
                 'activation of innate immune response',
                 'cellular response to tumor necrosis factor',
                 # 'cellular response to interferon-gamma',
                 'cellular response to interleukin-1',
                 'regulation of inflammatory response',
                 'cell cycle arrest',
                 'negative regulation of stem cell differentiation',
                 'negative regulation of cell growth',
                 # 'RNA splicing', 'RNA localization', 
                 'angiogenesis', 
                 # 'vascular endothelial growth factor production', 
                 # 'positive regulation of vasculature development',
                 'positive regulation of cell migration',
                 'negative regulation of cell adhesion', 
                 'extracellular matrix organization')
# sel.GO.MF.N <- c('oxidoreductase activity, acting on a heme group of donors',
#                  # 'structural constituent of ribosome', 
#                  'extracellular matrix structural constituent', 
#                'extracellular matrix binding', 'S100 protein binding')
sel.GO.MF.N <- c('extracellular matrix structural constituent', 
                 'structural constituent of ribosome',
                 'S100 protein binding')
sel.GO.MF.M <- c('extracellular matrix structural constituent')

sort.GO <- c('S100 protein binding',
             'cellular response to zinc ion', 'cellular response to copper ion',
             'cell aggregation', 'extracellular matrix organization',
             'extracellular matrix structural constituent', 
             'cartilage development', 'chondrocyte differentiation',
             'structural constituent of ribosome',
             'protein localization to endoplasmic reticulum', 
             'translational initiation', 
             'cellular respiration',
             'oxidative phosphorylation', 'electron transport chain', 
             'response to oxidative stress',
             'cell death in response to oxidative stress',
             'I-kappaB kinase/NF-kappaB signaling',
             'p38MAPK cascade', 
             'intrinsic apoptotic signaling pathway', 
             'extrinsic apoptotic signaling pathway',
             'response to unfolded protein', 
             'regulation of RNA stability',
             'activation of innate immune response',
             'regulation of inflammatory response',
             'cellular response to tumor necrosis factor',
             'cellular response to interleukin-1',
             'negative regulation of stem cell differentiation',
             'negative regulation of cell growth',
             'cell cycle arrest',
             'angiogenesis', 
             'positive regulation of cell migration',
             'negative regulation of cell adhesion')

sort.GO <- c('extracellular matrix structural constituent', 
             'electron transport chain', 
             'cellular response to zinc ion', 'cellular response to copper ion',
             'cell aggregation', 'S100 protein binding',
             'extracellular matrix organization',
             'oxidative phosphorylation', 
             'structural constituent of ribosome',
             'protein localization to endoplasmic reticulum', 
             'translational initiation', 
             'cellular respiration',
             'cartilage development', 'chondrocyte differentiation',
             'regulation of inflammatory response',
             'cellular response to tumor necrosis factor',
             'cellular response to interleukin-1',
             'angiogenesis', 
             'positive regulation of cell migration',
             'response to oxidative stress',
             'I-kappaB kinase/NF-kappaB signaling',
             'p38MAPK cascade', 
             'intrinsic apoptotic signaling pathway', 
             'extrinsic apoptotic signaling pathway',
             'response to unfolded protein', 
             'regulation of RNA stability',
             'activation of innate immune response',
             'negative regulation of cell adhesion',
             'cell death in response to oxidative stress',
             'negative regulation of stem cell differentiation',
             'negative regulation of cell growth',
             'cell cycle arrest')

terms <- c("CSC_Microtia_increase",
           "CSC_Microtia_decrease",
           "TC_Microtia_increase",
           "TC_Microtia_decrease",
           "C1_Microtia_increase",
           "C1_Microtia_decrease",
           "C2_Microtia_increase",
           "C2_Microtia_decrease")

df.plot <- data.frame(stringsAsFactors = F)
for (term in terms) {
    cell <- strsplit(term, split = '_')[[1]][1]
    status <- strsplit(term, split = '_')[[1]][3]
    sub.BP <- list.go.BP[[term]]
    rownames(sub.BP) <- sub.BP$Description
    sub.BP <- sub.BP[sub.BP$pvalue < 0.018,]
    if (status == 'increase') {
        sel.BP <- sub.BP[sel.GO.BP.M, c('Description', 'pvalue', 'geneID')]
        sel.BP$Description <- sel.GO.BP.M
    } else {
        sel.BP <- sub.BP[sel.GO.BP.N, c('Description', 'pvalue', 'geneID')]
        sel.BP$Description <- sel.GO.BP.N
    }
    sub.MF <- list.go.MF[[term]]
    rownames(sub.MF) <- sub.MF$Description
    sub.MF <- sub.MF[sub.MF$pvalue < 0.018,]
    if (status == 'increase') {
        sel.MF <- sub.MF[sel.GO.MF.M, c('Description', 'pvalue', 'geneID')]
        sel.MF$Description <- sel.GO.MF.M
    } else {
        sel.MF <- sub.MF[sel.GO.MF.N, c('Description', 'pvalue', 'geneID')]
        sel.MF$Description <- sel.GO.MF.N
    }
    sub.plot <- rbind(sel.BP, sel.MF)
    sub.plot$pvalue[is.na(sub.plot$pvalue)] <- 1
    sub.plot$CellType <- rep(cell, nrow(sub.plot))
    if (status == 'increase') {
        sub.plot$Status <- rep('Microtia', nrow(sub.plot))
        sub.plot$Coeff <- rep(1, nrow(sub.plot))
    } else {
        sub.plot$Status <- rep('Normal', nrow(sub.plot))
        sub.plot$Coeff <- rep(-1, nrow(sub.plot))
    }
    df.plot <- rbind(df.plot, sub.plot)
}
df.plot$log10Pval <- -log10(df.plot$pvalue)
df.plot$log10Pval[abs(df.plot$log10Pval) > 10] = 10
df.plot$log10Pval <- df.plot$log10Pval * df.plot$Coeff
df.plot$abs_log10Pval <- abs(df.plot$log10Pval)
# df.plot$Description <- factor(df.plot$Description, 
#                               levels = unique(c(rev(sel.GO.MF.N), 
#                                                 rev(sel.GO.BP.N), 
#                                                 rev(sel.GO.BP.M))))
df.plot$Description <- factor(df.plot$Description, levels = rev(sort.GO))
df.plot$CellType <- factor(df.plot$CellType, 
                           levels = c('CSC', 'TC', 'C1', 'C2'))

df.plot$col_name <- paste(df.plot$CellType, df.plot$Status, sep = '_')


mat.plot <- reshape2::dcast(df.plot, Description ~ col_name, value.var = 'log10Pval')
row.names(mat.plot) <- mat.plot$Description
mat.plot$Description <- NULL
mat.plot[is.na(mat.plot)] <- 0

# col annotation
annotation_col = data.frame(
    CellType = factor(c(rep('CSC', 2), rep('C1', 2), rep('C2', 2), rep('C0', 2)), 
                      levels = c('CSC', 'C0', 'C1', 'C2')), 
    Status = factor(rep(c('Microtia', 'Normal'), 4), levels = c('Normal', 'Microtia')),
    row.names = colnames(mat.plot)
)

cols <- c("CSC_Microtia", "TC_Microtia",
          "C1_Microtia", "C2_Microtia",
          "CSC_Normal", "TC_Normal",
          "C1_Normal", "C2_Normal")
mat.plot <- mat.plot[rev(rownames(mat.plot)), cols]
annotation_col <- annotation_col[cols,]

# ann_colors = list(
#     CellType = rep(c("#33A02C","#B2DF8A","#1F78B4","#A6CEE3"), 2),
#     Status = c(Microtia = rep("#F8766D",4), Normal = rep("#00BFC4",4)))
ann_colors = list(
    CellType = c(CSC="#EE9572",C0="#B2DF8A",C1="#A6CEE3",C2="#9999FF"),
    Status = c(Microtia = "#637FBF", Normal = "#6C6C6C")
)
cols.sort <- c("CSC_Normal", "TC_Normal",
          "C1_Normal", "C2_Normal",
          "CSC_Microtia", "TC_Microtia", 
          "C1_Microtia", "C2_Microtia")
mat.plot <- mat.plot[, cols.sort]
mat.plot <- apply(mat.plot, 2, as.numeric)
rownames(mat.plot) <- sort.GO

bk <- c(seq(-10,-0.1,by=0.01),seq(0,10,by=0.01))
# labels_row <- gsub('oxidoreductase activity, acting on a heme group of donors',
#                    'oxidoreductase activity', rownames(mat.plot))


bk <- c(seq(-10,-0.1,by=0.01),seq(0,10,by=0.01))
plot.heatmap <- pheatmap::pheatmap(mat.plot,fontsize = 6,fontsize_row = 6,fontsize_col = 6,
                                   cluster_rows = F, cluster_cols = F, scale = "none",
                                   display_numbers = F,border_color='grey',
                                   annotation_col = annotation_col ,annotation_colors = ann_colors,
                                   show_rownames = T, show_colnames = F, legend = TRUE, 
                                   gaps_col = c(4),lwd=0.3,
                                   color = c(colorRampPalette(colors = c("red","white"))(length(bk)/2),colorRampPalette(colors = c("white","blue"))(length(bk)/2)),
                                   legend_breaks=seq(-10,10,5),
                                   breaks=bk)

ggsave('/home/yzj/JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/Fig4C_heatmap_GO_enrich_sort.pdf',plot.GO,
       height = 10, width = 13, units = 'cm')


sel.genes <- c('PRDX2', 'SOD3', 'NDUFA13', 'COX7C',
               # 'PRDX2', 'PRDX4', 'SOD3', # 抗氧化
               # 'NDUFA1', 'NDUFA4', 'COX7C', 'COX5B', 
               # 'NDUFA11', 'NDUFA13', 'COX6C', 'COX8A', # 线粒体/呼吸链
               'SOD2', 'KDM6B', 'NFE2L2',  # 氧化应激
               'NFKB1', 'NFKB2',  # NF kappa
               'MAP2K3', 'MAPKAPK2', # MAPK 
               'HNRNPK', 'MCL1', 'PMAIP1', 'PPP1R15A',  # 凋亡
               'RPS28', 'RPS10', # 核糖体生成与翻译起始
               'TCP1', 'HSPA9', 'HNRNPC',  # 维持蛋白质/RNA稳定
               'TNFAIP3', 'CD44', 'HLA-C', 'IFNGR2', # 免疫
               'KMT2E', 'PPP2CA', 'CDKN1A', 'TSPYL2', # 抑制干细胞
               'PDPN', 'DLC1', 'CSF1', 'CLEC3A', # 细胞迁移与粘连
               'TNFAIP2', 'ESM1', 'ZC3H12A', # 血管生成
               'MMP3', # 基质解组装
               'LAMB3', 'FN1', # extracellular matrix
               'COL2A1', 'COL9A3', 'COL11A1', 'SPARC', 'MGP', 'ELN', 'ACAN',
               'FRZB', 'CTGF', 'VIT',  # 软骨发育
               'S100A1', 'S100B', # S100 protein
               'MT1X', 'MT1E') # 离子
sel.genes <- rev(sel.genes)
cells <- c('CSC', 'TC', 'C1', 'C2')
file.gsea.marker <- paste0(path.M123, 'marker_GSEA.Rdata')
list.marker.gsea <- readRDS(file.gsea.marker)

df.plot.gene <- data.frame(stringsAsFactors = F)
for (cell in cells) {
    sub.marker <- list.marker.gsea[[cell]][sel.genes, c('p_val', 'avg_logFC')]
    sub.seurat <- subset(seurat.child, subset = celltype == cell)
    sub.mat <- as.matrix(sub.seurat@assays$RNA@counts[sel.genes,])
    gene.prop <- (rowSums(sub.mat != 0)/ncol(sub.mat))*100
    sub.marker$Prop <- gene.prop
    sub.marker$CellType <- rep(cell, length(sel.genes))
    sub.marker$Gene <- sel.genes
    df.plot.gene <- rbind(df.plot.gene, sub.marker)
}
df.plot.gene$avg_logFC[is.na(df.plot.gene$avg_logFC)] <- 0
df.plot.gene$avg_logFC[df.plot.gene$avg_logFC < -1.5] <- -1.5
df.plot.gene$avg_logFC[df.plot.gene$avg_logFC > 1.5] <- 1.5
df.plot.gene$p_val[is.na(df.plot.gene$p_val)] <- 1
df.plot.gene$log_pval <- -log10(df.plot.gene$p_val)
df.plot.gene$log_pval[df.plot.gene$log_pval > 200] <- 200
df.plot.gene$log_log_pval <- log(df.plot.gene$log_pval + 1)
df.plot.gene$Gene <- factor(df.plot.gene$Gene, levels = sel.genes)
df.plot.gene$CellType <- factor(df.plot.gene$CellType, levels = cells)

plot.bubble <- 
    ggplot(data = df.plot.gene, 
           aes(x = CellType, y = Gene, size = Prop, color = avg_logFC)) + 
    geom_point(fill = 'cornsilk') + 
    scale_size_continuous(range = c(0,1)) +
    scale_colour_gradient2(low = 'blue', mid = "white", high = 'red',
                           breaks = c(-1, 0, 1)) + 
    scale_y_discrete(position = 'right') + 
    scale_x_discrete(breaks = cells,
                     labels = c('CSC', 'C0', 'C1', 'C2')) +
    labs(x = '', y = '', 
         # size = expression(paste("-log"[10], "(", italic("P"), "-value)")), 
         size = 'Proportion',
         color = 'logFC') + 
    theme(panel.background = element_rect(color = 'black', size = 0.5,fill = 'transparent'), 
          axis.ticks.length =  unit(0.05, "cm"),
          axis.text.x = element_text(size = 7, color = "black"),
          axis.text.y = element_text(size = 5, color = "black", face = 'italic'),
          legend.key.size = unit(0.2,'cm'),legend.key.height = unit(0.2,'cm'),
          legend.text = element_text(size=5),legend.title  = element_text(size=5),legend.key.width = unit(0.2,'cm'),)

ggsave('/home/yzj/JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/Fig4D_marker_genes.pdf',plot.bubble,
       height = 10, width = 6, units = 'cm')


