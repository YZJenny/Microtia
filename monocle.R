.libPaths(c("/home/yzj/R/x86_64-pc-linux-gnu-library/4.0","/home/zy/tools/R-4.0.0/library"))
library(monocle)
library(Seurat)
library(scales)

# ############
# ## 1. for chond lineage
# ############
# mycol <- c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF")
# pbmc_chond <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Chond.Rdata')
# 
# ###################
# ## Extract data, phenotype data, and feature data from the SeuratObject
# ###################
# data <- as(as.matrix(pbmc_chond@assays$RNA@counts), 'sparseMatrix')
# pd <- new('AnnotatedDataFrame', data = pbmc_chond@meta.data)
# fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
# fd <- new('AnnotatedDataFrame', data = fData)
# 
# ###################
# ## Construct monocle cds
# ###################
# monocle <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,
#                           expressionFamily = negbinomial.size())
# head(pData(monocle))
# names(pData(monocle))[names(pData(monocle))=='celltype']="Cluster"
# 
# monocle <- estimateSizeFactors(monocle)
# monocle <- estimateDispersions(monocle)
# 
# 
# ###################
# ## choose genes that define a cell's progress
# ###################
# # #try1: Seurat MK Gene
# MK.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Chond.RDS') 
# MK.CSC <- MK.lst$CSC
# MK.TC <- MK.lst$TC
# MK.C1 <- MK.lst$C1
# MK.C2 <- MK.lst$C2
# 
# pos_markers <- c(rownames(MK.CSC)[MK.CSC$avg_logFC > log(1.5) & MK.CSC$p_val_adj < 0.05 ],
#                  rownames(MK.TC)[MK.TC$avg_logFC > log(1.5) & MK.TC$p_val_adj < 0.05 ],
#                  rownames(MK.C1)[MK.C1$avg_logFC > log(1.5) & MK.C1$p_val_adj < 0.05 ],
#                  rownames(MK.C2)[MK.C2$avg_logFC > log(1.5) & MK.C2$p_val_adj < 0.05 ])
# 
# length(pos_markers)
# ordering_genes <- pos_markers
# monocle <- setOrderingFilter(monocle, ordering_genes)
# 
# pdf('/home/yzj/JingMA_NEW/res/Monocle/Chond/ordering_genes.pdf')
# plot_ordering_genes(monocle)
# dev.off()
# 
# ###################
# ## reduce data dimensionality
# ###################
# monocle <- reduceDimension(monocle, max_components = 2, reduction_method = "DDRTree")
# 
# ###################
# ## order cells along the trajectory
# ###################
# monocle <- orderCells(monocle)
# saveRDS(monocle,'/home/yzj/JingMA_NEW/res/Monocle/Chond/monocle.RDS')
# 
# pdf('/home/yzj/JingMA_NEW/res/Monocle/Chond/monocle.pdf')
# plot_cell_trajectory(monocle, color_by = "Cluster",state_number_size=0.2,cell_size=0.8,cell_link_size=2)+
#   theme(axis.line = element_line(size=15, colour = "black"),
#         axis.text = element_text(size = 15))
# 
# plot_cell_trajectory(monocle, color_by = "State",state_number_size=2,cell_size=0.8,cell_link_size=2)+
#   theme(axis.line = element_line(size=15, colour = "black"),
#         axis.text = element_text(size = 15))
# 
# plot_cell_trajectory(monocle, color_by = "Pseudotime",state_number_size=2,cell_size=0.8,cell_link_size=2)+
#   theme(axis.line = element_line(size=15, colour = "black"),
#         axis.text = element_text(size = 15))
# plot_cell_trajectory(monocle, color_by = "Cluster",cell_size=0.8) +
#   facet_wrap(~Cluster, nrow = 1)
# dev.off()
# 
# 
# ###################
# ## root_state: order cells along the trajectory
# ###################
# monocle <- readRDS('/home/yzj/JingMA_NEW/res/Monocle/Chond/monocle.RDS')
# root_state= 'CSC'
# monocle_root <- orderCells(monocle, root_state = root_state)
# saveRDS(monocle_root,'/home/yzj/JingMA_NEW/res/Monocle/Chond/monocle_root.RDS')
# 
# pdf('/home/yzj/JingMA_NEW/res/Monocle/Chond/monocle_root.pdf')
# plot_cell_trajectory(monocle_root, color_by = "Pseudotime",state_number_size=2,cell_size=1,cell_link_size=2)+
#   theme(axis.line = element_line(size=15, colour = "black"),
#         axis.text = element_text(size = 15))
# dev.off()



############
## 2. for stroma lineage
############
pbmc_stroma <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Stroma.Rdata')

data <- as(as.matrix(pbmc_stroma@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = pbmc_stroma@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
head(pData(monocle))
names(pData(monocle))[names(pData(monocle))=='celltype']="Cluster"

monocle <- estimateSizeFactors(monocle)
monocle <- estimateDispersions(monocle)


MK.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Stroma.RDS') 
MK.SSC <- MK.lst$SSC
MK.SC1 <- MK.lst$SC1
MK.SC2 <- MK.lst$SC2

pos_markers <- c(rownames(MK.SSC)[MK.SSC$avg_logFC > log(1.5) & MK.SSC$p_val_adj < 0.05 ],
                 rownames(MK.SC1)[MK.SC1$avg_logFC > log(1.5) & MK.SC1$p_val_adj < 0.05 ],
                 rownames(MK.SC2)[MK.SC2$avg_logFC > log(1.5) & MK.SC2$p_val_adj < 0.05 ])

length(pos_markers)
ordering_genes <- pos_markers
monocle <- setOrderingFilter(monocle, ordering_genes)

pdf('/home/yzj/JingMA_NEW/res/Monocle/Stroma/ordering_genes.pdf')
plot_ordering_genes(monocle)
dev.off()

monocle <- reduceDimension(monocle, max_components = 2, reduction_method = "DDRTree")

###################
## order cells along the trajectory
###################
monocle <- orderCells(monocle)
saveRDS(monocle,'/home/yzj/JingMA_NEW/res/Monocle/Stroma/monocle.RDS')

pdf('/home/yzj/JingMA_NEW/res/Monocle/Stroma/monocle.pdf')
plot_cell_trajectory(monocle, color_by = "Cluster",state_number_size=0.2,cell_size=0.8,cell_link_size=2)+
  theme(axis.line = element_line(size=15, colour = "black"),
        axis.text = element_text(size = 15))

plot_cell_trajectory(monocle, color_by = "State",state_number_size=2,cell_size=0.8,cell_link_size=2)+
  theme(axis.line = element_line(size=15, colour = "black"),
        axis.text = element_text(size = 15))

plot_cell_trajectory(monocle, color_by = "Pseudotime",state_number_size=2,cell_size=0.8,cell_link_size=2)+
  theme(axis.line = element_line(size=15, colour = "black"),
        axis.text = element_text(size = 15))
plot_cell_trajectory(monocle, color_by = "Cluster",cell_size=0.8) +
  facet_wrap(~Cluster, nrow = 1)
dev.off()


###################
## root_state: order cells along the trajectory
###################
monocle <- readRDS('/home/yzj/JingMA_NEW/res/Monocle/Stroma/monocle.RDS')
root_state= '1'
monocle_root <- orderCells(monocle, root_state = root_state)
saveRDS(monocle_root,'/home/yzj/JingMA_NEW/res/Monocle/Stroma/monocle_root.RDS')

pdf('/home/yzj/JingMA_NEW/res/Monocle/Stroma/monocle_root.pdf')
plot_cell_trajectory(monocle_root, color_by = "Pseudotime",state_number_size=2,cell_size=1,cell_link_size=2)+
  theme(axis.line = element_line(size=15, colour = "black"),
        axis.text = element_text(size = 15))
dev.off()

