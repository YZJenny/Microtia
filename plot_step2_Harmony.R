library(Seurat)
library(dplyr)
library(ggplot2)

pbmc_batch <- readRDS('/home/yzj/JingMA_NEW/res/QC/ALL/RDS/PBMC_ORI.RDS')

## 1.ElbowPlot/JackStrawPlot
pdf('/home/yzj/JingMA_NEW/res/QC/ALL/FIG/sFig_ElbowPlot.pdf')
ElbowPlot(pbmc_batch,ndims = 150)
dev.off()

tmp <- JackStraw(pbmc_batch,dims=150, num.replicate = 100)
tmp <- ScoreJackStraw(tmp, dims = 1:150)
# JackStrawPlot(tmp, dims = 1:20)
# print(sum(pbmc_batch@reductions$pca@feature.loadings[1:100])/sum(pbmc_batch@reductions$pca@feature.loadings))
#saveRDS(tmp,'/home/yzj/JingMA_NEW/res/QC/ALL/FIG/sFig.RDS')


## 2. 没有去除batches之前的batch分布
pdf('/home/yzj/JingMA_NEW/res/QC/ALL/FIG/sFig_UMAP_ORI.pdf',width = 10,height = 10)
DimPlot(pbmc_batch, reduction ='umap', group.by='batch')+
  theme(axis.line = element_line(size=1, colour = "black"))
dev.off()


