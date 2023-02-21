
.libPaths("/lustre/tianlab/yzj/software/R-4.0.0/lib64/R/library")
library(Seurat)
library(SCENIC)
library(foreach)

######################
## 只看主要细胞类型
#######################
pbmc <- readRDS('/lustre/tianlab/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
Idents(pbmc) <- pbmc$celltype
subpbmc <- subset(pbmc,idents = c('CSC','C','SSC','SC'))
VG <- VariableFeatures(subpbmc)
exprMat <- as.matrix(subpbmc@assays$RNA@counts)
setwd('/lustre/tianlab/yzj/JingMA_NEW/res/SCENIC_main')
######################


### Initialize settings 初始设置，导入评分数据库
data(defaultDbNames)
dbs <- defaultDbNames[['hgnc']]
scenicOptions <- initializeScenic(org="hgnc", dbDir="/lustre/tianlab/yzj/publicData/SCENIC",dbs =dbs, nCores=60)
saveRDS(scenicOptions, file='int/scenicOptions.Rds')

### 共表达网络
print('Step2')
scenicOptions <- readRDS('int/scenicOptions.Rds')
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept,]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)
save.image(file='int/step1_1.Rdata')

### Build and score the GRN 构建并给基因调控网络（GRN）打分
print('Step3')
exprMat_log <- log2(exprMat+1)
#scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions) ##不指定method
save.image(file='int/step1_2.Rdata')
print('DONE!')
##

.libPaths(c("/home/yzj/R/x86_64-pc-linux-gnu-library/4.0","/home/zy/tools/R-4.0.0/library" ))
library(Seurat)
library(SCENIC)
library(foreach)

##
print('Step4')
setwd('/home/yzj/JingMA_NEW/res/SCENIC_main')
load('int/step1_2.Rdata')

scenicOptions@settings$nCores <- 1
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
scenicOptions@settings$nCores <-25
runSCENIC_4_aucell_binarize(scenicOptions)
save.image(file='int/step1_4.Rdata')
print('DONE!')