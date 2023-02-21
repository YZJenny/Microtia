library("stringr")
library(tibble)
library(Seurat)
library(ggplot2)
source('/home/yzj/JingMA_ORI/res/CellPhoneDB/scr/cellphone_util.R')
source('/home/yzj/JingMA_ORI/res/CellPhoneDB/scr/cellphone_dotplot.R')
source('/home/yzj/JingMA_ORI/res/CellPhoneDB/scr/cellphone_heat.R')
source('/home/yzj/JingMA_ORI/res/CellPhoneDB/scr/cellphone_cytoscape.R')
source('/home/yzj/JingMA_ORI/res/CellPhoneDB/scr/cellphone_circle.R')

pbmc <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')

########
## 1. 正常和疾病分开跑，看celltype
########
Type='Microtia' # Normal
cells <- colnames(pbmc)[pbmc$type==Type]
subpbmc <- subset(pbmc,cells = cells)

## 1. prepare data
data.input  <- as.data.frame(subpbmc@assays$RNA@data)
data.input <- rownames_to_column(data.input,var = "Gene")

meta = data.frame(Cell = colnames(subpbmc),CellType =subpbmc$celltype)
write.table(meta,paste('JingMA_NEW/res/CellPhoneDB_alone_celltype/in/',Type,'_meta.txt',sep=''),row.names = FALSE,quote=FALSE,sep='\t')
write.table(data.input,paste('JingMA_NEW/res/CellPhoneDB_alone_celltype/in/',Type,'_count.txt',sep=''),row.names = FALSE,quote=FALSE,sep='\t')


## 2. run CellPhoneDB
# source /home/yjingjing/cpdb-venv/bin/activate
# cd /home/yzj/JingMA_NEW/res/CellPhoneDB_alone_celltype/
# nohup sh run.sh &


## 3. visival
Type='Microtia'
outDir=paste('/home/yzj/JingMA_NEW/res/CellPhoneDB_alone_celltype/out_',Type,'/',sep='')
setwd(outDir)

## Plot Heatmap
dot_plot(selected_rows = NULL,selected_columns = NULL,filename = paste(Type,'_dotplot1.pdf',sep=''),pvalue = 0.05,mean = 0,
         width = 30,height = 75,means_path = './means.txt',pvalues_path = './pvalues.txt', output_extension = '.pdf')

heatmaps_plot(meta_file = paste('../in/',Type,'_meta.txt',sep=''),pvalues_file = 'pvalues.txt',count_filename = paste(Type,'_heatmap_count.pdf',sep=''),log_filename = paste(Type,'_heatmap_log.pdf',sep=''),
              show_rownames = T,show_colnames = T,scale="none",cluster_cols = F,border_color='white',cluster_rows = F,fontsize_row=20,
              fontsize_col = 20, main = '',treeheight_row=0,family='Arial',treeheight_col = 0,col1 = "dodgerblue4",col2 = 'peachpuff',col3 = 'deeppink4',
              meta_sep='\t',pvalues_sep='\t',pvalue=0.05)

########
## 2. 正常和疾病分开跑，看cell subtype
########
Type='Microtia' # Normal
cells <- colnames(pbmc)[pbmc$type==Type]
subpbmc <- subset(pbmc,cells = cells)

## 1. prepare data
data.input  <- as.data.frame(subpbmc@assays$RNA@data)
data.input <- rownames_to_column(data.input,var = "Gene")

meta = data.frame(Cell = colnames(subpbmc),CellType =subpbmc$subcelltype)
write.table(meta,paste('JingMA_NEW/res/CellPhoneDB_alone_subcelltype/in/',Type,'_meta.txt',sep=''),row.names = FALSE,quote=FALSE,sep='\t')
write.table(data.input,paste('JingMA_NEW/res/CellPhoneDB_alone_subcelltype/in/',Type,'_count.txt',sep=''),row.names = FALSE,quote=FALSE,sep='\t')

## 2. run CellPhoneDB
# source /home/yjingjing/cpdb-venv/bin/activate
# cd /home/yzj/JingMA_NEW/res/CellPhoneDB_alone_subcelltype/
# nohup sh run.sh &


## 3. visival
Type='Microtia'
outDir=paste('/home/yzj/JingMA_NEW/res/CellPhoneDB_alone_subcelltype/out_',Type,'/',sep='')
setwd(outDir)

## Plot Heatmap
dot_plot(selected_rows = NULL,selected_columns = NULL,filename = paste(Type,'_dotplot1.pdf',sep=''),pvalue = 0.05,mean = 0,
         width = 30,height = 75,means_path = './means.txt',pvalues_path = './pvalues.txt', output_extension = '.pdf')

heatmaps_plot(meta_file = paste('../in/',Type,'_meta.txt',sep=''),pvalues_file = 'pvalues.txt',count_filename = paste(Type,'_heatmap_count.pdf',sep=''),log_filename = paste(Type,'_heatmap_log.pdf',sep=''),
              show_rownames = T,show_colnames = T,scale="none",cluster_cols = F,border_color='white',cluster_rows = F,fontsize_row=20,
              fontsize_col = 20, main = '',treeheight_row=0,family='Arial',treeheight_col = 0,col1 = "dodgerblue4",col2 = 'peachpuff',col3 = 'deeppink4',
              meta_sep='\t',pvalues_sep='\t',pvalue=0.05)


########
## 3. 正常和疾病和在一起跑，看celltype
########

## 1. prepare data
data.input  <- as.data.frame(pbmc@assays$RNA@data)
data.input <- rownames_to_column(data.input,var = "Gene")

meta = data.frame(Cell = colnames(pbmc),CellType =pbmc$celltype)
write.table(meta,'JingMA_NEW/res/CellPhoneDB_merge_celltype/in/meta.txt',row.names = FALSE,quote=FALSE,sep='\t')
write.table(data.input,'JingMA_NEW/res/CellPhoneDB_merge_celltype/in/count.txt',row.names = FALSE,quote=FALSE,sep='\t')

## 2. run CellPhoneDB
# source /home/yjingjing/cpdb-venv/bin/activate
# cd /home/yzj/JingMA_NEW/res/CellPhoneDB_merge_celltype/
# nohup sh run.sh &


## 3. visival

outDir='/home/yzj/JingMA_NEW/res/CellPhoneDB_merge_celltype/out/'
setwd(outDir)

## Plot Heatmap
dot_plot(selected_rows = NULL,selected_columns = NULL,filename = 'dotplot1.pdf',pvalue = 0.05,mean = 0,
         width = 30,height = 75,means_path = './means.txt',pvalues_path = './pvalues.txt', output_extension = '.pdf')

heatmaps_plot(meta_file = '../in/meta.txt',pvalues_file = 'pvalues.txt',count_filename = 'heatmap_count.pdf',log_filename = 'heatmap_log.pdf',
              show_rownames = T,show_colnames = T,scale="none",cluster_cols = F,border_color='white',cluster_rows = F,fontsize_row=20,
              fontsize_col = 20, main = '',treeheight_row=0,family='Arial',treeheight_col = 0,col1 = "dodgerblue4",col2 = 'peachpuff',col3 = 'deeppink4',
              meta_sep='\t',pvalues_sep='\t',pvalue=0.05)


########
## 4. 正常和疾病和在一起跑，看subcelltype
########

## 1. prepare data
data.input  <- as.data.frame(pbmc@assays$RNA@data)
data.input <- rownames_to_column(data.input,var = "Gene")

meta = data.frame(Cell = colnames(pbmc),CellType =pbmc$subcelltype)
write.table(meta,'JingMA_NEW/res/CellPhoneDB_merge_subcelltype/in/meta.txt',row.names = FALSE,quote=FALSE,sep='\t')
write.table(data.input,'JingMA_NEW/res/CellPhoneDB_merge_subcelltype/in/count.txt',row.names = FALSE,quote=FALSE,sep='\t')

## 2. run CellPhoneDB
# source /home/yjingjing/cpdb-venv/bin/activate
# cd /home/yzj/JingMA_NEW/res/CellPhoneDB_merge_subcelltype/
# nohup sh run.sh &


## 3. visival

outDir='/home/yzj/JingMA_NEW/res/CellPhoneDB_merge_subcelltype/out/'
setwd(outDir)

## Plot Heatmap
dot_plot(selected_rows = NULL,selected_columns = NULL,filename = 'dotplot1.pdf',pvalue = 0.05,mean = 0,
         width = 30,height = 75,means_path = './means.txt',pvalues_path = './pvalues.txt', output_extension = '.pdf')

heatmaps_plot(meta_file = '../in/meta.txt',pvalues_file = 'pvalues.txt',count_filename = 'heatmap_count.pdf',log_filename = 'heatmap_log.pdf',
              show_rownames = T,show_colnames = T,scale="none",cluster_cols = F,border_color='white',cluster_rows = F,fontsize_row=20,
              fontsize_col = 20, main = '',treeheight_row=0,family='Arial',treeheight_col = 0,col1 = "dodgerblue4",col2 = 'peachpuff',col3 = 'deeppink4',
              meta_sep='\t',pvalues_sep='\t',pvalue=0.05)


