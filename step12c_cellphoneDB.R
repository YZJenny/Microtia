library(ggplot2)
LR_Control<- read.table('/local/yzj/JingMA_NEW/res/CellPhoneDB_alone_subcelltype/out_Normal/LR_filter.txt',header = T,sep='\t')
LR_Microtia<- read.table('/local/yzj/JingMA_NEW/res/CellPhoneDB_alone_subcelltype/out_Microtia/LR_filter.txt',header = T,sep='\t')

subLR_Control <- LR_Control[LR_Control$clusters=]