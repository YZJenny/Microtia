### 用不同形状画点dotplot

library(ggplot2)
LR_Control<- read.table('/home/yzj/JingMA_NEW/res/CellPhoneDB/out_Control/LR_filter.txt',header = T,sep='\t')
LR_Microtia<- read.table('/home/yzj/JingMA_NEW/res/CellPhoneDB/out_Microtia/LR_filter.txt',header = T,sep='\t')

print(dim(LR_Control))
print(dim(LR_Microtia))

LR_C <- paste(LR_Control$pair,LR_Control$clusters,sep=':')
LR_M <- paste(LR_Microtia$pair,LR_Microtia$clusters,sep=':')

OL <- intersect(LR_C,LR_M)
uniqC <- setdiff(LR_C,LR_M)
uniqM <- setdiff(LR_M,LR_C)

length(OL)
length(uniqC)
length(uniqM)

LR_Control$Occur <- 'Common'
LR_Control$Occur[LR_C %in% uniqC] <- 'Unique'

LR_Microtia$Occur <- 'Common'
LR_Microtia$Occur[LR_M %in% uniqM] <- 'Unique'


Type='Microtia'
if(Type=='Control'){
  df <- LR_Control
}else if(Type=='Microtia'){
  df <- LR_Microtia
}
filename=paste('/home/yzj/JingMA_NEW/res/CellPhoneDB/out_',Type,'/',Type,'_dotplot.pdf',sep='')

my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
ggplot(df,aes(x=clusters,y=pair,shape=Occur)) +
  geom_point(aes(size=-log10(pvalue),color=log2mean)) +
  scale_color_gradientn('Log2 mean (Gene 1, Gene 2)', colors=my_palette) +
  theme_bw() +
  theme(
    axis.text=element_text(size=14, colour = "black",family = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, family = 'bold'),
    axis.text.y = element_text(size=12, colour = "black", family = 'bold'),
    axis.title=element_blank(),
    text = element_text(family='bold'),
    panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
ggsave(filename, width = 30, height = 75, device = cairo_pdf, limitsize=F)

