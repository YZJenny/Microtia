##################
###补充材料： 1. 计算CV
##################
library(ggplot2)
library(Seurat)
fun_CV <- function(x){
  cv <- sd(x)/mean(x)*100
  return(cv)
}

pbmc_C <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
pbmc_C$Phase <- factor(pbmc_C$Phase,levels = c('Children','Adults'))
celltypes <- levels(pbmc_C$celltype)
print(celltypes)

get_CV <- function(CT){
  pbmc_Children <- subset(pbmc_C,cells = colnames(pbmc_C)[pbmc_C$Phase=='Children' & pbmc_C$celltype==CT])
  pbmc_Adults <- subset(pbmc_C,cells = colnames(pbmc_C)[pbmc_C$Phase=='Adults' & pbmc_C$celltype==CT])
  
  EXP_Children <- pbmc_Children@assays$RNA@scale.data
  EXP_Adults <- pbmc_Adults@assays$RNA@scale.data

  N=min(ncol(EXP_Children),ncol(EXP_Adults))
  cell_Children <- sample(colnames(EXP_Children),N)
  cell_Adults <- sample(colnames(EXP_Adults),N)
  
  subEXP_Children <- EXP_Children[,colnames(EXP_Children) %in% cell_Children]
  subEXP_Adults <- EXP_Adults[,colnames(EXP_Adults) %in% cell_Adults]
  DIST=abs(subEXP_Children-subEXP_Adults)
  cv_value <- (apply(DIST,1,fun_CV))
  return(cv_value)
}

CV_CSC <- get_CV('CSC')
CV_TC <- get_CV('C0')
CV_C1 <- get_CV('C1')
CV_C2 <- get_CV('C2')

plot.df <- data.frame(gene=names(CV_CSC),CSC=CV_CSC,TC=CV_TC,C1=CV_C1,C2=CV_C2)
plot.df <- reshape2::melt(plot.df)
colnames(plot.df) <- c('gene','CellType','CV')

library(ggplot2)
library(ggpubr)
library(ggsignif)
plot.df <- na.omit(plot.df)
p <- ggplot(plot.df, aes(x=CellType, y=CV,fill=CellType)) + 
  geom_violin(trim=FALSE,color= 'black') + 
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.size = 0)+
  scale_fill_manual(values =  c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF"))+ 
  theme_classic()+ 
  theme(axis.text=element_text(size=6,colour="black"),
        axis.title=element_text(size = 6,colour="black"), 
        axis.ticks.y =element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 6,face = 'bold'),
        plot.background = element_rect(fill = "#EFEFEF",colour = '#EFEFEF'),
        panel.background = element_rect(fill = '#EFEFEF'))+
  labs(title="", x="", y="")+
  guides(fill=FALSE)
p
ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/sFig5C_CV.pdf',p,width = 6,height = 6,units = 'cm')

p <- ggplot(plot.df, aes(x=CellType, y=CV,fill=CellType)) + 
  geom_violin(trim=FALSE,color= 'black') + 
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.size = 0)+
  scale_fill_manual(values =  c("#EE9572","#B2DF8A" ,"#A6CEE3","#9999FF"))+ 
  theme_classic()+ 
  theme(axis.text=element_text(size=15,colour="black"),
        axis.title=element_text(size = 15,colour="black"), 
        axis.ticks.y =element_line(colour="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 15,face = 'bold'))+
  labs(title="Zoom", x="Chondral subtypes", y="Coeffient of variation(CV)")+
  ylim(0, 800)+
  guides(fill=FALSE)
p
ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/sFig5C_CV_Zoom.pdf',p,width = 6,height = 6)



##################
###补充材料： 2. 计算DEG数目
##################
library(xlsx)
library(Seurat)
library(ggplot2)

batch1='Children'
batch2='Adults'
resol='celltype'

pbmc_C <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
pbmc_C$Phase <- factor(pbmc_C$Phase,levels = c('Children','Adults'))
VG <- VariableFeatures(pbmc_C)

DEG.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inChildrenAdults.RDS')
names(DEG.lst) <- c('CSPC','EC','IC','LC')
wp=paste('JingMA_NEW/res/compControl/ChildrenvsAdults/')

FC <- 1.5
adjP <- 0.05

CellTypes=c('CSPC','EC','IC','LC')
upNum <- c()
dnNum <- c()
upGene.lst <- list()
dnGene.lst <- list()
for(i in 1:length(CellTypes)){
  CT=CellTypes[i]
  print(CT)

  data<- DEG.lst[[CT]]
  data <- tibble::rownames_to_column(data,'gene')
  ## 相对于成人来说
  data$avg_logFC <- -data$avg_logFC
  sigUP_gene <- data[data$p_val_adj  < adjP & data$avg_logFC > log(FC),]
  sigDN_gene <- data[data$p_val_adj  < adjP & data$avg_logFC < -log(FC),]
  upGene.lst[[CT]] <- sigUP_gene
  dnGene.lst[[CT]] <- sigDN_gene
  
  upNum <- c(upNum,nrow(sigUP_gene))
  dnNum <- c(dnNum,nrow(sigDN_gene))
}

Num <- data.frame(CellType=CellTypes,upGene=upNum,dnGene=dnNum)
Num$CellType <- factor(Num$CellType,levels = CellTypes)
library(reshape2)
NUM_meltdf <- melt(Num)

p <- ggplot(data=NUM_meltdf, mapping=aes(x=CellType,y=value,fill=variable))+
  geom_bar(stat="identity",width=0.7,position='dodge')+theme_classic()+
  theme(axis.text = element_text(size = 7,colour = 'black'),axis.title = element_text(size = 8,colour = 'black'),
        legend.key.height = unit(0.5,'cm'),legend.key.width = unit(0.5,'cm'))+
  scale_fill_manual(values=c('#B15E72','#7F99CE'))+
  labs(x="Cell type", y="DEG numbers")

ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/sFig_DEGnum_CA.pdf',p,height = 6,width = 10,units = 'cm')


##################
###补充材料：3. Venn plot of DEG
##################
library(VennDiagram)
library(RColorBrewer);
mycolor <- brewer.pal(9,"Set1");

## Upregulated gene
p <- venn.diagram(list(CSPC=upGene.lst$CSPC$gene,LC=upGene.lst$LC$gene,EC=upGene.lst$EC$gene,IC=upGene.lst$IC$gene),
             resolution = 300,alpha=c(0.5,0.5,0.5,0.5),lwd=rep(2,4),lty=rep(1,4),cat.cex = 1,cex =1,
             fill=c("#EE9572","#9999FF","#B2DF8A" ,"#A6CEE3"), col=c("#EE9572","#9999FF","#B2DF8A" ,"#A6CEE3"), 
             cat.fontface=2,
             main="Upregulated DEGs",imagetype = "tiff",
             height = 3,width = 3,units = 'cm',
             filename = NULL)

pdf("JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/sFig_upDEGveen_CA.pdf",height = 2,width = 2)
grid.draw(p)
dev.off()

## Downregulated gene
p <- venn.diagram(list(CSPC=dnGene.lst$CSPC$gene,LC=dnGene.lst$LC$gene,EC=dnGene.lst$EC$gene,IC=dnGene.lst$IC$gene),
                  resolution = 300,alpha=c(0.5,0.5,0.5,0.5),lwd=rep(2,4),lty=rep(1,4),cat.cex = 1,cex =1,
                  fill=c("#EE9572","#9999FF","#B2DF8A" ,"#A6CEE3"), col=c("#EE9572","#9999FF","#B2DF8A" ,"#A6CEE3"), 
                  cat.fontface=2,
                  main="Downregulated DEGs",imagetype = "tiff",
                  height = 3,width = 3,units = 'cm',
                  filename = NULL)

pdf("JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/sFig_dnDEGveen_CA.pdf",height = 2,width = 2)
grid.draw(p)
dev.off()


##################
###补充材料：DEG里maker gene的ratio,不符合预期！
##################
MK.lst <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/Markers_celltype_Chond.RDS') #log(FC)>log(1.5), p_val_adj<0.05
names(MK.lst) <- c('CSPC','EC','IC','LC')

ratio <- c()
ratio.up <- c()
ratio.dn <- c()
for(i in 1:length(MK.lst)){
  CT=names(MK.lst)[i]
  MK <- MK.lst[[CT]]
  MK <- rownames(MK)[MK$avg_logFC > 0.25 & MK$p_val_adj < 0.05]
 
  dnDEG <- dnGene.lst[[CT]][,'gene']
  upDEG <- upGene.lst[[CT]][,'gene']
  DEG <- c(dnDEG,upDEG)
  
  DEG_MK <- intersect(MK,DEG)
  ratio <- c(length(DEG_MK)/length(MK),ratio)
  
  upDEG_MK <- intersect(MK,upDEG)
  ratio.up <- c(length(upDEG_MK)/length(MK),ratio.up)
  
  dnDEG_MK <- intersect(MK,dnDEG)
  ratio.dn <- c(length(dnDEG_MK)/length(MK),ratio.dn)
}
print(ratio)
print(ratio.up)
print(ratio.dn)


##################
###补充材料：4.不在GenAge and Aging Atlas里的代表性基因热图
##################
### gene associated aging
GA_mtx <- read.csv('publicData/GenAge/human_genes/genage_human.csv');Aging.db <- GA_mtx[,2]
library(xlsx)
Aging.BIG_df <- read.xlsx('publicData/GeneAtlas/aging_list_V2.0.xlsx',sheetIndex = 1)
Aging.BIG_df <- as.data.frame(Aging.BIG_df);Aging.BIG <- Aging.BIG_df$Symbol

Aging.genes <- union(Aging.db,Aging.BIG)
length(Aging.genes)

## Upregulated
dnValues_mtx <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_DNmtx.RDS')
dnValues_mtx <- as.data.frame(dnValues_mtx)
colnames(dnValues_mtx)[2] <- 'C0'

dnValues_mtx_noAging <-dnValues_mtx[!rownames(dnValues_mtx)%in%Aging.genes,] 
up_sum <- apply(dnValues_mtx_noAging,1,sum)
print(names(up_sum)[up_sum==4])


## Downregulated
upValues_mtx <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_UPmtx.RDS')
upValues_mtx <- as.data.frame(upValues_mtx)
colnames(upValues_mtx)[2] <- 'C0'

upValues_mtx_noAging <-upValues_mtx[!rownames(upValues_mtx)%in%Aging.genes,] 
dn_sum <- apply(upValues_mtx_noAging,1,sum)
print(names(dn_sum)[up_sum==4])

## 导入 GO term
GO <- read.gmt('/home/yzj/publicData/GMT/c5.all.v6.2.symbols.gmt')
#ECM_as <- GO[GO$term == 'GO_EXTRACELLULAR_MATRIX_ASSEMBLY','gene'] 
ECM_dis <- GO[GO$term=='GO_EXTRACELLULAR_MATRIX_DISASSEMBLY','gene']
ECM_comp <- GO[GO$term == 'GO_EXTRACELLULAR_MATRIX_COMPONENT','gene']
chond_diff <- GO[GO$term %in% c('GO_CHONDROCYTE_DIFFERENTIATION','GO_CHONDROCYTE_DEVELOPMENT'),'gene']

ECM_comp.gene <- intersect(ECM_comp,rownames(upValues_mtx_noAging))
chond_diff.gene <- intersect(chond_diff,rownames(upValues_mtx_noAging))

# pick.dnOL <- intersect(ECM_comp,rownames(upValues_mtx_noAging))
# pick.upOL <- intersect(ECM_dis,rownames(dnValues_mtx_noAging))

DEG.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inChildrenAdults.RDS')
names(DEG.lst) <- c('CSPC','EC','IC','LC')

pick.gene <- c(ECM_comp.gene,chond_diff.gene)
FC_df <- c()
for(i in 1:length(pick.gene)){
  gene=pick.gene[i]
  CSPC.FC <- DEG.lst$CSPC;CSPC.FC <- CSPC.FC$avg_logFC[rownames(CSPC.FC)==gene]
  EC.FC <- DEG.lst$EC;EC.FC <- EC.FC$avg_logFC[rownames(EC.FC)==gene]
  IC.FC <- DEG.lst$IC;IC.FC <- IC.FC$avg_logFC[rownames(IC.FC)==gene]
  LC.FC <- DEG.lst$LC;LC.FC <- LC.FC$avg_logFC[rownames(LC.FC)==gene]
  fc <- c(CSPC.FC,EC.FC,IC.FC,LC.FC)
  FC_df <- rbind(FC_df,fc)
}
colnames(FC_df) <- c('CSPC','EC','IC','LC'); rownames(FC_df) <- pick.gene
FC_df  <- -(FC_df)
FC_df[FC_df > -log(1.5,2)] <- 0

bk <- c(seq(-3,-0.1,by=0.01),seq(0,3,by=0.01))
plot.heatmap <- pheatmap::pheatmap(t(FC_df),fontsize = 6,fontsize_row = 6,fontsize_col = 6,
                                   cluster_rows = F, cluster_cols = F, scale = "none",
                                   display_numbers = F,border_color='grey',
                                   show_rownames = T, show_colnames = T, legend = TRUE, 
                                   gaps_col = c(length(ECM_comp.gene)),lwd=0.5,
                                   color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
                                   legend_breaks=seq(-3,3,1),
                                   breaks=bk)

ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/sFig_pickDEG_FC.pdf',plot.heatmap,width = 19 ,height =6 ,units = 'cm')


######
## 补充材料：5.关注共同的基因的富集分析
######
library(ggplot2)
upValues_mtx <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_UPmtx.RDS')
dnValues_mtx <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_DNmtx.RDS')

up_sum <- apply(upValues_mtx,1,sum)
commonUP <- names(up_sum)[which(up_sum > 1)]
up_df <- upValues_mtx[-(which(up_sum>1)),]
uniqUP_CSC <- rownames(up_df)[which(up_df[,1]==1)]
uniqUP_C0 <- rownames(up_df)[which(up_df[,2]==1)]
uniqUP_C1 <- rownames(up_df)[which(up_df[,3]==1)]
uniqUP_C2 <- rownames(up_df)[which(up_df[,4]==1)]
UPlst <- list(common=commonUP,CSC=uniqUP_CSC,C0=uniqUP_C0,C1=uniqUP_C1,C2=uniqUP_C2)

dn_sum <- apply(dnValues_mtx,1,sum)
commonDN <- names(dn_sum)[which(dn_sum > 1)]
dn_df <- dnValues_mtx[-(which(dn_sum>1)),]
uniqDN_CSC <- rownames(dn_df)[which(dn_df[,1]==1)]
uniqDN_C0 <- rownames(dn_df)[which(dn_df[,2]==1)]
uniqDN_C1 <- rownames(dn_df)[which(dn_df[,3]==1)]
uniqDN_C2 <- rownames(dn_df)[which(dn_df[,4]==1)]
DNlst <- list(common=commonDN,CSC=uniqDN_CSC,C0=uniqDN_C0,C1=uniqDN_C1,C2=uniqDN_C2)


### gene>1 富集分析
library(xlsx)
data(c2BroadSets)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(AnnotationHub)
library(org.Hs.eg.db)   #人类注释数据库
library(clusterProfiler)
library(DOSE)
library(dplyr)
library(tidyverse)
library(reshape2)

### for UPlst
for (i in 1:length(UPlst)){
  testgenes<- UPlst[[i]]
  CT=names(UPlst)[i]
  symbol2id=mapIds(org.Hs.eg.db,testgenes,"ENTREZID",'SYMBOL')
  id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID
  
  if(length(id) > 0){
    #GO BP 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
    ego_res <- as.data.frame(ego)
    if(nrow(ego_res)>0){
      write.xlsx(ego_res,"/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_1/UP_BP.xlsx",row.names = FALSE,sheetName = CT,append = TRUE)
    }
    
    #GO CC 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "MF", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
    ego_res <- as.data.frame(ego)
    if(nrow(ego_res)>0){
      write.xlsx(ego_res,"/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_1/UP_MF.xlsx",row.names = FALSE,sheetName = CT,append = TRUE)
    } 
    
    #GO MF 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "CC", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
    ego_res <- as.data.frame(ego)
    if(nrow(ego_res)>0){
      write.xlsx(ego_res,"/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_1/UP_CC.xlsx",row.names = FALSE,sheetName = CT,append = TRUE)
    } 
    
    #KEGG分析#
    ekk <- enrichKEGG(gene= id,organism  = 'hsa')	 #KEGG富集分析
    ekk_res <- as.data.frame(ekk)
    if(nrow(ekk_res)>0){
      write.xlsx(ekk_res,"/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_1/UP_KEGG.xlsx",row.names = FALSE,sheetName = CT,append = TRUE)
    }
  }
}

### for DNlst
for(i in 1:length(DNlst)){
  testgenes<- DNlst[[i]]
  CT=names(DNlst)[i]
  symbol2id=mapIds(org.Hs.eg.db,testgenes,"ENTREZID",'SYMBOL')
  id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID
  
  if(length(id) > 0){
    #GO BP 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
    ego_res <- as.data.frame(ego)
    if(nrow(ego_res)>0){
      write.xlsx(ego_res,"/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_1/DN_BP.xlsx",row.names = FALSE,sheetName = CT,append = TRUE)
    }
    
    #GO MF 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "CC", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
    ego_res <- as.data.frame(ego)
    if(nrow(ego_res)>0){
      write.xlsx(ego_res,"/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_1/DN_CC.xlsx",row.names = FALSE,sheetName = CT,append = TRUE)
    } 
    
    #GO CC 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "MF", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
    ego_res <- as.data.frame(ego)
    if(nrow(ego_res)>0){
      write.xlsx(ego_res,"/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_1/DN_MF.xlsx",row.names = FALSE,sheetName = CT,append = TRUE)
    } 
    
    #KEGG分析#
    ekk <- enrichKEGG(gene= id,organism  = 'hsa')	 #KEGG富集分析
    ekk_res <- as.data.frame(ekk)
    if(nrow(ekk_res)>0){
      write.xlsx(ekk_res,"/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_1/DN_KEGG.xlsx",row.names = FALSE,sheetName = CT,append = TRUE)
    }
  }
}


### 可视化 pick term for Gene>1
DEG.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inChildrenAdults.RDS')
MK.CSC <- DEG.lst$CSC; MK.C0 <- DEG.lst$TC; MK.C1 <- DEG.lst$C1; MK.C2 <- DEG.lst$C2

pickterm_UP <- c("extracellular matrix organization","collagen fibril organization",'ribosome assembly','translational initiation',
                 'ribosomal small subunit biogenesis')

pickterm_DN <- c("aging","response to oxidative stress","reactive oxygen species metabolic process","reactive oxygen species biosynthetic process",
                 "cell death in response to oxidative stress",
                 "response to interleukin-6",'regulation of inflammatory response',"intrinsic apoptotic signaling pathway",
                 'response to tumor necrosis factor')

type='DN'
if(type=='DN'){
  pickterm=pickterm_DN
  lfc.max=3
  lfc.min=0
}else if(type=='UP'){
  pickterm=pickterm_UP
  lfc.max=0
  lfc.min=-3
}
Common_BP <- read.xlsx(paste('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_1/',type,'_BP.xlsx',sep=''),sheetName = 'common')
Common_BP <- Common_BP[Common_BP$p.adjust < 0.05,]
index <- which(Common_BP$Description %in% pickterm)
pickCommon_BP <- Common_BP[index,]

pickCommon <- pickCommon_BP

pickGene <- c()
for(i in 1:nrow(pickCommon)){
  genes <- unlist(strsplit(pickCommon$geneID[i],'/'))
  pickGene <- c(pickGene,genes)
}
pickGene <- unique(pickGene)

pick_MK.CSC <- MK.CSC[rownames(MK.CSC) %in% pickGene,]
pick_MK.CSC <- data.frame(gene=rownames(pick_MK.CSC),avg_logFC=pick_MK.CSC$avg_logFC)
pick_MK.C0 <- MK.C0[rownames(MK.C0) %in% pickGene,]
pick_MK.C0 <- data.frame(gene=rownames(pick_MK.C0),avg_logFC=pick_MK.C0$avg_logFC)
pick_MK.C1 <- MK.C1[rownames(MK.C1) %in% pickGene,]
pick_MK.C1 <- data.frame(gene=rownames(pick_MK.C1),avg_logFC=pick_MK.C1$avg_logFC)
pick_MK.C2 <- MK.C2[rownames(MK.C2) %in% pickGene,]
pick_MK.C2 <- data.frame(gene=rownames(pick_MK.C2),avg_logFC=pick_MK.C2$avg_logFC)

pickGene_df <- merge(pick_MK.CSC,merge(pick_MK.C0,merge(pick_MK.C1,pick_MK.C2,by='gene'),by='gene'),by='gene')
pickGene_df <- data.frame(genes=pickGene_df$gene,logFC=apply(pickGene_df[,-1],1,mean))
pickGene_df$logFC <- -(pickGene_df$logFC)

data <- data.frame()
for(i in 1:nrow(pickCommon)){
  term=pickCommon$Description[i]
  genes=unlist(strsplit(pickCommon$geneID[i],'/'))
  df <- data.frame(
    term=rep(term,length(genes)),
    genes=genes,
    count=rep(pickCommon$Count[i],length(genes)))
  data <- rbind(data,df)
}
data <- merge(data,pickGene_df,by='genes')
data <- data[,c('term','genes','logFC')]

require("RColorBrewer")
library(GOplot)
library(ggplot2)
mycol=brewer.pal(length(unique(data$term)), "Set3")[1:length(unique(data$term))]
chord <- chord_dat(data,pickGene_df , unique(data$term))
Chord.plot <- GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5,
                      process.label = 5,lfc.max = lfc.max,lfc.min =lfc.min,ribbon.col = mycol)+
  theme(plot.margin = unit(c(0,0,-0.1,-0.1), "cm"),
        legend.key.size = unit(0.05,'cm'),legend.key.height = unit(0.05,'cm'),legend.key.width = unit(0.05,'cm'))
ggsave(paste('JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_1/sFig_Chord_',type,'.pdf',sep=''),Chord.plot,
    width = 22,height = 24,units = 'cm')

## barplot
pickCommon$log10Pval <- -log10(pickCommon$pvalue)
pickCommon <- pickCommon[order(pickCommon$log10Pval,decreasing = F),]
pickCommon$Description <- factor(pickCommon$Description, levels = pickCommon$Description)

if(type=='DN'){
  p <- ggplot(pickCommon, aes(x = Description, y = log10Pval)) + 
    geom_bar(stat = 'identity', color = '#B15E72', fill = '#B15E72',width = 0.8) + 
    theme_classic() + coord_flip() +
    labs(x='',y = expression(paste("-log"[10], "(", italic("P"), "-value)"))) +
    theme(axis.title = element_text(size = 8, colour = 'black'), 
          axis.text.y = element_text(size = 8, colour = 'black'), 
          axis.text.x = element_text(size = 8, colour = 'black'))
  ggsave(paste('JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_1/sFig_Barplot_',type,'.pdf',sep=''),p,
         height = 4, width = 10, units = 'cm')
}else{
  p <- ggplot(pickCommon, aes(x = Description, y = log10Pval)) + 
    geom_bar(stat = 'identity', color = '#7F99CE', fill = '#7F99CE',width = 0.8) + 
    theme_classic() + coord_flip() +
    labs(x='',y = expression(paste("-log"[10], "(", italic("P"), "-value)"))) +
    theme(axis.title = element_text(size = 8, colour = 'black'), 
          axis.text.y = element_text(size = 8, colour = 'black'), 
          axis.text.x = element_text(size = 8, colour = 'black'))
  ggsave(paste('JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_1/sFig_Barplot_',type,'.pdf',sep=''),p,
         height = 3, width = 8, units = 'cm')
}



######
### gene=4 富集分析
######
commonUP <- names(up_sum)[which(up_sum == 4)]
commonDN <- names(dn_sum)[which(dn_sum == 4)]
commonLst <- list(UP=commonUP,DN=commonDN)
for(i in 1:length(commonLst)){
  type=names(commonLst)[i]
  testgenes<- commonLst[[i]]
  CT=names(commonLst)[i]
  symbol2id=mapIds(org.Hs.eg.db,testgenes,"ENTREZID",'SYMBOL')
  id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID
  
  if(length(id) > 0){
    #GO BP 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
    ego_res <- as.data.frame(ego)
    if(nrow(ego_res)>0){
      write.xlsx(ego_res,"/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_4/Common_BP.xlsx",row.names = FALSE,sheetName = type,append = TRUE)
    }
    
    #GO MF 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "CC", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
    ego_res <- as.data.frame(ego)
    if(nrow(ego_res)>0){
      write.xlsx(ego_res,"/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_4/Common_CC.xlsx",row.names = FALSE,sheetName = type,append = TRUE)
    } 
    
    #GO CC 富集分析#
    ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "MF", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
    ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
    ego_res <- as.data.frame(ego)
    if(nrow(ego_res)>0){
      write.xlsx(ego_res,"/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_4/Common_MF.xlsx",row.names = FALSE,sheetName = type,append = TRUE)
    } 
    
    #KEGG分析#
    ekk <- enrichKEGG(gene= id,organism  = 'hsa')	 #KEGG富集分析
    ekk_res <- as.data.frame(ekk)
    if(nrow(ekk_res)>0){
      write.xlsx(ekk_res,"/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_4/Common_KEGG.xlsx",row.names = FALSE,sheetName = type,append = TRUE)
    }
  }
}



### 可视化 pick term for Gene=4
DEG.lst <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/DEGs_inChond_inChildrenAdults.RDS')
MK.CSC <- DEG.lst$CSC; MK.C0 <- DEG.lst$C0; MK.C1 <- DEG.lst$C1; MK.C2 <- DEG.lst$C2

pickterm_UP <- c("extracellular matrix structural constituent","extracellular matrix binding")

pickterm_DN <- c("aging","response to oxidative stress","reactive oxygen species metabolic process","reactive oxygen species biosynthetic process",
                 "cell death in response to oxidative stress",
                 "response to interleukin-6",'regulation of inflammatory response',"intrinsic apoptotic signaling pathway",
                 'cellular response to tumor necrosis factor')

type='UP'
if(type=='DN'){
  pickterm=pickterm_DN
  lfc.max=3
  lfc.min=0
}else if(type=='UP'){
  pickterm=pickterm_UP
  lfc.max=0
  lfc.min=-3
}
Common_BP <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_4/Common_BP_all.xlsx',sheetName = type)
Common_BP <- Common_BP[Common_BP$p.adjust < 0.05,]
index <- which(Common_BP$Description %in% pickterm)
pickCommon_BP <- Common_BP[index,]

Common_MF <- read.xlsx('/home/yzj/JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_4/Common_MF_all.xlsx',sheetName = type)
Common_MF <- Common_MF[Common_MF$p.adjust < 0.05,]
index <- which(Common_MF$Description %in% pickterm)
pickCommon_MF<- Common_MF[index,]

pickCommon <- rbind(pickCommon_BP,pickCommon_MF)

pickGene <- c()
for(i in 1:nrow(pickCommon)){
  genes <- unlist(strsplit(pickCommon$geneID[i],'/'))
  pickGene <- c(pickGene,genes)
}
pickGene <- unique(pickGene)

pick_MK.CSC <- MK.CSC[rownames(MK.CSC) %in% pickGene,]
pick_MK.CSC <- data.frame(gene=rownames(pick_MK.CSC),avg_logFC=pick_MK.CSC$avg_logFC)
pick_MK.C0 <- MK.C0[rownames(MK.C0) %in% pickGene,]
pick_MK.C0 <- data.frame(gene=rownames(pick_MK.C0),avg_logFC=pick_MK.C0$avg_logFC)
pick_MK.C1 <- MK.C1[rownames(MK.C1) %in% pickGene,]
pick_MK.C1 <- data.frame(gene=rownames(pick_MK.C1),avg_logFC=pick_MK.C1$avg_logFC)
pick_MK.C2 <- MK.C2[rownames(MK.C2) %in% pickGene,]
pick_MK.C2 <- data.frame(gene=rownames(pick_MK.C2),avg_logFC=pick_MK.C2$avg_logFC)

pickGene_df <- merge(pick_MK.CSC,merge(pick_MK.C0,merge(pick_MK.C1,pick_MK.C2,by='gene'),by='gene'),by='gene')
pickGene_df <- data.frame(genes=pickGene_df$gene,logFC=apply(pickGene_df[,-1],1,mean))
pickGene_df$logFC <- -(pickGene_df$logFC)

data <- data.frame()
for(i in 1:nrow(pickCommon)){
  term=pickCommon$Description[i]
  genes=unlist(strsplit(pickCommon$geneID[i],'/'))
  df <- data.frame(
    term=rep(term,length(genes)),
    genes=genes,
    count=rep(pickCommon$Count[i],length(genes)))
  data <- rbind(data,df)
}
data <- merge(data,pickGene_df,by='genes')
data <- data[,c('term','genes','logFC')]

require("RColorBrewer")
mycol=brewer.pal(length(unique(data$term)), "Set3")[1:length(unique(data$term))]
chord <- chord_dat(data,pickGene_df , unique(data$term))
Chord.plot <- GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 8,
                      process.label = 8,lfc.max = lfc.max,lfc.min =lfc.min,ribbon.col = mycol)
pdf(paste('JingMA_NEW/res/compControl/ChildrenvsAdults/ClusterPro_sep/FC1.5_adjP0.05/Gene_4/Chord_',type,'.pdf',sep=''),width = 13,height = 14)
print(Chord.plot)
dev.off()


##################
###补充材料6：Geneset score的boxplot
##################
library(xlsx)
library(ggsignif)
library(ggpubr)
pbmc_C <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
pbmc_C$Phase <- factor(pbmc_C$Phase,levels = c('Children','Adults'))

geneSet <- read.xlsx('JingMA_NEW/res/compControl/ChildrenvsAdults/geneset_skinAging.xlsx',sheetIndex = 4)
geneSet <- geneSet[,c('DNA.repair.genes','DNA.damage.genes',
                      'Inflammatory.response.genes','ROS.genes',
                      'NF.κB.pathway.genes.')]


for(i in 1:length(colnames(geneSet))){
  set.name=colnames(geneSet)[i]
  print(set.name)
  pbmc.score <- AddModuleScore(pbmc_C,features = list(na.omit(geneSet[,set.name])),name = 'score')
  dt <- data.frame(cells = colnames(pbmc.score),score=pbmc.score$score1,celltype=pbmc.score$celltype,phase=pbmc.score$Phase)
  
  p <- ggplot(dt, aes(x=celltype, y=score,fill=phase)) + 
    geom_boxplot(width=0.5,position=position_dodge(0.9),outlier.size = 0.05)+
    scale_fill_manual(values = c( "#6C6C6C","#DC143C"))+ 
    theme_bw()+
    theme(axis.text.x=element_blank(),axis.ticks.x =element_blank(),
          axis.text.y=element_text(size=5,colour="black"),axis.title.y=element_text(size = 5,colour="black"), 
          axis.ticks.y =element_line(colour="black",size = 0.01),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 5),
          plot.margin = unit(c(0.1,0,-0.3,-0.2), "cm"))+ 
    ylab("")+xlab("")+ labs(title = set.name)+
    facet_wrap(~celltype,ncol = 4,scales= "free_x")+
    theme(strip.background = element_rect(color = "black", fill = "white",size = 0.5),
          strip.text.x = element_text(size = 5, color = "black",face = 'bold'),
          panel.grid = element_blank(),panel.border = element_rect(color = 'black',size = 0.3))+
    stat_compare_means(method = 'wilcox.test',label = "p.signif",label.x=0.1)+guides(fill=FALSE)
  ggsave(paste('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/sFig_',set.name,'.pdf',sep=''),p,width = 3,height = 1.5)
}


SASP <- read.gmt('publicData/GMT/SASP_reactome.gmt')[,'gene']
pbmc.score <- AddModuleScore(pbmc_C,features = list(SASP),name = 'score')
dt <- data.frame(cells = colnames(pbmc.score),score=pbmc.score$score1,celltype=pbmc.score$celltype,phase=pbmc.score$Phase)

p <- ggplot(dt, aes(x=celltype, y=score,fill=phase)) + 
  geom_boxplot(width=0.5,position=position_dodge(0.9),outlier.size = 0.05)+
  scale_fill_manual(values = c( "#6C6C6C","#DC143C"))+ 
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks.x =element_blank(),
        axis.text.y=element_text(size=5,colour="black"),axis.title.y=element_text(size = 5,colour="black"), 
        axis.ticks.y =element_line(colour="black",size = 0.01),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 5),
        plot.margin = unit(c(0.1,0,-0.3,-0.2), "cm"))+ 
  ylab("")+xlab("")+ labs(title = 'SASP')+
  facet_wrap(~celltype,ncol = 4,scales= "free_x")+
  theme(strip.background = element_rect(color = "black", fill = "white",size = 0.5),
        strip.text.x = element_text(size = 5, color = "black",face = 'bold'),
        panel.grid = element_blank(),panel.border = element_rect(color = 'black',size = 0.3))+
  stat_compare_means(method = 'wilcox.test',label = "p.signif",label.x=0.1)+guides(fill=FALSE)
ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG_new/sFig_SASP.pdf',p,width = 3,height = 1.5)


c5_GO <- read.gmt('publicData/GMT/c5.all.v6.2.symbols.gmt')
Aging_GO <- c5_GO$gene[c5_GO$term=='GO_AGING']
pbmc.score <- AddModuleScore(pbmc_C,features = list(Aging_GO),name = 'score')
dt <- data.frame(cells = colnames(pbmc.score),score=pbmc.score$score1,celltype=pbmc.score$celltype,phase=pbmc.score$Phase)
dt$celltype <- as.character(dt$celltype)
dt$CT=dt$celltype
dt$celltype[dt$CT=='CSC']='CSPC'
dt$celltype[dt$CT=='C0']='EC'
dt$celltype[dt$CT=='C1']='IC'
dt$celltype[dt$CT=='C2']='LC'


p <- ggplot(dt, aes(x=celltype, y=score,fill=phase)) + 
  geom_boxplot(width=0.5,position=position_dodge(0.9),outlier.size = 0.05)+
  scale_fill_manual(values = c( "#6C6C6C","#DC143C"))+ 
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks.x =element_blank(),
        axis.text.y=element_text(size=5,colour="black"),axis.title.y=element_text(size = 5,colour="black"), 
        axis.ticks.y =element_line(colour="black",size = 0.01),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 5),
        plot.margin = unit(c(0.1,0,-0.3,-0.2), "cm"))+ 
  ylab("")+xlab("")+ labs(title = 'Aging')+
  facet_wrap(~celltype,ncol = 4,scales= "free_x")+
  theme(strip.background = element_rect(color = "black", fill = "white",size = 0.5),
        strip.text.x = element_text(size = 5, color = "black",face = 'bold'),
        panel.grid = element_blank(),panel.border = element_rect(color = 'black',size = 0.3))+
  stat_compare_means(method = 'wilcox.test',label = "p.signif",label.x=0.1)+guides(fill=FALSE)
ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG_new/sFig_Aging.pdf',p,width = 3,height = 1.5)



##################
###补充材料7：CEBPB的共表达
##################
####### 补充材料: 计算CEBPB与其他下调的TF的共表达
library(ggplot2)
library(SCENIC)
library(AUCell)
library(ggrepel)
library(Seurat)

regulon2Targets <- readRDS('/home/yzj/JingMA_NEW/res/SCENIC_main/int/2.5_regulonTargetsInfo.Rds')
regulonAUC <- readRDS(file='/home/yzj/JingMA_NEW/res/SCENIC_main/int/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

TF <- apply(as.matrix(rownames(regulonAUC)),1,function(x) {unlist(strsplit(x,' '))[1]})
TF <- apply(as.matrix(TF),1,function(x) {unlist(strsplit(x,'_extend'))[1]})

dnValues_mtx <- readRDS('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/DEGHeatmap_DNmtx.RDS')
dnGene <- rownames(dnValues_mtx)

dnTF <- intersect(TF,dnGene)
dnTF <- setdiff(dnTF,'CEBPB')
dnTF <- unique(regulon2Targets$TF[regulon2Targets$TF %in% dnTF & regulon2Targets$highConfAnnot =='TRUE'])

pbmc_C <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
pbmc_C$celltype <- factor(pbmc_C$celltype,levels = c('CSC','C0','C1','C2'))

children.cell <- colnames(pbmc_C)[pbmc_C$Phase=='Children']
adults.cell <- colnames(pbmc_C)[pbmc_C$Phase=='Adults']

## EXP
mat.exp <- as.data.frame(pbmc_C@assays$RNA@data)
children.exp <- mat.exp[,colnames(mat.exp) %in% children.cell]
adults.exp <- mat.exp[,colnames(mat.exp) %in% adults.cell]

CEBPB.children.exp <- as.numeric(children.exp[rownames(children.exp)=='CEBPB',])
CEBPB.adults.exp <- as.numeric(adults.exp[rownames(adults.exp)=='CEBPB',])

children.cor <- c()
adults.cor <- c()
for(i in 1:length(dnTF)){
  tf=dnTF[i]
  tf.children.exp <- as.numeric(children.exp[rownames(children.exp)==tf,])
  tf.adults.exp <- as.numeric(adults.exp[rownames(adults.exp)==tf,])
  
  cor <- cor.test(CEBPB.children.exp,tf.children.exp,method = 'spearman')
  children.cor[i] <- cor$estimate
  
  cor <- cor.test(CEBPB.adults.exp,tf.adults.exp,method = 'spearman')
  adults.cor[i] <- cor$estimate
}

mat.cor <- data.frame(TF=dnTF,children.cor=children.cor,adults.cor=adults.cor)
mat.cor$cor=mat.cor$children.cor+mat.cor$adults.cor
mat.cor[order(mat.cor$cor,decreasing = T),]

plot.cor <- ggplot(mat.cor) +
  geom_point(aes(children.cor, adults.cor,size=cor), color = 'firebrick3') +
  geom_text_repel(aes(children.cor, adults.cor, label = TF)) +
  theme_classic()+
  scale_size_continuous(range = c(0,3)) +
  theme(axis.title = element_text(size=6,colour = 'black'),axis.text = element_text(size=6,colour = 'black'),
        legend.text = element_text(size=6,colour = 'black'),legend.title  = element_text(size=6,colour = 'black'),
        legend.key.size = unit(0.1,'cm'),legend.key.height = unit(0.1,'cm'),legend.key.width = unit(0.1,'cm'))+
  labs(x='Correlation in children',y='Correlation in adults')
ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/sFig_Correaltion_EXP.pdf',plot.cor,width = 11,height = 9,units = 'cm')



## AUC activity
mat.auc <- getAUC(regulonAUC)[,colnames(pbmc_C)]
children.auc <- mat.auc[,colnames(mat.auc) %in% children.cell]
adults.auc <- mat.auc[,colnames(mat.auc) %in% adults.cell]

CEBPB.children.auc <- as.numeric(children.auc[rownames(children.auc)=='CEBPB (258g)',])
CEBPB.adults.auc <- as.numeric(adults.auc[rownames(adults.auc)=='CEBPB (258g)',])

dnTF.regulon <- c()
for(i in 1:length(dnTF)){
  tf=dnTF[i]
  dnTF.regulon <- c(dnTF.regulon,rownames(mat.auc)[grep(tf,rownames(mat.auc))])
}

children.cor <- c()
adults.cor <- c()
for(i in 1:length(dnTF.regulon)){
  tf=dnTF.regulon[i]
  tf.children.auc <- as.numeric(children.auc[rownames(children.auc)==tf,])
  tf.adults.auc <- as.numeric(adults.auc[rownames(adults.auc)==tf,])
  
  cor <- cor.test(CEBPB.children.auc,tf.children.auc,method = 'spearman')
  children.cor[i] <- cor$estimate
  
  cor <- cor.test(CEBPB.adults.auc,tf.adults.auc,method = 'spearman')
  adults.cor[i] <- cor$estimate
}

mat.cor <- data.frame(TF=dnTF.regulon,children.cor=children.cor,adults.cor=adults.cor)
mat.cor$cor=mat.cor$children.cor+mat.cor$adults.cor
mat.cor[order(mat.cor$cor,decreasing = T),]

plot.cor <- ggplot(mat.cor) +
  geom_point(aes(children.cor, adults.cor,size=cor), color = 'firebrick3') +
  geom_text_repel(aes(children.cor, adults.cor, label = TF)) +
  theme_classic()+
  scale_size_continuous(range = c(0,3)) +
  theme(axis.title = element_text(size=6,colour = 'black'),axis.text = element_text(size=6,colour = 'black'),
        legend.text = element_text(size=6,colour = 'black'),legend.title  = element_text(size=6,colour = 'black'),
        legend.key.size = unit(0.1,'cm'),legend.key.height = unit(0.1,'cm'),legend.key.width = unit(0.1,'cm'))+
  labs(x='Correlation in children',y='Correlation in adults')
ggsave('JingMA_NEW/res/compControl/ChildrenvsAdults/FIG/sFig_Correaltion_AUC.pdf',plot.cor,width = 11,height = 9,units = 'cm')




##################
###补充材料： 计算UMI/Gene
##################
pbmc_C <- readRDS('JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype_Control_Chond.Rdata')
data <- data.frame(celltype=pbmc_C$celltype,UMI=pbmc_C$nCount_RNA,Gene=pbmc_C$nCount_RNA)
library(ggplot2)
p.UMI <- ggplot(data=data,aes(x=celltype,y=UMI,fill=celltype))+geom_boxplot()+
  theme_classic()

p.Gene <- ggplot(data=data,aes(x=celltype,y=Gene,fill=celltype))+geom_boxplot()+
  theme_classic()

library(gridExtra)
p <- grid.arrange(p.UMI,p.Gene,nrow = 2, ncol = 1)
ggsave('JingMA_NEW/res/Harmony/ALL/FIG/sFig_UMI_Gene_bycelltype.pdf',p,width = 6,height = 6)
