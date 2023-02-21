############
## 3. TF Enrichment
############
.libPaths(c("/home/yzj/R/x86_64-pc-linux-gnu-library/4.0","/home/zy/tools/R-4.0.0/library"))
library(pheatmap)
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
library(xlsx)

type='main'
op=paste('/home/yzj/JingMA_NEW/res/SCENIC_',type,sep='')
regulonGeneSet <-  readRDS(paste('/home/yzj/JingMA_NEW/res/SCENIC_',type,'/int/2.6_regulons_asGeneSet.Rds',sep=''))
TF <- readRDS(paste('/home/yzj/JingMA_NEW/res/SCENIC_',type,'/FIG/heatmap.RDS',sep=''))
TF <- apply(as.matrix(rownames(TF)),1,function(x) unlist(strsplit(x," "))[1])

for(i in 1:length(regulonGeneSet)){
  genes <- regulonGeneSet[[i]]
  tf <- names(regulonGeneSet)[i]
  print(tf)
  if(tf %in% TF){
    testgenes<- genes
    symbol2id=mapIds(org.Hs.eg.db,testgenes,"ENTREZID",'SYMBOL')
    id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID
    
    if(length(id) > 0){
      #GO BP 富集分析#
      ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
      ego <- clusterProfiler::simplify(ego,cutoff=0.8,by="p.adjust",select_fun=min,measure="Wang")
      ego_res <- as.data.frame(ego)
      if(nrow(ego_res)>0){
        write.xlsx(ego_res,paste0(op,'/GO/',tf,"_BP.xlsx",sep=''),row.names = FALSE, sheetName = tf)
      }
     
      #KEGG分析#
      ekk <- enrichKEGG(gene= id,organism  = 'hsa')	 #KEGG富集分析
      ekk_res <- as.data.frame(ekk)
      if(nrow(ekk_res) > 0){
        write.xlsx(ekk_res,paste0(op,'/KEGG/',tf,"_KEGG.xlsx"),row.names = FALSE, sheetName = tf )
      }
    }
  }
}


#########
###Fig 2k 挑选term组合画图
#########
library(xlsx)
library(ggplot2)

FOS <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/FOS_extended_BP.xlsx',sheetIndex = 1)
index <- c(17,19,32,51,58,81,180,218,241)
pickFOS <- FOS[index,]
pickFOS$Regulon <- 'FOS'
print(pickFOS[,c(2,8)])


KLF10 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/KLF10_BP.xlsx',sheetIndex = 1)
index <- c(12,25)
pickKLF10 <- KLF10[index,]
pickKLF10$Regulon <- 'KLF10'
print(pickKLF10[,c(2,8)])


DDIT3 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/DDIT3_extended_BP.xlsx',sheetIndex = 1)
index <- c(15,16,35,50,129,131,133,156,161,235,280,281,300,376)
pickDDIT3 <- DDIT3[index,]
pickDDIT3$Regulon <- 'DDIT3'
print(pickDDIT3[,c(2,8)])


FOXC1 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/FOXC1_BP.xlsx',sheetIndex = 1)
index <- c(2,4,9,21,25,31,64,72,76,87,146,151,271,362,530,584)
pickFOXC1 <- FOXC1[index,]
pickFOXC1$Regulon <- 'FOXC1'
print(pickFOXC1[,c(2,8)])

FOXO1 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/FOXO1_BP.xlsx',sheetIndex = 1)
index <- c(4,6,12,15,33,55,76,85,109,110,125,176,188,234,369)
pickFOXO1 <- FOXO1[index,]
pickFOXO1$Regulon <- 'FOXO1'
print(pickFOXO1[,c(2,8)])


HES1 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/HES1_BP.xlsx',sheetIndex = 1)
index <- c(12,26,30,85,103)
pickHES1 <- HES1[index,]
pickHES1$Regulon <- 'HES1'
print(pickHES1[,c(2,8)])

SOX5 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/SOX5_BP.xlsx',sheetIndex = 1)
index <- c(2,5,21,22,26,72,98,99,105)
pickSOX5 <- SOX5[index,]
pickSOX5$Regulon <- 'SOX5'
print(pickSOX5[,c(2,8)])


JUN <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/JUN_BP.xlsx',sheetIndex = 1)
index <- c(25)
pickJUN <- JUN[index,]
pickJUN$Regulon <- 'JUN'
print(pickJUN[,c(2,8)])


ATF3 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/ATF3_BP.xlsx',sheetIndex = 1)
index <- c(36,49,60)
pickATF3<- ATF3[index,]
pickATF3$Regulon <- 'ATF3'
print(pickATF3[,c(2,8)])


JUNB <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/JUNB_BP.xlsx',sheetIndex = 1)
index <- c(6,21,87,186)
pickJUNB <- JUNB[index,]
pickJUNB$Regulon <- 'JUNB'
print(pickJUNB[,c(2,8)])


EGR1 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/EGR1_BP.xlsx',sheetIndex = 1)
index <- c(4,6,10,36,62,67,215,270,344,403,497)
pickEGR1 <- EGR1[index,]
pickEGR1$Regulon <- 'EGR1'
print(pickEGR1[,c(2,8)])


JUND <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/JUND_BP.xlsx',sheetIndex = 1)
index <- c(1,7)
pickJUND <- JUND[index,]
pickJUND$Regulon <- 'JUND'
print(pickJUND[,c(2,8)])


SMAD7 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/SMAD7_BP.xlsx',sheetIndex = 1)
index <- c(1,3,5,13,45,68)
pickSMAD7 <- SMAD7[index,]
pickSMAD7$Regulon <- 'SMAD7'
print(pickSMAD7[,c(2,8)])

SOX8 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/SOX8_BP.xlsx',sheetIndex = 1)
index <- c(4,10,18,19,22,23,87,91,112,157,160,212)
pickSOX8 <- SOX8[index,]
pickSOX8$Regulon <- 'SOX8'
print(pickSOX8[,c(2,8)])

FOXL1 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/FOXL1_extended_BP.xlsx',sheetIndex = 1)
index <- c(6,7,9,13,59,65,74,79,106,107,137,214,324)
pickFOXL1 <- FOXL1[index,]
pickFOXL1$Regulon <- 'FOXL1'
print(pickFOXL1[,c(2,8)])


NFATC2 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/NFATC2_BP.xlsx',sheetIndex = 1)
index <- c(1,3,7,8,27,29,33,86,165,188)
pickNFATC2 <- NFATC2[index,]
pickNFATC2$Regulon <- 'NFATC2'
print(pickNFATC2[,c(2,8)])


FOXA2 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/FOXA2_extended_BP.xlsx',sheetIndex = 1)
index <- c(5,12,15,20,23,40,49,50,57,120,125,135,147,157,278)
pickFOXA2 <- FOXA2[index,]
pickFOXA2$Regulon <- 'FOXA2'
print(pickFOXA2[,c(2,8)])


FOXC2 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/FOXC2_BP.xlsx',sheetIndex = 1)
index <- c(6,13,14,23,24,41,50,53,93,104,105,206,207,255,427)
pickFOXC2 <- FOXC2[index,]
pickFOXC2$Regulon <- 'FOXC2'
print(pickFOXC2[,c(2,8)])


FOXD1 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/FOXD1_extended_BP.xlsx',sheetIndex = 1)
index <- c(9,15,17,31,40,84,103,109,147,200,224,230,264,270,287,302,308,355,460)
pickFOXD1 <- FOXD1[index,]
pickFOXD1$Regulon <- 'FOXD1'
print(pickFOXD1[,c(2,8)])


HMGB2 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/HMGB2_BP.xlsx',sheetIndex = 1)
index <- c(15,18,19,24,73,97,163,171,209,251,316,332,418,423,452,530)
pickHMGB2 <- HMGB2[index,]
pickHMGB2$Regulon <- 'HMGB2'
print(pickHMGB2[,c(2,8)])


PITX1 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/PITX1_extended_BP.xlsx',sheetIndex = 1)
index <- c(14:16,48,91,267)
pickPITX1 <- PITX1[index,]
pickPITX1$Regulon <- 'PITX1'
print(pickPITX1[,c(2,8)])


PPARG <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/PPARG_BP.xlsx',sheetIndex = 1)
index <- c(2)
pickPPARG <- PPARG[index,]
pickPPARG$Regulon <- 'PPARG'
print(pickPPARG[,c(2,8)])


SMAD3 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/SMAD3_BP.xlsx',sheetIndex = 1)
index <- c(1,2)
pickSMAD3 <- SMAD3[index,]
pickSMAD3$Regulon <- 'SMAD3'
print(pickSMAD3[,c(2,8)])


HIVEP3 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/HIVEP3_BP.xlsx',sheetIndex = 1)
index <- c(5,38,77,78,116,136,228,281,330,352)
pickHIVEP3 <- HIVEP3[index,]
pickHIVEP3$Regulon <- 'HIVEP3'
print(pickHIVEP3[,c(2,8)])


KLF3 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/KLF3_BP.xlsx',sheetIndex = 1)
index <- c(1,12,23,25,26)
pickKLF3 <- KLF3[index,]
pickKLF3$Regulon <- 'KLF3'
print(pickKLF3[,c(2,8)])

ELF3 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/ELF3_BP.xlsx',sheetIndex = 1)
index <- c(1:3,19,59,85,237,277)
pickELF3<- ELF3[index,]
pickELF3$Regulon <- 'ELF3'
print(pickELF3[,c(2,8)])

HMGA2 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/HMGA2_BP.xlsx',sheetIndex = 1)
index <- c(2,14,45,61,71,145,156,208,285)
pickHMGA2 <- HMGA2[index,]
pickHMGA2$Regulon <- 'HMGA2'
print(pickHMGA2[,c(2,8)])


SOX9 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/SOX9_BP.xlsx',sheetIndex = 1)
index <- c(1,5,10)
pickSOX9 <- SOX9[index,]
pickSOX9$Regulon <- 'SOX9'
print(pickSOX9[,c(2,8)])


MSX1 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/MSX1_extended_BP.xlsx',sheetIndex = 1)
index <- c(3,10,94,97,112,131,148,160,178)
pickMSX1 <- MSX1[index,]
pickMSX1$Regulon <- 'MSX1'
print(pickMSX1[,c(2,8)])


NFIL3 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/NFIL3_extended_BP.xlsx',sheetIndex = 1)
index <- c(12,17,60,151,164,173)
pickNFIL3 <- NFIL3[index,]
pickNFIL3$Regulon <- 'NFIL3'
print(pickNFIL3[,c(2,8)])


CREB3L1 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/CREB3L1_extended_BP.xlsx',sheetIndex = 1)
index <- c(1,4,5,71,130,215,221,294,378)
pickCREB3L1 <- CREB3L1[index,]
pickCREB3L1$Regulon <- 'CREB3L1'
print(pickCREB3L1[,c(2,8)])


BHLHE22 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/BHLHE22_extended_BP.xlsx',sheetIndex = 1)
index <- c(1:4,6,11,28,29,30)
pickBHLHE22 <- BHLHE22[index,]
pickBHLHE22$Regulon <- 'BHLHE22'
print(pickBHLHE22[,c(2,8)])


EN2 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/EN2_BP.xlsx',sheetIndex = 1)
index <- c(14,32,52,61,77,97,118,181,188,207,234,239,242,327,363)
pickEN2 <- EN2[index,]
pickEN2$Regulon <- 'EN2'
print(pickEN2[,c(2,8)])


FOXP1 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/FOXP1_BP.xlsx',sheetIndex = 1)
index <- c(1,2,10,132,133,189,215,345,443)
pickFOXP1 <- FOXP1[index,]
pickFOXP1$Regulon <- 'FOXP1'
print(pickFOXP1[,c(2,8)])


IRX2 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/IRX2_BP.xlsx',sheetIndex = 1)
index <- c(5,11,13,19,45,158,169)
pickIRX2 <- IRX2[index,]
pickIRX2$Regulon <- 'IRX2'
print(pickIRX2[,c(2,8)])


LMO2 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/LMO2_extended_BP.xlsx',sheetIndex = 1)
index <- c(1,2,8,22,24,32,72,95,96,99)
pickLMO2 <- LMO2[index,]
pickLMO2$Regulon <- 'LMO2'
print(pickLMO2[,c(2,8)])


REL <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/REL_BP.xlsx',sheetIndex = 1)
index <- c(1,4,5,10)
pickREL <- REL[index,]
pickREL$Regulon <- 'REL'
print(pickREL[,c(2,8)])


CEBPB <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/CEBPB_BP.xlsx',sheetIndex = 1)
index <- c(1,2,4,6,55,90,276,305)
pickCEBPB <- CEBPB[index,]
pickCEBPB$Regulon <- 'CEBPB'
print(pickCEBPB[,c(2,8)])


EZH2 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/EZH2_extended_BP.xlsx',sheetIndex = 1)
index <- c(1,2,6)
pickEZH2 <- EZH2[index,]
pickEZH2$Regulon <- 'EZH2'
print(pickEZH2[,c(2,8)])

TWIST1 <- read.xlsx('JingMA_NEW/res/SCENIC_main/GO/TWIST1_BP.xlsx',sheetIndex = 1)
index <- c(1:4,18,23,111,114,179,214)
pickTWIST1 <- TWIST1[index,]
pickTWIST1$Regulon <- 'TWIST1'
print(pickTWIST1[,c(2,8)])

TF.lst <- c()

############
## Fig 2K
############
TFlst <- rev(c('SOX8','NFATC2','SOX9','FOXP1'))
Term <- rbind(pickSOX8[1:3,],pickNFATC2[1:3,],pickSOX9[1:2,],pickFOXP1[1:3,])
#Term <- Term[!duplicated(Term$Description),]
Term$log_adjP <- -log(Term$p.adjust)

sdata=split(Term,Term$Regulon)
result=lapply(sdata,function(x) x[order(x[,11]),])
result <- rbind(result$FOXP1,result$SOX9,result$NFATC2,result$SOX8)
result$Regulon <- factor(result$Regulon,levels = TFlst)
result$newTerm <- paste(result$Regulon,result$Description,sep='_')
result$newTerm <- factor(result$newTerm,levels = result$newTerm)

require("RColorBrewer")
library(ggplot2)
mycol <- brewer.pal(n = length(result), name = "Set3")
plot.bar <- 
  ggplot(data = result, aes(x = newTerm, y = log_adjP, color = Regulon, fill = Regulon)) + 
  geom_bar(stat = 'identity',width = 0.8) + 
  theme_classic() + coord_flip() +
  scale_x_discrete(breaks = paste(result$Regulon, result$Description, sep = '_'), 
                   labels = result$Description) + 
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol) +
  labs(title = 'Enriched GO terms',x='', y = expression(paste("-log"[10], "(", italic("P"), "-value)"))) +
  theme(axis.title = element_text(size = 8,color = 'black'), 
        axis.text = element_text(size = 8,color = 'black'),
        legend.text = element_text(size = 8,color = 'black'),
        legend.key.size = unit(0.1,'cm'),legend.key.width  = unit(0.1,'cm'),legend.key.height = unit(0.1,'cm'),
        legend.title = element_text(size = 8, color = 'black'),legend.position = 'bottom')

ggsave('JingMA_NEW/res/SCENIC_main/FIG/Fig2K_GO_batplot.pdf',plot.bar,width = 8,height = 7,unit='cm')




############
## 补充材料
############
### Barplot no SOX8 and FOXP1
TFlst <- rev(c('EGR1','HES1','SOX5','HMGA2','HIVEP3','PPARG',
               'TWIST1','CEBPB','ETS2'))
Term <- rbind(pickJUN[1:3,],pickHMGA2[1:3,],pickCEBPB[1:3,],
              pickTWIST1[1:3,],pickPPARG[1:3,],pickHIVEP3[1:3,],pickSOX5[1:3,],
              pickHES1[1:3,],pickEGR1[1:3,])
#Term <- Term[!duplicated(Term$Description),]
Term$log_adjP <- -log(Term$p.adjust)

sdata=split(Term,Term$Regulon)
result=lapply(sdata,function(x) x[order(x[,11]),])
result <- rbind(result$JUN,result$HMGA2,result$CEBPB,result$TWIST1,result$PPARG,result$HIVEP3,
                result$SOX5,result$HES1,result$EGR1)
result$Regulon <- factor(result$Regulon,levels = TFlst)
result$newTerm <- paste(result$Regulon,result$Description,sep='_')
result$newTerm <- factor(result$newTerm,levels = result$newTerm)

require("RColorBrewer")
library(ggplot2)
mycol <- brewer.pal(n = length(result), name = "Set3")
plot.bar <- 
  ggplot(data = result, aes(x = newTerm, y = log_adjP, color = Regulon, fill = Regulon)) + 
  geom_bar(stat = 'identity') + 
  theme_classic() + coord_flip() +
  scale_x_discrete(breaks = paste(result$Regulon, result$Description, sep = '_'), 
                   labels = result$Description) + 
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol) +
  labs(x = 'Enriched GO terms', y = expression(paste("-log"[10], "(", italic("P"), "-value)"))) +
  theme(axis.title = element_text(size = 25, face='bold',color = 'black'), 
        axis.text.y = element_text(size = 25, face='bold',color = 'black'), 
        axis.text.x = element_text(size = 25, face='bold',color = 'black'),
        legend.text = element_text(size = 20, face='bold',color = 'black'),
        legend.title = element_text(size = 20, face='bold',color = 'black'))

pdf('JingMA_NEW/res/SCENIC_main/FIG/sFig_GO_batplot.pdf',width = 12,height = 12)
print(plot.bar)
dev.off()

