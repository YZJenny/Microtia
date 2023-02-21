setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/zy/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(patchwork)
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(tidyverse)
library(SCENIC)
library(igraph)
library(ggraph)
library(tidygraph)

path.TF.net <- '/home/disk/drizzle/wgk/TF_net/'

fc.cutoff <- 0.4
path.M123 <- '/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/'
path.cutoff <- paste0(path.M123, 'cutoff_', fc.cutoff, '/')
file.go.BP <- paste0(path.cutoff, 'GO_BP_all.Rdata')
list.go.BP <- readRDS(file.go.BP)
file.go.MF <- paste0(path.cutoff, 'GO_MF_all.Rdata')
list.go.MF <- readRDS(file.go.MF)

#####################
####### down
#####################
file_CSC_down <- paste0(path.TF.net, 'igraph_CSC_down.RDS')
igraph_down_2 <- readRDS(file_CSC_down)

sel.MF_SCS_down <- c('antioxidant activity',
                     'extracellular matrix structural constituent',
                     'S100 protein binding')
sel.BP_SCS_down <- c('cellular response to zinc ion',
                     'cartilage development')

BP_SCS_down <- list.go.BP[['CSC_Microtia_decrease']]
rownames(BP_SCS_down) <- BP_SCS_down$Description
MF_SCS_down <- list.go.MF[['CSC_Microtia_decrease']]
rownames(MF_SCS_down) <- MF_SCS_down$Description
df_GO <- rbind(BP_SCS_down[sel.BP_SCS_down, c('Description', 'geneID')],
               MF_SCS_down[sel.MF_SCS_down, c('Description', 'geneID')])

node_type <- V(igraph_down)$name
list_gene_GO <- list()
for (gene in node_type) {
    vec_GO <- c()
    for (term in df_GO$Description) {
        gene.set <- strsplit(df_GO$geneID[df_GO$Description == term], split = '/')[[1]]
        if (gene %in% gene.set) {
            vec_GO <- c(vec_GO, term)
        }
    }
    list_gene_GO[[gene]] <- vec_GO
}


ggraph_CSC_down <- as_tbl_graph(igraph_down_2)

plot_CSC_down <-
    # ggraph(ggraph_CSC_down, layout = 'centrality',cent=page_rank) +
    ggraph(ggraph_CSC_down, layout = 'stress') +
    geom_edge_link(aes(edge_width=weight, alpha = weight),color="gray",
                   arrow = arrow(length = unit(2, 'mm')),
                   end_cap = circle(1.5, 'mm'), 
                   start_cap = circle(3, 'mm')) +
    scale_edge_width(range=c(0.5,1)) +
    scale_edge_alpha(range=c(0.2,1)) + 
    scale_size_continuous(range = c(2,10)) +
    geom_node_point(aes(size = page_rank, fill = group),
                    shape=21, color = 'transparent') +
    # scale_color_manual(values = c('#0000CD', '#DC143C', 'dimgray')) +
    scale_fill_manual(values = c('#4169E1', '#FF4500', 'gray')) +
    #scale_alpha_manual(values = c(1, 1, 0.2)) +
    geom_node_text(aes(filter = group == 1,label=name),size=3) +
    geom_node_text(aes(filter = group == 2,label=name),size=2, repel = T) +
    theme_void() + theme(legend.position = 'none')
ggsave(filename = '/home/yzj/JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/Fig5A_TF_CSC_down.pdf',plot_CSC_down,
       height = 10, width = 10, units = 'cm')


#####################
####### down
#####################
file_CSC_up <- paste0(path.TF.net, 'igraph_CSC_up.RDS')
igraph_up_2 <- readRDS(file_CSC_up)

sel.BP_SCS_up <- c('cell cycle arrest',
                   'negative regulation of stem cell differentiation',
                   'response to oxidative stress',
                   'response to unfolded protein',
                   'intrinsic apoptotic signaling pathway',
                   'p38MAPK cascade',
                   'I-kappaB kinase/NF-kappaB signaling',
                   'cell death in response to oxidative stress',
                   'regulation of RNA stability',
                   'activation of innate immune response',
                   'cellular response to tumor necrosis factor', 
                   'negative regulation of cell death')

BP_SCS_up <- list.go.BP[['CSC_Microtia_increase']]
rownames(BP_SCS_up) <- BP_SCS_up$Description
# MF_SCS_up <- list.go.MF[['CSC_Microtia_increase']]
# rownames(MF_SCS_up) <- MF_SCS_up$Description
# df_GO <- rbind(BP_SCS_up[sel.BP_SCS_up, c('Description', 'geneID')],
#                MF_SCS_up[sel.MF_SCS_up, c('Description', 'geneID')])
df_GO <- BP_SCS_up[sel.BP_SCS_up, c('Description', 'geneID')]

node_type <- V(igraph_up_2)$name
list_gene_GO_up <- list()
for (gene in node_type) {
    vec_GO <- c()
    for (term in df_GO$Description) {
        gene.set <- strsplit(df_GO$geneID[df_GO$Description == term], split = '/')[[1]]
        if (gene %in% gene.set) {
            vec_GO <- c(vec_GO, term)
        }
    }
    list_gene_GO_up[[gene]] <- vec_GO
}


ggraph_CSC_up <- as_tbl_graph(igraph_up_2)

plot_CSC_up <- 
    ggraph(ggraph_CSC_up, layout = 'stress') +
    geom_edge_link(aes(edge_width=weight),color="gray",
                   arrow = arrow(length = unit(1.5, 'mm')),
                   end_cap = circle(1, 'mm'), 
                   start_cap = circle(2, 'mm')) +
    scale_edge_width(range=c(0.5,1)) +
    scale_edge_alpha(range=c(0.2,1)) + 
    scale_size_continuous(range = c(2,6)) +
    geom_node_point(aes(size = page_rank, fill = group),
                    shape=21, color = 'transparent') +
    # scale_color_manual(values = c('#0000CD', '#DC143C', 'dimgray')) +
    scale_fill_manual(values = c( 'firebrick3','#4169E1', '#BEBEBE')) +
    # scale_alpha_manual(values = c(0.8, 1, 0.6)) + 
    geom_node_text(aes(filter = group == 1,label=name),size=3) + 
    geom_node_text(aes(filter = group == 2,label=name), 
                   # nudge_x = 0.08, nudge_y = 0.005,
                   size=2, repel = T) + 
    theme_void() + theme(legend.position = 'none')
ggsave(filename = '/home/yzj/JingMA_NEW/res/compMicrotia/MicrotiavsNormal_inChildren/FIG/Fig5B_TF_CSC_up.pdf',plot_CSC_up,
       height = 10, width = 10, units = 'cm')






