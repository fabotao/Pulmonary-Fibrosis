setwd('/home/rstudio/Projects/Lung_PD1/')

library(Seurat)
library(magrittr)
library(SingleR)
library(ggplot2)
library(dplyr)
library(org.Mm.eg.db)
library(ggpubr)
library(CellChat)


## Cell-cell communication


sce.PBS.CIP.subtype.mod <- readRDS('processed_data/sce_PBS_CIP_subtype_mod.RDS')
sce.PBS.CIP.subtype.mod$subtype_mod[sce.PBS.CIP.subtype.mod$subtype_mod=='unknown' & sce.PBS.CIP.subtype.mod$cell_type=='Monocytes'] <- 'Fcn1-high Mono'
sce.PBS.CIP.subtype.mod$subtype_mod[sce.PBS.CIP.subtype.mod$subtype_mod=='unknown' & sce.PBS.CIP.subtype.mod$cell_type=='T cells'] <- 'Trm'



cellchat <- createCellChat(object = sce.PBS.CIP.subtype.mod,
                           meta = sce.PBS.CIP.subtype.mod@meta.data,
                           group.by = "subtype_mod")
cellchat

CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
# set the used database in the object
cellchat@DB <- CellChatDB.use


# This step is necessary even if using the whole database
cellchat <- subsetData(cellchat)
# do parallel ，根据配置设置
future::plan("multisession", workers = 1)

#识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
#识别过表达配体受体对
cellchat <- identifyOverExpressedInteractions(cellchat)

#project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat@data.project[1:4,1:4]


cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)


#all the inferred cell-cell communications at the level of ligands/receptors
df.net <- subsetCommunication(cellchat)

#access the the inferred communications at the level of signaling pathways
df.net1 <- subsetCommunication(cellchat,slot.name = "netP")

#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
levels(cellchat@idents)
df.net2 <- subsetCommunication(cellchat, sources.use = c("Epi"), targets.use = c("Fibroblast" ,"T"))

#gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
df.net3 <- subsetCommunication(cellchat, signaling = c("CCL", "TGFb"))

#计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞类型之间通信结果
cellchat <- aggregateNet(cellchat)

#saveRDS(cellchat, file = 'processed_data/cellchat.RDS')


# Read RDS
cellchat <- readRDS('processed_data/cellchat.RDS')

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#指定受体-配体细胞类型
netVisual_bubble(cellchat, sources.use = c(1:10), targets.use = c(1:10),remove.isolate = FALSE)


par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")




source('/home/rstudio/cell_type/HCA_BM/netvisual_heatmap_mod.R', encoding = 'utf-8')


valid.paths <- c("TGFb","CCL","CXCL","IL2","CSF","IFN-I","TNF","VISFATIN","COMPLEMENT","ANNEXIN",
                 "GALECTIN")



# Focus on the Fibro,Neutrophil
idents <- c('CAF','Fibroblast','Lipofibroblast','Myofibroblast','Pericytes','Aged','Mature','Immature','Pro/pre')
p.net.circle.1 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize[idents],
                                   sources.use = idents,
                                   targets.use = idents,
                                   remove.isolate = T,
                                   weight.scale = T, label.edge= F, title.name = "Number of interactions")
p.net.circle.2 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize[idents],
                                   sources.use = idents,
                                   targets.use = idents,
                                   remove.isolate = T,
                                   weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
p.net.circle <- ggarrange(p.net.circle.1, p.net.circle.2, ncol=2)


pairLR.use <- extractEnrichedLR(cellchat, signaling = valid.paths)
p.pairs.Neu.to.B <- netVisual_bubble(cellchat, targets.use = c('CAF','Fibroblast','Lipofibroblast','Myofibroblast','Pericytes'), remove.isolate = T, pairLR.use = pairLR.use,
                                     angle.x = 45,
                                     sources.use = c('Aged','Mature','Immature','Pro/pre')) +
  scale_color_gradient(low='white', high='red')

p.pairs.B.to.Neu <- netVisual_bubble(cellchat, targets.use = c('Aged','Mature','Immature','Pro/pre'), remove.isolate = T, pairLR.use = pairLR.use,
                                     angle.x = 45,
                                     sources.use = c('CAF','Fibroblast','Lipofibroblast','Myofibroblast','Pericytes')) +
  scale_color_gradient(low='white', high='red')

ggsave(p.pairs.Neu.to.B, filename = 'output/Neu_to_Fibro_pairLR.pdf', width = 8, height = 4)
ggsave(p.pairs.B.to.Neu, filename = 'output/Fibro_to_Neu_pairLR.pdf', width = 8, height = 6)




