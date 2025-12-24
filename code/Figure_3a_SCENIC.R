library(Seurat)
library(SCENIC)
library(doParallel)
library(SCopeLoomR)
library(magrittr)
library(ggplot2)
library(dplyr)
library(org.Mm.eg.db)
library(ggrepel)
library(patchwork)
library(stringi)


# --- Load Input Datasets ------------------------------------------------------
# All data files should be placed in the appropriate subfolders.
# Example folder structure:
#   ./processed_data/           - for main input files
#   ./output/         - for saving plots and results



Neu.obj <- readRDS('processed_data/Neu_obj.RDS')

neu.mk1 <- c('Cd34','Sox4','Rpl12','Elane','Cebpa','Prtn3','Mpo','Chil3','Tuba1b','Fcnb','Ltf','Ngp','Camp','Retnlg','Cxcl2','Mmp8','Ccl6','Gm5483','Stfa2l1','Isg15','Rsad2','Ifit3','Fgl2','Gm2a','Gngt2')
DotPlot(Neu.obj, features = neu.mk1, cluster.idents = T) + coord_flip()

Neu.obj$subtype_mod <- c('Aged','Mature','Mature','Immature','Mature','Mature','Mature','Immature','Aged','Pro/pre','Mature','Mature',
                         'Mature','Mature','Aged','Aged','Pro/pre','Aged','Mature','Pro/pre','Pro/pre','Pro/pre')[match(Neu.obj$seurat_clusters, c(0:21))]




## SCENIC
library(SCENIC)

exprMat <- as.matrix(Neu.obj@assays$RNA@data)
dim(exprMat)
exprMat[1:4,1:4]
cellInfo <-  Neu.obj@meta.data[,c('time','stim','subtype_mod')]
cellInfo$subtype <- paste0(cellInfo$stim,'_',cellInfo$subtype_mod)
cellInfo <- cellInfo[,c('time', 'stim', 'subtype')]
head(cellInfo)

### Initialize settings
# 保证cisTarget_databases 文件夹下面有下载好2个1G的文件
scenicOptions <- initializeScenic(org="mgi",
                                  dbDir="/home/rstudio/cell_type/SCENIC_db",
                                  nCores=30)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings`

## Save output
export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

rm(list = ls())
library(Seurat)
library(SCENIC)
library(doParallel)

scenicOptions=readRDS(file="int/scenicOptions.Rds")

### Exploring output
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
as.data.frame(sort(table(motifEnrichment_selfMotifs_wGenes$highlightedTFs),decreasing = T))

## Take IRF7 as example
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Irf7"]
viewMotifs(tableSubset)

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="IRF7" & highConfAnnot==TRUE]
viewMotifs(tableSubset)


rm(list = ls())
gc()
library(Seurat)
library(SCENIC)
library(doParallel)
library(SCopeLoomR)
scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
#loom <- open_loom(scenicLoomPath)
# Read information from loom file:
#regulons_incidMat <- get_regulons(loom)
#regulons <- regulonsToGeneLists(regulons_incidMat)
#regulonsAUC <- get_regulons_AUC(loom)
#regulonsAucThresholds <- get_regulon_thresholds(loom)
#embeddings <- get_embeddings(loom)



tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

# Show TF expression:
AUCell::AUCell_plotTSNE(tSNE_scenic$Y,
                        exprMat,
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Nfkb2", "Irf7", "Irf9","Nfkb1")],], plots="Expression")
pdf(file="output/AUCell_plotTSNE.pdf")
par(mfrow=c(4,5))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()


regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$subtype),
                                     function(cells) rowMeans(AUCell:::getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))


pdf(file="output/regulonActivity_by_SubType_Scaled.pdf")
col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#003366", "white", "#D32F2f"))
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",col=col_fun)
dev.off()


regulonActivity_byCellType_Scaled_Subset <- t(scale(t(regulonActivity_byCellType[c('Nfkb1 (66g)',
                                                                                   'Smad3 (19g)',
                                                                                   'Stat1 (31g)',
                                                                                   'Irf9 (31g)',
                                                                                   'Junb (16g)',
                                                                                   'Myc (21g)',
                                                                                   'Fos_extended (22g)',
                                                                                   'Ets1_extended (14g)',
                                                                                   'Nfe2_extended (21g)',
                                                                                   'Klf2_extended (23g)',
                                                                                   'Cebpb_extended (414g)'),
                                                                                 c('PBS_Pro/pre', 'BLM_Pro/pre', 'CIP_Pro/pre',
                                                                                   'PBS_Immature', 'BLM_Immature', 'CIP_Immature',
                                                                                   'PBS_Mature', 'BLM_Mature', 'CIP_Mature',
                                                                                   'PBS_Aged', 'BLM_Aged', 'CIP_Aged')]), center = T, scale=T))

pdf(file="output/regulonActivity_by_SubType_Scaled_Subset.pdf")
col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#003366", "white", "#D32F2f"))
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled_Subset, name="Regulon activity",col=col_fun, cluster_columns = F)
dev.off()




topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
head(topRegulators)


minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$subtype),
                                               function(cells) {
                                                 rowMeans(binaryRegulonActivity[,cells, drop=FALSE])
                                               })
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]

pdf(file="output/binaryActPerc_subset.pdf")
ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Regulon activity (%)", col = c("white","pink","red"))
dev.off()


topRegulators <- reshape2::melt(regulonActivity_byCellType_Binarized)
colnames(topRegulators) <- c("Regulon", "Subtype", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>minPerc),]
head(topRegulators)


rss <- calcRSS(AUC=AUCell:::getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "subtype"])
pdf(file="output/rss.pdf")
plotRSS_oneSet(rss, setName = "Aged")
dev.off()

