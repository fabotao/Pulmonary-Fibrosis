library(Seurat)
library(edgeR)
library(ggplot2)
library(ggpubr)

# --- Load Input Datasets ------------------------------------------------------
# All data files should be placed in the appropriate subfolders.
# Example folder structure:
#   ./processed_data/           - for main input files
#   ./output/         - for saving plots and results


## Compute cell type fraction in conditions
sce.all <- readRDS('processed_data/sce_all_celltype.RDS')

group.map <- data.frame(sample=c('C56','C58','C60','C62','C63','C64','C66','C67','C69','C70'),
                        stim=c('CIP','BLM','CIP','BLM','BLM','PBS','PBS','CIP','CIP','BLM'),
                        time=c('D28','D28','D10','D28','D10','D10','D10','D10','D28','D10'))
sce.all$stim <- group.map$stim[match(sce.all$orig.ident, group.map$sample)]
sce.all$time <- group.map$time[match(sce.all$orig.ident, group.map$sample)]


count.bulk <- cbind(rowSums(sce.all@assays$RNA@counts[,sce.all$stim=='PBS']),
                    rowSums(sce.all@assays$RNA@counts[,sce.all$stim=='BLM' & sce.all$time=='D10']),
                    rowSums(sce.all@assays$RNA@counts[,sce.all$stim=='BLM' & sce.all$time=='D28']))
colnames(count.bulk) <- c('PBS','BLM10','BLM28')


## DE analysis and volcano plot
DE.analysis <- function(id.neg, id.pos, raw.count){

  count <- raw.count

  countData <- count[,c(which(colnames(count) %in% id.neg), which(colnames(count) %in% id.pos))]
  condition <- c(rep('control', length(id.neg)), rep('pos', length(id.pos)))
  #colData <- data.frame(row.names=colnames(countData), condition)

  library(edgeR)

  exprSet <- DGEList(counts=countData, group = condition)
  keep <- rowSums(cpm(exprSet)>1) >= 1
  exprSet <- exprSet[keep, , keep.lib.sizes=FALSE]

  exprSet$samples$lib.size <- colSums(exprSet$counts)
  exprSet <- calcNormFactors(object = exprSet)
  exprSet <- estimateDisp(y = exprSet)

  bcv = 0.1#设置bcv为0.1
  et <- exactTest(exprSet, dispersion=bcv^2)
  DEG_edgeR=as.data.frame(topTags(et, n = nrow(exprSet$counts)))
  head(DEG_edgeR)

  return(DEG_edgeR)
}


DE.10 <- DE.analysis(id.neg = 'PBS', id.pos = 'BLM10', raw.count = count.bulk)
DE.28 <- DE.analysis(id.neg = 'PBS', id.pos = 'BLM28', raw.count = count.bulk)

#saveRDS(DE.10, file = 'output/PseudoBulk_DE_BLM10_vs_PBS.RDS')
#saveRDS(DE.28, file = 'output/PseudoBulk_DE_BLM28_vs_PBS.RDS')

DE.10.bulk <- readRDS('output/PseudoBulk_DE_BLM10_vs_PBS.RDS')
DE.28.bulk <- readRDS('output/PseudoBulk_DE_BLM28_vs_PBS.RDS')



type.diff <- function(obj, type, time){
  obj$group <- ifelse(obj$time==time & obj$cell_type==type & obj$stim=='BLM', 'BLM',
                      ifelse(obj$time=='D10' & obj$cell_type==type & obj$stim=='PBS', 'PBS', 'Oth'))
  Idents(obj) <- 'group'

  mk.type <- FindMarkers(obj, ident.1 = 'BLM', ident.2 = 'PBS')

  mk.type

}


## Exclude cell types with few cells
exclude.type <- c('Adipocytes', 'Basophils')
uni.types <- unique(sce.all$cell_type)


D10.DE.type <- list()
D28.DE.type <- list()
for(type in uni.types){
  if(!type %in% exclude.type) {
    print(type)
    D10.DE.type[[type]] <- type.diff(obj = sce.all, type=type, time='D10')
    D28.DE.type[[type]] <- type.diff(obj = sce.all, type=type, time='D28')
  }
}



prop.PBS <- prop.table(table(sce.all$cell_type[sce.all$time=='D10' & sce.all$stim=='PBS']))
prop.BLM.D10 <- prop.table(table(sce.all$cell_type[sce.all$time=='D10' & sce.all$stim=='BLM']))
prop.BLM.D28 <- prop.table(table(sce.all$cell_type[sce.all$time=='D28' & sce.all$stim=='BLM']))

all(names(prop.PBS)==names(prop.BLM.D10))
all(names(prop.PBS)==names(prop.BLM.D28))

Compute.score <- function(DE.type, DE.bulk, prop.BLM, prop.PBS){

  get.type.score <- function(type){
    DE.type.res <- DE.type[[type]]
    DE.type.res$gene <- rownames(DE.type.res)
    bulk.up.g <- rownames(DE.bulk)[order(DE.bulk$logFC, decreasing = T)[1:100]]
    bulk.down.g <- rownames(DE.bulk)[order(DE.bulk$logFC, decreasing = F)[1:100]]

    DE.type.up <- DE.type.res[order(DE.type.res$avg_log2FC, decreasing = T)[1:100],]
    DE.type.down <- DE.type.res[order(DE.type.res$avg_log2FC, decreasing = F)[1:100],]

    DE.type.up.cross <- DE.type.up[DE.type.up$gene %in% bulk.up.g,,drop=F]
    DE.type.down.cross <- DE.type.down[DE.type.down$gene %in% bulk.down.g,,drop=F]

    prop.type <- abs((prop.BLM[type]-prop.PBS[type]))

    get.FC.score <- function(DE.type.cross){
      if(dim(DE.type.cross)[1]==0){
        return(0)
      }else{
        return(sqrt(abs(DE.type.cross$avg_log2FC) * prop.type))
      }
    }

    FC.score <- c(get.FC.score(DE.type.up.cross), get.FC.score(DE.type.down.cross))
    return(mean(FC.score))
  }

  res <- unlist(lapply(names(DE.type), function(x){get.type.score(x)}))
  names(res) <- names(DE.type)
  res
}

BLM10.score <- Compute.score(DE.type = D10.DE.type, DE.bulk = DE.10.bulk, prop.BLM = prop.BLM.D10, prop.PBS = prop.PBS)
BLM28.score <- Compute.score(DE.type = D28.DE.type, DE.bulk = DE.28.bulk, prop.BLM = prop.BLM.D28, prop.PBS = prop.PBS)


bar.plot <- function(score, title){
  df <- data.frame(type=names(score),
                   score=score)
  df <- df[order(df$score, decreasing = T),]
  df$type <- factor(df$type, levels=df$type, ordered=T)
  p <- ggplot(data=df) + geom_bar(aes(x=type, y=score), stat = 'identity') + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust=1, colour = 'black'),
          axis.text.y = element_text(colour = 'black')) +
    labs(y='Contribution', title=title)
  p
}

p.BLM.D10 <- bar.plot(BLM10.score, title='Day10')
p.BLM.D28 <- bar.plot(BLM28.score, title='Day28')

p.contribution <- ggarrange(p.BLM.D10, p.BLM.D28,  nrow=1, align='h')
ggsave(p.contribution, filename = 'output/celltype_contribution_BLM.pdf', width=10, height = 4)


score.df <- data.frame(type=names(BLM10.score),
                       score_D10=BLM10.score,
                       score_D28=BLM28.score)
#write.csv(score.df, file = 'output/cell_type_contribution_score_BLM.csv', row.names = F)
