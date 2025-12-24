library(Seurat)
library(ggplot2)
library(ggpubr)
library(org.Mm.eg.db)

# --- Load Input Datasets ------------------------------------------------------
# All data files should be placed in the appropriate subfolders.
# Example folder structure:
#   ./processed_data/           - for main input files
#   ./output/         - for saving plots and results

sce = readRDS('processed_data/sce_all_celltype.RDS')

group.map <- data.frame(sample=c('C56','C58','C60','C62','C63','C64','C66','C67','C69','C70'),
                        stim=c('CIP','BLM','CIP','BLM','BLM','PBS','PBS','CIP','CIP','BLM'),
                        time=c('D28','D28','D10','D28','D10','D10','D10','D10','D28','D10'))
sce$stim <- group.map$stim[match(sce$orig.ident, group.map$sample)]
sce$time <- group.map$time[match(sce$orig.ident, group.map$sample)]



get.res <- function(stim1, time1){

  obj <- subset(sce, subset=(stim == stim1 & time == time1) | stim == 'PBS')
  obj$stim_type = paste0(obj$stim, '_', obj$cell_type)
  Idents(obj) = 'stim_type'

  type.uniq = sort(setdiff(unique(obj$cell_type), c('Basophils','Adipocytes')))
  DE.list <- lapply(type.uniq, function(x){
    mk <- FindMarkers(obj, ident.1 = paste0(stim1,'_',x), ident.2 = paste0('PBS_',x))
    return(mk)
  })
  names(DE.list) <- type.uniq

  return(DE.list)
}


out.BLM.D10 = get.res(stim1 = 'BLM', time1 = 'D10')
out.BLM.D28 = get.res(stim1 = 'BLM', time1 = 'D28')
out.CIP.D10 = get.res(stim1 = 'CIP', time1 = 'D10')
out.CIP.D28 = get.res(stim1 = 'CIP', time1 = 'D28')


gsea.ana <- function(DE.res, category){

  m_t2g <- msigdbr::msigdbr(species = "Mus musculus", category = category) %>%
    dplyr::select(gs_name, entrez_gene)

  genename <- rownames(DE.res)
  gene_map <- AnnotationDbi::select(org.Mm.eg.db, keys=genename, keytype="SYMBOL", columns=c("SYMBOL","ENTREZID"))
  colnames(gene_map)[1]<-"Gene"

  genelist_input <- data.frame(Gene=rownames(DE.res), logFC=DE.res$avg_log2FC)

  aaa<-inner_join(gene_map,genelist_input,by = "Gene")
  aaa<-aaa[,-1]
  aaa<-na.omit(aaa)
  aaa$logFC<-sort(aaa$logFC,decreasing = T)

  geneList = aaa[,1]
  names(geneList) = as.character(aaa[,1])
  geneList
  if(length(geneList) >= 1){
    #Go_gseresult <- GSEA(geneList, TERM2GENE=m_t2g)

    Go_gseresult <- clusterProfiler::enricher(geneList, TERM2GENE=m_t2g)
    #Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
    #KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
    #gseaplot2(Go_gseresult,1:5,pvalue_table = TRUE) #输出第1个结果
    go.res <- Go_gseresult@result
    go.res <- go.res[go.res$p.adjust < 0.01,]
    #return(go.res)
    return(summary(Go_gseresult))

  }else{
    return(NA)
  }

}





enrich.fun <- function(DE.list, up=T){
  types = names(DE.list)
  enrich.list <- lapply(DE.list, function(x){
    if(up){
      DE.res <- x[(x$avg_log2FC > 0.25 & x$p_val_adj < 0.01),]
    }else{
      DE.res <- x[(x$avg_log2FC < -0.25 & x$p_val_adj < 0.01),]
    }
    enrich.res.x = gsea.ana(DE.res = DE.res, category = 'C5')
    return(enrich.res.x)
  })
  names(enrich.list) = types
  return(enrich.list)
}


go.res.BLM.D10.up = enrich.fun(out.BLM.D10)
go.res.BLM.D28.up = enrich.fun(out.BLM.D28)
go.res.CIP.D10.up = enrich.fun(out.CIP.D10)
go.res.CIP.D28.up = enrich.fun(out.CIP.D28)


go.res.BLM.D10.down = enrich.fun(out.BLM.D10, up=F)
go.res.BLM.D28.down = enrich.fun(out.BLM.D28, up=F)
go.res.CIP.D10.down = enrich.fun(out.CIP.D10, up=F)
go.res.CIP.D28.down = enrich.fun(out.CIP.D28, up=F)


output.pathway <- function(go.res){

  path.rbind = do.call(rbind, lapply(names(go.res), function(x){
    path.x = go.res[[x]]
    if(!is.data.frame(path.x)){
      path.x = blank.df
    }else if(dim(path.x)[1]==0){
      path.x = blank.df
    }

    path.x$type = x
    return(path.x)
  }))

  return(path.rbind)
}

write.csv(output.pathway(go.res.BLM.D10.up), file = 'output/Go_res_BLM_D10_up.csv')
write.csv(output.pathway(go.res.BLM.D28.up), file = 'output/Go_res_BLM_D28_up.csv')
write.csv(output.pathway(go.res.CIP.D10.up), file = 'output/Go_res_CIP_D10_up.csv')
write.csv(output.pathway(go.res.CIP.D28.up), file = 'output/Go_res_CIP_D28_up.csv')
write.csv(output.pathway(go.res.BLM.D10.down), file = 'output/Go_res_BLM_D10_down.csv')
write.csv(output.pathway(go.res.BLM.D28.down), file = 'output/Go_res_BLM_D28_down.csv')
write.csv(output.pathway(go.res.CIP.D10.down), file = 'output/Go_res_CIP_D10_down.csv')
write.csv(output.pathway(go.res.CIP.D28.down), file = 'output/Go_res_CIP_D28_down.csv')



blank.df <- data.frame(ID=NA,
                       Description=NA,
                       GeneRatio=NA,
                       BgRatio=NA,
                       pvalue=NA,
                       p.adjust=NA,
                       qvalue=NA,
                       geneID=NA,
                       Count=NA)

get.top.path <- function(BLM.enrich, CIP.enrich){

  types = names(BLM.enrich)

  path.res <- do.call(rbind, lapply(types, function(x){
    print(x)
    BLM.path = BLM.enrich[[x]]

    if(!is.data.frame(BLM.path)){
      BLM.path.top = blank.df
    }else if(dim(BLM.path)[1] == 0){
      BLM.path.top = blank.df
    }else{
      BLM.path.top = BLM.path[1:min(c(10, dim(BLM.path)[1])),]
    }
    BLM.path.top$stim = 'BLM'
    BLM.path.top$type = x

    CIP.path = CIP.enrich[[x]]
    if(!is.data.frame(CIP.path)){
      CIP.path.top = blank.df
    }else if(dim(CIP.path)[1] ==0){
      CIP.path.top = blank.df
    }else{
      CIP.path.top = CIP.path[1:min(c(10, dim(CIP.path)[1])),]
    }

    CIP.path.top$stim = 'CIP'
    CIP.path.top$type = x

    BLM.CIP.dif.path = setdiff(BLM.path.top$ID, CIP.path.top$ID)
    CIP.BLM.dif.path = setdiff(CIP.path.top$ID, BLM.path.top$ID)


    if(length(BLM.CIP.dif.path)==0){
      BLM.CIP.dif = blank.df
    }else if(length(BLM.CIP.dif.path) > 0 & !is.na(BLM.CIP.dif.path)){
      BLM.CIP.dif = BLM.path[BLM.path$ID %in% BLM.CIP.dif.path,]
    }else{
      BLM.CIP.dif = blank.df
    }
    BLM.CIP.dif$stim = 'BLM'
    BLM.CIP.dif$type = x

    if(length(CIP.BLM.dif.path)==0){
      CIP.BLM.dif = blank.df
    }else if(length(CIP.BLM.dif.path) > 0 & !is.na(CIP.BLM.dif.path)){
      CIP.BLM.dif = CIP.path[CIP.path$ID %in% CIP.BLM.dif.path,]
    }else{
      CIP.BLM.dif = blank.df
    }
    CIP.BLM.dif$stim = 'CIP'
    CIP.BLM.dif$type = x

    return(rbind(BLM.path.top, BLM.CIP.dif, CIP.path.top, CIP.BLM.dif))
  }))

  return(path.res)
}


D10.path.up = get.top.path(BLM.enrich=go.res.BLM.D10.up, CIP.enrich=go.res.CIP.D10.up)
D28.path.up = get.top.path(BLM.enrich=go.res.BLM.D28.up, CIP.enrich=go.res.CIP.D28.up)
D10.path.down = get.top.path(BLM.enrich=go.res.BLM.D10.down, CIP.enrich=go.res.CIP.D10.down)
D28.path.down = get.top.path(BLM.enrich=go.res.BLM.D28.down, CIP.enrich=go.res.CIP.D28.down)



plot.path <- function(path.info, col){
  types = sort(unique(path.info$type))
  plot.list <- lapply(types, function(x){
    path.type = path.info[path.info$type == x,]
    path.type <- path.type[!is.na(path.type$ID),]
    path.type$Ratio = as.numeric(sapply(strsplit(path.type$GeneRatio, '/', fixed = T), '[', 1))/as.numeric(sapply(strsplit(path.type$GeneRatio, '/', fixed = T), '[', 2))
    path.type$log10qvalue = -log10(path.type$qvalue)
    path.type$Pathway <- gsub('_', ' ', path.type$ID, fixed = T)
    path.type$Pathway <- stringr::str_to_sentence(path.type$Pathway)
    p <- ggplot(data = path.type) + geom_point(aes(x=factor(stim), y=factor(Pathway), size=Ratio, color=log10qvalue)) +
      scale_color_gradient(low='grey', high=col) + labs(x='Enrich ratio', y=NULL, title = x) +
      theme(axis.text = element_text(colour = 'black'))
    return(p)
  })

  p.comb <- ggpubr::ggarrange(plotlist = plot.list, ncol=2, nrow=6)
  return(p.comb)
}

ggsave(plot.path(D10.path.up, col='#D32F2F'), filename = 'output/D10_up_pathway.pdf', width = 15, height = 26)
ggsave(plot.path(D28.path.up, col='#D32F2F'), filename = 'output/D28_up_pathway.pdf', width = 15, height = 26)
ggsave(plot.path(D10.path.down, col='#003366'), filename = 'output/D10_down_pathway.pdf', width = 12, height = 26)
ggsave(plot.path(D28.path.down, col='#003366'), filename = 'output/D28_down_pathway.pdf', width = 14, height = 26)





up.shared.BLM.CIP = c('GOBP_ACTIVATION_OF_IMMUNE_RESPONSE',
                      'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_ANTIGEN',
                      'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_ANTIGEN',
                      'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN',
                      'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION',
                      'GOBP_BLOOD_VESSEL_ENDOTHELIAL_CELL_MIGRATION',
                      'GOBP_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
                      'GOBP_CELL_CHEMOTAXIS',
                      'GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS',
                      'GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY',
                      'GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_IMMUNE_RESPONSE',
                      'GOBP_DENDRITIC_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION',
                      'GOBP_GRANULOCYTE_CHEMOTAXIS',
                      'GOBP_GRANULOCYTE_MIGRATION',
                      'GOBP_IMMUNE_RESPONSE_REGULATING_CELL_SURFACE_RECEPTOR_SIGNALING_PATHWAY',
                      'GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY',
                      'GOBP_IMMUNOGLOBULIN_PRODUCTION_INVOLVED_IN_IMMUNOGLOBULIN_MEDIATED_IMMUNE_RESPONSE',
                      'GOBP_INTERLEUKIN_6_PRODUCTION',
                      'GOBP_INTERLEUKIN_8_PRODUCTION',
                      'GOBP_LEUKOCYTE_CHEMOTAXIS',
                      'GOBP_LEUKOCYTE_MEDIATED_IMMUNITY',
                      'GOBP_LEUKOCYTE_MIGRATION',
                      'GOBP_LEUKOCYTE_PROLIFERATION',
                      'GOBP_POSITIVE_REGULATION_OF_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION',
                      'GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE',
                      'GOBP_REGULATION_OF_LEUKOCYTE_DIFFERENTIATION',
                      'GOBP_REGULATION_OF_LEUKOCYTE_MIGRATION',
                      'GOBP_REGULATION_OF_LEUKOCYTE_PROLIFERATION',
                      'GOBP_REGULATION_OF_NEUTROPHIL_CHEMOTAXIS',
                      'GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR',
                      'GOBP_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY',
                      'GOBP_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION')


down.shared.BLM.CIP <- c('GOBP_ACTIVATION_OF_IMMUNE_RESPONSE',
                         'GOBP_ACUTE_INFLAMMATORY_RESPONSE',
                         'GOBP_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY',
                         'GOBP_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
                         'GOBP_CELL_CHEMOTAXIS',
                         'GOBP_CELLULAR_DEFENSE_RESPONSE',
                         'GOBP_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS',
                         'GOBP_CELLULAR_RESPONSE_TO_INTERFERON_GAMMA',
                         'GOBP_CELLULAR_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES',
                         'GOBP_CHAPERONE_COFACTOR_DEPENDENT_PROTEIN_REFOLDING',
                         'GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY',
                         'GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE',
                         'GOBP_ENDOTHELIAL_CELL_MIGRATION',
                         'GOBP_EPITHELIAL_CELL_PROLIFERATION',
                         'GOBP_GRANULOCYTE_ACTIVATION',
                         'GOBP_GRANULOCYTE_CHEMOTAXIS',
                         'GOBP_GRANULOCYTE_MIGRATION',
                         'GOBP_IMMUNE_RESPONSE_REGULATING_CELL_SURFACE_RECEPTOR_SIGNALING_PATHWAY',
                         'GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY',
                         'GOBP_INTERFERON_GAMMA_PRODUCTION',
                         'GOBP_INTERLEUKIN_6_PRODUCTION',
                         'GOBP_INTERLEUKIN_8_PRODUCTION',
                         'GOBP_LEUKOCYTE_CELL_CELL_ADHESION',
                         'GOBP_LEUKOCYTE_CHEMOTAXIS',
                         'GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY',
                         'GOBP_LEUKOCYTE_MEDIATED_IMMUNITY',
                         'GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE',
                         'GOBP_LEUKOCYTE_MIGRATION',
                         'GOBP_LEUKOCYTE_PROLIFERATION',
                         'GOBP_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
                         'GOBP_MACROPHAGE_ACTIVATION',
                         'GOBP_MYELOID_LEUKOCYTE_ACTIVATION',
                         'GOBP_MYELOID_LEUKOCYTE_MIGRATION')



D10.up.only.BLM.CIP = c('GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_PEPTIDE_ANTIGEN',
                        'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I',
                        'GOBP_DEVELOPMENTAL_CELL_GROWTH',
                        'GOBP_DEVELOPMENTAL_GROWTH_INVOLVED_IN_MORPHOGENESIS',
                        'GOBP_GRANULOCYTE_MACROPHAGE_COLONY_STIMULATING_FACTOR_PRODUCTION',
                        'GOBP_IMMUNOGLOBULIN_PRODUCTION',
                        'GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS',
                        'GOBP_INTERLEUKIN_17_PRODUCTION',
                        'GOBP_LEUKOCYTE_ADHESION_TO_VASCULAR_ENDOTHELIAL_CELL',
                        'GOBP_MEMBRANE_RAFT_ORGANIZATION',
                        'GOBP_NEGATIVE_REGULATION_OF_CELL_SUBSTRATE_ADHESION',
                        'GOBP_NEGATIVE_REGULATION_OF_CELL_SUBSTRATE_JUNCTION_ORGANIZATION',
                        'GOBP_NEGATIVE_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS',
                        'GOBP_NEGATIVE_REGULATION_OF_VASCULAR_PERMEABILITY',
                        'GOBP_POSITIVE_REGULATION_OF_ALPHA_BETA_T_CELL_PROLIFERATION',
                        'GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE',
                        'GOBP_POSITIVE_REGULATION_OF_TYPE_2_IMMUNE_RESPONSE',
                        'GOBP_POSITIVE_REGULATION_OF_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS',
                        'GOBP_POSITIVE_REGULATION_OF_VASCULAR_ASSOCIATED_SMOOTH_MUSCLE_CELL_MIGRATION',
                        'GOBP_POSITIVE_REGULATION_OF_VASCULAR_ASSOCIATED_SMOOTH_MUSCLE_CELL_PROLIFERATION',
                        'GOBP_POSITIVE_THYMIC_T_CELL_SELECTION',
                        'GOBP_REGULATION_OF_ALPHA_BETA_T_CELL_DIFFERENTIATION')


D28.up.only.BLM.CIP = c('GOBP_AEROBIC_RESPIRATION',
                        'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN',
                        'GOBP_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS',
                        'GOBP_INFLAMMATORY_RESPONSE_TO_WOUNDING',
                        'GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_DNA_DAMAGE_BY_P53_CLASS_MEDIATOR',
                        'GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_DNA_DAMAGE',
                        'GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_OXIDATIVE_STRESS',
                        'GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN',
                        'GOBP_MONOCYTE_CHEMOTAXIS',
                        'GOBP_NEGATIVE_REGULATION_OF_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS',
                        'GOBP_NEGATIVE_REGULATION_OF_OSSIFICATION',
                        'GOBP_NEGATIVE_REGULATION_OF_OXIDATIVE_STRESS_INDUCED_CELL_DEATH',
                        'GOBP_POSITIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY',
                        'GOBP_POSITIVE_REGULATION_OF_JNK_CASCADE',
                        'GOBP_POSITIVE_REGULATION_OF_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_PRODUCTION')


D10.down.only.BLM.CIP = c('GOBP_ACUTE_PHASE_RESPONSE',
                          'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION',
                          'GOBP_ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE',
                          'GOBP_CD4_POSITIVE_OR_CD8_POSITIVE_ALPHA_BETA_T_CELL_LINEAGE_COMMITMENT',
                          'GOBP_CHEMOKINE_PRODUCTION',
                          'GOBP_GLUTATHIONE_METABOLIC_PROCESS',
                          'GOBP_HUMORAL_IMMUNE_RESPONSE',
                          'GOBP_HYDROGEN_PEROXIDE_METABOLIC_PROCESS',
                          'GOBP_I_KAPPAB_KINASE_NF_KAPPAB_SIGNALING',
                          'GOBP_INTRACELLULAR_STEROL_TRANSPORT',
                          'GOBP_NEGATIVE_REGULATION_OF_FIBROBLAST_PROLIFERATION',
                          'GOBP_NEGATIVE_REGULATION_OF_INTERLEUKIN_1_BETA_PRODUCTION',
                          'GOBP_NEGATIVE_REGULATION_OF_LONG_TERM_SYNAPTIC_POTENTIATION',
                          'GOBP_NEGATIVE_REGULATION_OF_MYELOID_CELL_DIFFERENTIATION',
                          'GOBP_NEGATIVE_REGULATION_OF_NF_KAPPAB_TRANSCRIPTION_FACTOR_ACTIVITY',
                          'GOBP_NEGATIVE_REGULATION_OF_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION',
                          'GOBP_RESPONSE_TO_INTERFERON_BETA',
                          'GOBP_RESPONSE_TO_PEPTIDOGLYCAN')


D28.down.only.BLM.CIP <- c('GOBP_AMYLOID_PRECURSOR_PROTEIN_CATABOLIC_PROCESS',
                           'GOBP_AUTOPHAGIC_CELL_DEATH',
                           'GOBP_CALCIUM_ION_TRANSPORT_INTO_CYTOSOL',
                           'GOBP_CALCIUM_ION_TRANSPORT',
                           'GOBP_CDC42_PROTEIN_SIGNAL_TRANSDUCTION',
                           'GOBP_CELL_CELL_ADHESION_MEDIATED_BY_INTEGRIN',
                           'GOBP_CONNECTIVE_TISSUE_DEVELOPMENT',
                           'GOBP_ER_NUCLEUS_SIGNALING_PATHWAY',
                           'GOBP_ESTABLISHMENT_OF_LYMPHOCYTE_POLARITY',
                           'GOBP_GAMMA_DELTA_T_CELL_DIFFERENTIATION',
                           'GOBP_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS',
                           'GOBP_INNATE_IMMUNE_RESPONSE_ACTIVATING_SIGNAL_TRANSDUCTION',
                           'GOBP_NEGATIVE_REGULATION_OF_LEUKOCYTE_MEDIATED_IMMUNITY',
                           'GOBP_NEGATIVE_REGULATION_OF_LYMPHOCYTE_MEDIATED_IMMUNITY',
                           'GOBP_NEGATIVE_REGULATION_OF_LYMPHOCYTE_MIGRATION',
                           'GOBP_NEGATIVE_REGULATION_OF_MUSCLE_CELL_DIFFERENTIATION',
                           'GOBP_NEGATIVE_REGULATION_OF_PHAGOCYTOSIS',
                           'GOBP_NEGATIVE_REGULATION_OF_T_HELPER_CELL_DIFFERENTIATION')




go.res.BLM.D10.up.all <- read.csv('output/Go_res_BLM_D10_up.csv')
go.res.BLM.D28.up.all <- read.csv('output/Go_res_BLM_D28_up.csv')
go.res.CIP.D10.up.all <- read.csv('output/Go_res_CIP_D10_up.csv')
go.res.CIP.D28.up.all <- read.csv('output/Go_res_CIP_D28_up.csv')
go.res.BLM.D10.down.all <- read.csv('output/Go_res_BLM_D10_down.csv')
go.res.BLM.D28.down.all <- read.csv('output/Go_res_BLM_D28_down.csv')
go.res.CIP.D10.down.all <- read.csv('output/Go_res_CIP_D10_down.csv')
go.res.CIP.D28.down.all <- read.csv('output/Go_res_CIP_D28_down.csv')



visualize.path <- function(path.info, path.BLM, path.CIP, up=T){
  path.BLM$stim = 'BLM'
  path.CIP$stim = 'CIP'

  path.BLM.subset = path.BLM[path.BLM$ID %in% path.info,c('ID','p.adjust','type')]
  path.CIP.subset = path.CIP[path.CIP$ID %in% path.info,c('ID','p.adjust','type')]

  path.BLM.subset.df <- reshape2::dcast(data=path.BLM.subset, formula = ID ~ type, value.var = 'p.adjust')
  path.CIP.subset.df <- reshape2::dcast(data=path.CIP.subset, formula = ID ~ type, value.var = 'p.adjust')

  path.BLM.subset.mat <- as.matrix(path.BLM.subset.df[,2:dim(path.BLM.subset.df)[2]])
  rownames(path.BLM.subset.mat) = path.BLM.subset.df$ID

  path.CIP.subset.mat <- as.matrix(path.CIP.subset.df[,2:dim(path.CIP.subset.df)[2]])
  rownames(path.CIP.subset.mat) = path.CIP.subset.df$ID

  all(rownames(path.BLM.subset.mat)==rownames(path.CIP.subset.mat))

  common.path = intersect(rownames(path.BLM.subset.mat), rownames(path.CIP.subset.mat))

  path.comb <- cbind(path.BLM.subset.mat[common.path,], path.CIP.subset.mat[common.path,])

  if(up){
    high.val = 6
    col_fun = circlize::colorRamp2(c(0, high.val), c("white", "#D32F2f"))
  }else{
    high.val = 5
    col_fun = circlize::colorRamp2(c(0, high.val), c("white", "#003366"))
  }

  ht = (ComplexHeatmap::Heatmap(-log10(path.comb), cluster_rows = F, cluster_columns = F, row_names_side = 'left', column_split = c(rep('BLM', dim(path.BLM.subset.mat)[2]),
                                                                                                                                   rep('CIP', dim(path.CIP.subset.mat)[2])),
                               name = '-log10p', col = col_fun))
  return(ggplotify::as.ggplot(ht))
}


p1 = visualize.path(path.info = up.shared.BLM.CIP,
               path.BLM = go.res.BLM.D10.up.all,
               path.CIP = go.res.CIP.D10.up.all,
               up=T)
p2= visualize.path(path.info = up.shared.BLM.CIP,
               path.BLM = go.res.BLM.D28.up.all,
               path.CIP = go.res.CIP.D28.up.all,
               up=T)
p3 = visualize.path(path.info = down.shared.BLM.CIP,
               path.BLM = go.res.BLM.D10.down.all,
               path.CIP = go.res.CIP.D10.down.all,
               up=F)
p4 = visualize.path(path.info = down.shared.BLM.CIP,
               path.BLM = go.res.BLM.D28.down.all,
               path.CIP = go.res.CIP.D28.down.all,
               up=F)


p.comb.common = ggarrange(p1,p2,p3,p4, ncol=2, nrow=2)
ggsave(p.comb.common, filename = 'output/Pathway_by_celltype_common_heatmap.pdf', width = 12, height = 14)




p.d10.up.only = visualize.path(path.info = D10.up.only.BLM.CIP,
                               path.BLM = go.res.BLM.D10.up.all,
                               path.CIP = go.res.CIP.D10.up.all,
                               up=T)
p.d28.up.only= visualize.path(path.info = D28.up.only.BLM.CIP,
                              path.BLM = go.res.BLM.D28.up.all,
                              path.CIP = go.res.CIP.D28.up.all,
                              up=T)
p.d10.down.only = visualize.path(path.info = D10.down.only.BLM.CIP,
                                 path.BLM = go.res.BLM.D10.down.all,
                                 path.CIP = go.res.CIP.D10.down.all,
                                 up=F)
p.d28.down.only = visualize.path(path.info = D28.down.only.BLM.CIP,
                                 path.BLM = go.res.BLM.D28.down.all,
                                 path.CIP = go.res.CIP.D28.down.all,
                                 up=F)

p.comb.only = ggarrange(p.d10.up.only,
                        p.d28.up.only,
                        p.d10.down.only,
                        p.d28.down.only, ncol=2, nrow=2)
ggsave(p.comb.only, filename = 'output/Pathway_by_celltype_only_heatmap.pdf', width = 12, height = 14)



