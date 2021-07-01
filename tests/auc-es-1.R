library(splatter)
library(Seurat)
library(plotROC)
library(zoo)

auc_list <- list()
time0 <- Sys.time()
for (ITER in 1:10) {
  print(ITER)
  time1 <- Sys.time()
  print(time1 - time0)
  time0 <- time1
  #sim <- splatSimulate(nGenes = 1000,
  #                     batchCells = c(100, 200, 300, 400, 500, 600), batch.facLoc = 2., batch.facScale = .5,
  #                     group.prob = c(0.3, 0.7),
  #                     method = "groups", verbose = FALSE, seed=2020 + ITER)
  
  sim <- splatSimulate(nGenes = 1000,
                       batchCells = c(100, 200, 300, 400, 500, 600), batch.facLoc = 2., batch.facScale = .8,
                       group.prob = c(0.3, 0.7),
                       de.facLoc = 0.1,
                       bcv.common = 0.5,
                       method = "groups", verbose = FALSE, seed=2020 + ITER)
  
  gene_info = as.data.frame(rowData(sim))
  gene_info
  
  true_markers = rownames(gene_info)[(abs(gene_info$DEFacGroup1 - 1) > 1e-6) | (abs(gene_info$DEFacGroup2 - 1) > 1e-6)]
  length(true_markers)
  
  meta.data <- as.data.frame(colData(sim))
  batch_names <- unique(meta.data$Batch)
  for (i in batch_names){
    meta.data[[i]] <- as.numeric(meta.data$Batch == i)
  }
  meta.data
  
  obj <- CreateSeuratObject(counts(sim), meta.data = meta.data)
  obj <- NormalizeData(obj)
  Idents(obj) = obj$Group
  
  
  gene_auc = function(true_gene, pred_gene){
    fpr = rep(0, length(pred_gene))
    tpr = rep(0, length(pred_gene))
    
    n_good = 0
    n_bad = 0
    
    n_p = length(true_gene)
    n_n = length(pred_gene) - n_p
    for (i in 1:length(pred_gene)){
      if (pred_gene[i] %in% true_gene) {
        n_good = n_good + 1
      } else {
        n_bad = n_bad + 1
      }
      
      tpr[i] = n_good / n_p
      fpr[i] = n_bad / n_n
    }
    sum(diff(c(0, fpr)) * rollmean(c(0, tpr), 2))
  }
  
  tests = c('MAST')
  names(tests) = c('MAST')
  markers_mast = list()
  markers_mast$MAST <-  FindMarkers(obj, test.use = 'MAST', 
                                    logfc.threshold = -Inf, min.pct = -Inf, min.cells.feature = 0,
                                    ident.1 = 'Group1', ident.2 = 'Group2')
  
  markers_mast$MAST_with_covariate <- FindMarkers(obj, test.use = 'MAST', latent.vars = batch_names,
                                                  logfc.threshold = -Inf, min.pct = -Inf, min.cells.feature = 0,
                                                  ident.1 = 'Group1', ident.2 = 'Group2')
  
  tests = c('wilcox', 'bimod', 'roc', 't', 'negbinom', 'poisson', 'LR', 'DESeq2')
  names(tests) = c('Wilcoxon', 'Bimod', 'ROC_Analysis', 't-test', 'Negative_Binomial', 'Poisson', 'LR', 'DESeq2')
  
  markers_orig <- lapply(tests,
                         function(test.use) {
                           FindMarkers(obj, test.use = test.use, slot="counts",
                                       logfc.threshold = -Inf, min.pct = -Inf, min.cells.feature = 0,
                                       ident.1 = 'Group1', ident.2 = 'Group2')
                         })
  
  # Note: MAST only accept numeric latent var
  tests = c('negbinom', 'poisson', 'LR')
  names(tests) = c('Negative_Binomial_with_covariate', 'Poisson_with_covariate', 'LR_with_covariate')
  
  markers_covar <- lapply(tests,
                          function(test.use) {
                            FindMarkers(obj, test.use = test.use, latent.vars = batch_names, slot="counts",
                                        logfc.threshold = -Inf, min.pct = -Inf, min.cells.feature = 0,
                                        ident.1 = 'Group1', ident.2 = 'Group2')
                          })
  
  source('../R/stratified_test_seurat_rename.R')
  
  tests = c('VE')
  names(tests) = c('vanElteren')
  
  markers_ve <- lapply(tests,
                       function(test.use) {
                         FindMarkers2.Seurat(obj, test.use = test.use, latent.vars = 'Batch', slot="counts",
                                             logfc.threshold = -Inf, min.pct = -Inf, min.cells.feature = 0,
                                             ident.1 = 'Group1', ident.2 = 'Group2')
                       })
  
  tests = c('DESeq2')
  names(tests) = c('DESeq2 +')
  
  markers_deseq2 <- lapply(tests,
                           function(test.use) {
                             FindMarkers2.Seurat(obj, test.use = test.use, latent.vars = 'Batch', slot="counts",
                                                 logfc.threshold = -Inf, min.pct = -Inf, min.cells.feature = 0,
                                                 ident.1 = 'Group1', ident.2 = 'Group2')
                           })
  
  roc_auc = function(true_markers, markers, test) {
    df <- data.frame(D = as.numeric(rownames(markers) %in% true_markers), 
                     M = -log10(markers$p_val),  
                     Test = test)
    plt = ggplot(df, aes(m = M, d = D)) + geom_roc() + style_roc()
    list(df = df, auc = calc_auc(plt)$AUC)
  }
  
  roc_df = function(true_markers, markers, test) {
    if ('p_val' %in% colnames(markers)) {
      df <- data.frame(D = as.numeric(rownames(markers) %in% true_markers), 
                       M = -log10(markers$p_val),  
                       Test = test)
    }
    else {
      df <- data.frame(D = as.numeric(rownames(markers) %in% true_markers), 
                       M = markers$power,  
                       Test = test)
    }
    df
  }
  
  ##########################################
  # SigEMD
  ##########################################
  library(aod)
  library(arm)
  library(fdrtool)
  library(lars)
  source("SigEMD-master/FunImpute.R")
  source("SigEMD-master/SigEMDHur.R")
  source("SigEMD-master/SigEMDnonHur.R")
  source("SigEMD-master/plot_sig.R")
  data <- as.matrix(obj@assays$RNA@data)
  condition <- obj$Group

  databinary <- databin(data)
  Hur_gene<- idfyImpgene(data,databinary,condition)
  #genes_use<- idfyUsegene(data,databinary,condition,ratio = 0.8) # ratio is default as 0.8, but it can be set by users to obtain appropriate number of genes_use.
  #{sink("/dev/null"); datimp <- FunImpute(object = data, genes_use = (genes_use), genes_fit = (Hur_gene),dcorgene = NULL); sink();}
  #data<-datimp$alldat
  data<-data[rowSums(abs(data)) > 0, ]
  results <- calculate_single(data =  data,condition =  condition,Hur_gene = Hur_gene, binSize=0.2,nperm=200)

  markers_sigemd <- list(SigEMD=results$emdall)
  colnames(markers_sigemd$SigEMD) <- c('emd', 'p_val', 'p_val_adj')

  markers_sigemd$SigEMD <- as.data.frame(markers_sigemd$SigEMD)
  markers_sigemd$SigEMD <- markers_sigemd$SigEMD[order(markers_sigemd$SigEMD$p_val), ]
  
  ##########################################
  # DESingle
  ##########################################
  library(DEsingle)
  results <- DEsingle(counts = obj@assays$RNA@counts, group = obj$Group, parallel = TRUE)

  markers_desingle <- results[c('pvalue', 'pvalue.adj.FDR')]
  colnames(markers_desingle) <- c('p_val', 'p_val_adj')
  markers_desingle = list(DEsingle=markers_desingle)
  
  
  ##########################################
  # Summarize
  ##########################################
  
  all_markers <- c(markers_orig, markers_mast, markers_covar, markers_ve[-2], 
                   markers_sigemd, markers_desingle, markers_deseq2)
  all_roc_dfs <- list()
  
  for (i in names(all_markers)){
    all_roc_dfs[[i]] <- roc_df(true_markers, all_markers[[i]], i)
  }
  
  
  all_auc = rep(0, length(all_markers))
  names(all_auc) = names(all_markers)
  for (i in names(all_markers)){
    all_auc[i] = gene_auc(true_markers, rownames(all_markers[[i]]))
  }
  
  auc_list[[ITER]] <- all_auc
}

saveRDS(object = auc_list, file = "auc-es-0.RDS")
