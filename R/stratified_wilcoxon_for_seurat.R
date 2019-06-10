# Modified from Seurat:::FindMarkers.default
# Per GNU GPL v3.0, modification and distribution is allowed.

ve_test_core <- function (x, y, correct = TRUE, ...) 
{
  alternative = "two.sided"
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  if (!is.numeric(y)) 
    stop("'y' must be numeric")
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  if (length(x) < 1L) 
    return(list(U = 0, V = 0))
  #stop("not enough (finite) 'x' observations")
  if (length(y) < 1L) 
    return(list(U = 0, V = 0))
  #stop("not enough (finite) 'y' observations")
  
  METHOD <- "Wilcoxon rank sum test"
  r <- rank(c(x, y))
  n.x <- as.double(length(x))
  n.y <- as.double(length(y))
  
  STATISTIC <- c(W = sum(r[seq_along(x)]) - n.x * (n.x + 1)/2)
  TIES <- (length(r) != length(unique(r)))
  NTIES <- table(r)
  z <- STATISTIC - n.x * n.y/2
  SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) -  sum(NTIES^3 - NTIES)/((n.x + n.y) * (n.x + n.y - 1))))
  if (correct) {
    CORRECTION <- switch(alternative, two.sided = sign(z) * 0.5)
    METHOD <- paste(METHOD, "with continuity correction")
  }
  # need to do something else here
  list(U = z - CORRECTION, V = SIGMA)
}

ve_test <- function (x, label, strata, correct = TRUE) 
{
  Us = c()
  Vs = c()
  Ms = c()
  Ns = c()
  Ss = c()
  for (s in unique(strata))
  {
    sub_x <- x[strata == s]
    sub_label <- label$group[strata == s]
    #print(sub_label)
    temp <- ve_test_core(sub_x[sub_label == "Group1"], sub_x[sub_label == "Group2"])
    Us = c(Us, temp$U)
    Vs = c(Vs, temp$V)
    Ms = c(Ms, sum(sub_label == "Group1"))
    Ns = c(Ns, sum(sub_label == "Group2"))
    Ss = c(Ss, s)
  }
  #print(data.frame(Ss, Us, Vs, Ms, Ns))
  list(p.value = pchisq(sum(Us / (Ms + Ns + 1)) ^ 2 / sum((Vs / (Ms + Ns + 1)) ^ 2), df = 1, lower.tail = F))
  #list(p.value = 1 - pnorm(abs(sum(Us / (Ms + Ns + 1))) / sqrt(sum((Vs / (Ms + Ns + 1)) ^ 2))))
}

library(pbapply)
library(future.apply)

VeDETest <- function (data.use, cells.1, cells.2, verbose = TRUE, ...) 
{
  strata <- list(...)$strata
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  data.use <- data.use[, rownames(x = group.info), drop = FALSE]
  my.sapply <- ifelse(test = verbose && Seurat:::PlanThreads() == 1, yes = pbsapply, no = future_sapply)
  strata.info = strata[rownames(x = group.info)]
  
  p_val <- my.sapply(X = 1:nrow(x = data.use), FUN = function(x) {
    return(ve_test(data.use[x, ], group.info, strata.info)$p.value)
  })
  return(data.frame(p_val, row.names = rownames(x = data.use)))
}



NewFindMarkers.default <- function (object, cells.1 = NULL, cells.2 = NULL, features = NULL, 
                                    logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                                    min.diff.pct = -Inf, verbose = TRUE, only.pos = FALSE, max.cells.per.ident = Inf, 
                                    random.seed = 1, latent.vars = NULL, min.cells.feature = 3, 
                                    min.cells.group = 3, pseudocount.use = 1, ...) 
{
  features <- Seurat:::`%||%`(features, rownames(x = object))
  methods.noprefiliter <- c("DESeq2", "zingeR")
  if (test.use %in% methods.noprefiliter) {
    features <- rownames(object)
    min.diff.pct <- -Inf
    logfc.threshold <- 0
  }
  if (length(x = cells.1) == 0) {
    stop("Cell group 1 is empty - no cells with identity class ", 
         cells.1)
  }
  else if (length(x = cells.2) == 0) {
    stop("Cell group 2 is empty - no cells with identity class ", 
         cells.2)
    return(NULL)
  }
  else if (length(x = cells.1) < min.cells.group) {
    stop("Cell group 1 has fewer than ", min.cells.group, 
         " cells")
  }
  else if (length(x = cells.2) < min.cells.group) {
    stop("Cell group 2 has fewer than ", min.cells.group, 
         " cells")
  }
  else if (any(!cells.1 %in% colnames(x = object))) {
    bad.cells <- colnames(object)[which(!as.character(x = cells.1) %in% 
                                          colnames(object))]
    stop("The following cell names provided to cells.1 are not present: ", 
         paste(bad.cells, collapse = ", "))
  }
  else if (any(!cells.2 %in% colnames(x = object))) {
    bad.cells <- colnames(object)[which(!as.character(x = cells.2) %in% 
                                          colnames(object))]
    stop("The following cell names provided to cells.2 are not present: ", 
         paste(bad.cells, collapse = ", "))
  }
  thresh.min <- 0
  pct.1 <- round(x = Matrix::rowSums(x = object[features, cells.1, 
                                                drop = FALSE] > thresh.min)/length(x = cells.1), digits = 3)
  pct.2 <- round(x = Matrix::rowSums(x = object[features, cells.2,
                                                drop = FALSE] > thresh.min)/length(x = cells.2), digits = 3)
  data.alpha <- cbind(pct.1, pct.2)
  colnames(x = data.alpha) <- c("pct.1", "pct.2")
  alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
  names(x = alpha.min) <- rownames(x = data.alpha)
  features <- names(x = which(x = alpha.min > min.pct))
  if (length(x = features) == 0) {
    stop("No features pass min.pct threshold")
  }
  alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, 
                                  FUN = min)
  features <- names(x = which(x = alpha.min > min.pct & alpha.diff > 
                                min.diff.pct))
  if (length(x = features) == 0) {
    stop("No features pass min.diff.pct threshold")
  }
  data.1 <- apply(X = object[features, cells.1, drop = FALSE], 
                  MARGIN = 1, FUN = function(x) {
                    return(log(x = mean(x = expm1(x = x)) + pseudocount.use))
                  })
  data.2 <- apply(X = object[features, cells.2, drop = FALSE], 
                  MARGIN = 1, FUN = function(x) {
                    return(log(x = mean(x = expm1(x = x)) + pseudocount.use))
                  })
  total.diff <- (data.1 - data.2)
  features.diff <- if (only.pos) {
    names(x = which(x = total.diff > logfc.threshold))
  }
  else {
    names(x = which(x = abs(x = total.diff) > logfc.threshold))
  }
  features <- intersect(x = features, y = features.diff)
  if (length(x = features) == 0) {
    stop("No features pass logfc.threshold threshold")
  }
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    if (length(x = cells.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(x = cells.2) > max.cells.per.ident) {
      cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
    }
    if (!is.null(x = latent.vars)) {
      latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
    }
  }
  if (!(test.use %in% c("negbinom", "poisson", "MAST", "LR")) && 
      !is.null(x = latent.vars)) {
    warning("'latent.vars' is only used for 'negbinom', 'poisson', 'LR', and 'MAST' tests")
  }
  de.results <- switch(EXPR = test.use,
                       ve = VeDETest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], 
                                     cells.1 = cells.1, cells.2 = cells.2, verbose = verbose, ...), 
                       wilcox = Seurat:::WilcoxDETest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], 
                                                      cells.1 = cells.1, cells.2 = cells.2, verbose = verbose, ...), 
                       bimod = Seurat:::DiffExpTest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], 
                                                    cells.1 = cells.1, cells.2 = cells.2, verbose = verbose), 
                       roc = Seurat:::MarkerTest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], 
                                                 cells.1 = cells.1, cells.2 = cells.2, verbose = verbose), 
                       t = Seurat:::DiffTTest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], 
                                              cells.1 = cells.1, cells.2 = cells.2, verbose = verbose), 
                       negbinom = Seurat:::GLMDETest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], 
                                                     cells.1 = cells.1, cells.2 = cells.2, min.cells = min.cells.feature, 
                                                     latent.vars = latent.vars, test.use = test.use, verbose = verbose), 
                       poisson = Seurat:::GLMDETest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], 
                                                    cells.1 = cells.1, cells.2 = cells.2, min.cells = min.cells.feature, 
                                                    latent.vars = latent.vars, test.use = test.use, verbose = verbose), 
                       MAST = Seurat:::MASTDETest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], 
                                                  cells.1 = cells.1, cells.2 = cells.2, latent.vars = latent.vars, verbose = verbose), 
                       DESeq2 = Seurat:::DESeq2DETest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], 
                                                      cells.1 = cells.1, cells.2 = cells.2, verbose = verbose), 
                       LR = Seurat:::LRDETest(data.use = object[features, c(cells.1, cells.2), drop = FALSE], 
                                              cells.1 = cells.1, cells.2 = cells.2, latent.vars = latent.vars, verbose = verbose), 
                       stop("Unknown test: ", test.use))
  de.results[, "avg_logFC"] <- total.diff[rownames(x = de.results)]
  de.results <- cbind(de.results, data.alpha[rownames(x = de.results), 
                                             , drop = FALSE])
  de.results$p_val_adj = p.adjust(p = de.results$p_val, method = "bonferroni", 
                                  n = nrow(object))
  if (test.use == "roc") {
    de.results <- de.results[order(-de.results$power, -de.results$avg_logFC), ]
  }
  else {
    de.results <- de.results[order(de.results$p_val, -de.results$avg_logFC), ]
  }
  if (only.pos) {
    de.results <- subset(x = de.results, subset = avg_logFC >  0)
  }
  return(de.results)
}


NewFindMarkers.Seurat <- function (object, ident.1 = NULL, ident.2 = NULL, group.by = NULL, 
                                   subset.ident = NULL, assay = NULL, slot = "data", reduction = NULL, 
                                   features = NULL, logfc.threshold = 0.25, test.use = "wilcox", 
                                   min.pct = 0.1, min.diff.pct = -Inf, verbose = TRUE, only.pos = FALSE, 
                                   max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL, 
                                   min.cells.feature = 3, min.cells.group = 3, pseudocount.use = 1, 
                                   ...) 
{
  if (!is.null(x = group.by)) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    Idents(object = object) <- group.by
  }
  if (!is.null(x = assay) & !is.null(x = reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  data.slot <- ifelse(test = test.use %in% c("negbinom", "poisson", 
                                             "DESeq2"), yes = "counts", no = slot)
  if (is.null(x = reduction)) {
    assay <- Seurat:::`%||%`(assay, DefaultAssay(object = object))
    data.use <- GetAssayData(object = object[[assay]], slot = data.slot)
  }
  else {
    if (data.slot == "counts") {
      stop("The following tests cannot be used when specifying a reduction as they assume a count model: negbinom, poisson, DESeq2")
    }
    data.use <- t(x = Embeddings(object = object, reduction = reduction))
  }
  if (is.null(x = ident.1)) {
    stop("Please provide ident.1")
  }
  else if (ident.1 == "clustertree" || is(object = ident.1, 
                                          class2 = "phylo")) {
    if (is.null(x = ident.2)) {
      stop("Please pass a node to 'ident.2' to run FindMarkers on a tree")
    }
    tree <- if (is(object = ident.1, class2 = "phylo")) {
      ident.1
    }
    else {
      Tool(object = object, slot = "BuildClusterTree")
    }
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' or pass an object of class 'phylo' as 'ident.1'")
    }
    ident.1 <- tree$tip.label[GetLeftDescendants(tree = tree, 
                                                 node = ident.2)]
    ident.2 <- tree$tip.label[GetRightDescendants(tree = tree, 
                                                  node = ident.2)]
  }
  if (length(x = as.vector(x = ident.1)) > 1 && any(as.character(x = ident.1) %in% 
                                                    colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(x = !as.character(x = ident.1) %in% 
                                                colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.1 are not present in the object: ", 
                  paste(bad.cells, collapse = ", ")))
    }
  }
  else {
    ident.1 <- WhichCells(object = object, idents = ident.1)
  }
  if (length(x = as.vector(x = ident.2)) > 1 && any(as.character(x = ident.2) %in% 
                                                    colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(!as.character(x = ident.2) %in% 
                                                colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.2 are not present in the object: ", 
                  paste(bad.cells, collapse = ", ")))
    }
  }
  else {
    if (is.null(x = ident.2)) {
      ident.2 <- setdiff(x = colnames(x = data.use), y = ident.1)
    }
    else {
      ident.2 <- WhichCells(object = object, idents = ident.2)
    }
  }
  if (!is.null(x = latent.vars)) {
    latent.vars <- FetchData(object = object, vars = latent.vars, 
                             cells = c(ident.1, ident.2))
  }
  de.results <- NewFindMarkers.default(object = data.use, cells.1 = ident.1, 
                                       cells.2 = ident.2, features = features, logfc.threshold = logfc.threshold, 
                                       test.use = test.use, min.pct = min.pct, min.diff.pct = min.diff.pct, 
                                       verbose = verbose, only.pos = only.pos, max.cells.per.ident = max.cells.per.ident, 
                                       random.seed = random.seed, latent.vars = latent.vars, 
                                       min.cells.feature = min.cells.feature, min.cells.group = min.cells.group, 
                                       pseudocount.use = pseudocount.use, ...)
  return(de.results)
}