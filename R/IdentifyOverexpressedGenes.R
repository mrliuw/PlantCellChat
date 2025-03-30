#' @title Identify over-expressed signaling genes associated with each cell group.
#'
#' @description
#' This function is modified according to the function `IdentifyOverExpressedGenes` from
#'     the package 'CellChat' (https://github.com/sqjin/CellChat/blob/master/R/utilities.R).
#'
#' @param pcc_obj A PlantCellChat object.
#' @param input.data A sparse matrix, using `data.signaling` slot as default.
#' @param input.ident The identifier to be set as the default label, using `idents` slot as default.
#' @param thresh.pc Threshold of the percent of cells expressed in one cluster
#' @param thresh.fc Threshold of Log Fold Change
#' @param thresh.p Threshold of p-values
#'
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom stats sd wilcox.test
#' @importFrom stats p.adjust
#' @importFrom Matrix rowSums
#'
#' @return A data frame recording the results of differential expression analysis
#'     with each cell group and A character vector recording the names of the over-expression
#'         genes which saved in the `diffexp` slot of the PlantCellChat object.
#'
#' @export
#'
IdentifyOverExpressedGenes <- function(pcc_obj,input.data = NULL,input.ident = NULL,thresh.pc = 0,thresh.fc = 0,thresh.p = 0.05) {
  if (is.null(input.data)) {
    data <- pcc_obj@data.signaling
    if (nrow(data) < 3) {
      stop("Please check `data.signaling` slot and make sure run `ExtractSignalingData` first")
    }
  } else {
    if ("dgcMatrix" %in% class(input.data)) {
      data <- input.data
    }
  }

  features <- row.names(data)

  if (is.null(input.ident)) {
    idents <- factor(pcc_obj@idents)
  } else {
    idents <- factor(pcc_obj@meta[[input.ident]])
  }

  celltype <- levels(idents)[levels(idents) %in% unique(idents)]

  numCluster <- length(celltype)

  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )

  mean.fxn <- function(x) {
    return(log(mean(expm1(x)) + 1))
  }
  idents <- as.character(idents)
  genes <- vector("list", length = numCluster)
  for (i in 1:numCluster) {
    cell.use1 <- which(idents == celltype[i])
    cell.use2 <- base::setdiff(1:length(idents), cell.use1)

    thresh.min <- 0
    pct.1 <- round(Matrix::rowSums(data[features, cell.use1, drop = FALSE] > thresh.min) / length(cell.use1),
                   digits = 3)
    pct.2 <- round(Matrix::rowSums(data[features, cell.use2, drop = FALSE] > thresh.min) / length(cell.use2),
                   digits = 3)
    data.alpha <- cbind(pct.1, pct.2)
    colnames(data.alpha) <- c("pct.1", "pct.2")
    alpha.min <- apply(data.alpha, MARGIN = 1, FUN = max)
    names(alpha.min) <- rownames(data.alpha)
    features <- names(which(alpha.min > thresh.pc))
    if (length(features) == 0) {
      next
    }

    data.1 <- apply(data[features, cell.use1, drop = FALSE], MARGIN = 1, FUN = mean.fxn)
    data.2 <- apply(data[features, cell.use2, drop = FALSE], MARGIN = 1, FUN = mean.fxn)
    FC <- (data.1 - data.2)
    features.diff <- names(which(abs(FC) > thresh.fc))

    features <- intersect(features,features.diff)
    if (length(features) == 0) {
      next
    }

    data1 <- data[features, cell.use1, drop = FALSE]
    data2 <- data[features, cell.use2, drop = FALSE]

    pvalues <- unlist(
      x = my.sapply(
        X = 1:nrow(x = data1),
        FUN = function(x) {
          return(wilcox.test(data1[x, ], data2[x, ])$p.value)
        }
      )
    )

    padj = stats::p.adjust(
      p = pvalues,
      method = "bonferroni",
      n = nrow(data)
    )
    genes[[i]] <- data.frame(clusters = celltype[i], features = as.character(rownames(data1)), pvalues = pvalues, logFC = FC[features], data.alpha[features,, drop = F],pvalues.adj = padj, stringsAsFactors = FALSE)
  }

  marker <- data.frame()
  for (i in 1:numCluster) {
    deg <- genes[[i]]
    if (!is.null(deg)) {
      deg <- deg[order(deg$pvalues, -deg$logFC), ]
      deg <- subset(deg, subset = pvalues < thresh.p)
      if (nrow(deg) > 0) {
        marker <- rbind(marker, deg)
      }
    }
  }

  pcc_obj@diffexp$features <- as.character(marker$features)
  pcc_obj@diffexp$features.info <- marker

  message(paste0(length(marker$features)," over-expressed genes have been identified successfully. The differential expression analysis result saved in the `diffexp` slot."))

  return(pcc_obj)
}
