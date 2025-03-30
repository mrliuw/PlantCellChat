#' @title Calculate the average expression.
#' @description Providing a pipeline to calculate the average gene expression in
#'     different cells,including three methods(average,median,quartile).The quartile
#'         method using formula "1/2Q2+1/4(Q1+Q3)" to calculate the average expression.
#'
#' @param pcc_obj A PlantCellChat object.
#' @param input.data A sparse matrix, using `data.signaling` slot as default.
#' @param methods The method(average,median,quartile) used to calculate the average expression.
#' @param return.object Logical, whether to return the updated PlantCellChat object.
#'
#' @return The updated PlantCellChat object with the average expression in the 'data.average' slot.
#' @export
#'
#' @importFrom Matrix rowMeans
#' @importFrom MatrixGenerics rowMedians
#' @importFrom MatrixGenerics rowQuantiles
#'
#' @examples pcc_obj <- CalculateAvgExp(pcc_obj,methods = "average")
CalculateAvgExp <- function(pcc_obj,input.data = NULL,methods = NULL,return.object = TRUE) {
  if (is.null(input.data)) {
    if (is(pcc_obj,"PlantCellChat")) {
      if (is.null(pcc_obj@data.signaling) || ncol(pcc_obj@data.signaling) == 0) {
        stop("The `data.signaling` slot is empty, please run `ExtractSignalingData` first.")
      }
      data <- pcc_obj@data.signaling
    } else {
      stop("Please provide a PlantCellChat object.")
    }
  } else {
    if (is(input.data,"dgCMatrix")) {
      data <- input.data
    } else {
      stop("Please provide a dgCMatrix matrix.")
    }
  }

  colnames(data) <- PccIdents(pcc_obj)
  celltype <- as.character(unique(pcc_obj@idents))
  exp <- data.frame(gene = rownames(data))

  if (!(methods %in% c("average","median","quartile"))) {
    stop("Please select a valid method to calculate the average expression.")
  }

  if (is.null(methods) || length(methods) == 0) {
    stop("Please select a method to calculate the average expression.")
  }

  if (methods == "average") {
    for (i in 1:length(celltype)) {
      subset <- data[,which(colnames(data) == celltype[i])]
      avgexp <- as.data.frame(Matrix::rowMeans(subset,na.rm = T))
      colname <- paste("AvgExp", celltype[i], sep = "_")
      exp[colname] <- avgexp
    }
    if (return.object == TRUE) {
      message("The average expression has been calculated using the `average` method. The expression matrix saved in the `data.average` slot.")
    } else {
      message("The average expression has been calculated using the `average` method.")
    }
  }

  if(methods == "median") {
    for (i in 1:length(celltype)) {
      subset <- data[,which(colnames(data) == celltype[i])]
      avgexp <- as.data.frame(MatrixGenerics::rowMedians(subset,na.rm = T))
      colname <- paste("AvgExp", celltype[i], sep = "_")
      exp[colname] <- avgexp
    }
    if (return.object == TRUE) {
      message("The average expression has been calculated using the `median` method. The expression matrix saved in the `data.average` slot.")
    } else {
      message("The average expression has been calculated using the `median` method.")
    }
  }

  if (methods == "quartile") {
    for (i in 1:length(celltype)) {
      subset <- data[,which(colnames(data) == celltype[i])]
      Q <- as.data.frame(MatrixGenerics::rowQuantiles(subset,na.rm = T))[,c(2,3,4)]
      colnames(Q) <- c("Q1","Q2","Q3")
      avgexp <- 1/2 * Q[,2] + 1/4 * (Q[,1] + Q[,3])
      colname <- paste("AvgExp", celltype[i], sep = "_")
      exp[colname] <- avgexp
    }
    if (return.object == TRUE) {
      message("The average expression has been calculated using the `quartile` method. The expression matrix saved in the `data.average` slot.")
    } else {
      message("The average expression has been calculated using the `quartile` method.")
    }
  }

  exp <- exp[,-1]
  exp_matrix <- matrix(as.numeric(as.matrix(exp)),ncol = length(celltype),nrow = nrow(data),dimnames = list(rownames(data),celltype))
  sorted_exp_matrix <- exp_matrix[,levels(PccIdents(pcc_obj))]

  if (return.object == TRUE) {
    pcc_obj@data.average <- sorted_exp_matrix
    return(pcc_obj)
  } else {
    return(sorted_exp_matrix)
  }
}


