#' @title Check if Seurat package is installed.
#'
#' @return The result of judgement.
#' @export
#'
#' @examples error_if_no_Seurat()
.error_if_no_Seurat <- function(){
  if(!requireNamespace("Seurat",quietly = T)){
    stop("The Seurat package is not installed!")
  }
}

#' @title Class definitions
#' @importFrom methods setClassUnion
#' @importClassesFrom Matrix dgCMatrix

methods::setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))
methods::setClassUnion(name = 'AnyFactor', members = c("factor", "list"))

#' @title Create PlantCellChat class
#'
#' @slot data.raw raw count data matrix
#' @slot data normalized data matrix for PlantCellChat analysis (Genes should be in rows and cells in columns)
#' @slot data.average A matrix of average expression of each cell group
#' @slot data.signaling a subset of normalized matrix only containing signaling genes
#' @slot net Three-dimensional matrix (K x K x N) representing the communication strength between different cell group for each ligand-receptor pairs
#' @slot netSignal Three-dimensional(K x K x M) array recording the communication strength between different cell group for signals.
#' @slot database ligand-receptor interaction database used in the analysis
#' @slot lrpairs subset data of ligand-receptor interaction database for signaling
#' @slot meta data frame storing the information associated with each cell
#' @slot idents the default label defining the cell identity used for all analysis
#' @slot diffexp A list: one element is a vector consisting of the identified over-expressed signaling genes; one element is a data frame returned from the differential expression analysis
#' @slot spatial description A list containing spatial information for cells, such as spatial coordinates.
#'
#' @exportClass PlantCellChat
#' @importFrom methods setClass
PlantCellChat <- methods::setClass("PlantCellChat",
                                   slots = c(data.raw = 'AnyMatrix',
                                             data = 'AnyMatrix',
                                             data.average = 'AnyMatrix',
                                             data.signaling = "AnyMatrix",
                                             net = "list",
                                             netSignal = "list",
                                             meta = "data.frame",
                                             idents = "AnyFactor",
                                             database = "list",
                                             lrpairs = "list",
                                             diffexp = "list",
                                             spatial = "list"))

#' @title Create PlantCellChat object to analysis cell communication.
#'
#' @description This function creates a PlantCellChat object for cell communication analysis,
#'     modified from the package `CellChat` https://github.com/sqjin/CellChat/blob/master/R/CellChat_class.R
#'
#' @param object A normalized data matrix or Seurat object.
#' @param meta Cell identity information.
#' @param input.ident A character name of the variable in meta data, defining cell groups.
#' @param assay Assay to use when the input is a Seurat object.
#'
#' @return A Plantcellchat object
#'
#' @export
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat GetAssayData
#' @importFrom methods new
#'
CreatePlantCellChat <- function(object,
                                meta = NULL,
                                input.ident = NULL,
                                assay = NULL,
                                spatial.coord = NULL) {
  
  ## -------------------------------
  ## 1. Matrix input
  ## -------------------------------
  if (inherits(x = object, what = c("matrix", "Matrix", "dgCMatrix"))) {
    message("Create a PlantCellChat object from a data matrix.\n")
    data <- object
  }
  
  ## -------------------------------
  ## 2. Seurat input
  ## -------------------------------
  if (is(object, "Seurat")) {
    
    .error_if_no_Seurat()
    
    if (is.null(spatial.coord)) {
      message("Creating a PlantCellChat object from a Seurat object.\n")
    } else {
      message("Creating a PlantCellChat object from a Seurat object based on the spatial coordinates.\n")
    }
    
    ## assay
    if (is.null(assay)) {
      assay <- Seurat::DefaultAssay(object)
      message("The default assay `", assay, "` is used in the data layer.\n")
    }
    
    if (assay == "integrated") {
      stop("The data in the `integrated` assay is not suitable for PlantCellChat analysis! ",
           "Please use the 'RNA' or 'SCT' assay!")
    }
    
    ## -------------------------------
    ## Seurat v4 / v5 compatibility
    ## -------------------------------
    is_v5 <- packageVersion("SeuratObject") >= "5.0.0"
    
    if (is_v5) {
      # Seurat v5: use layer
      data <- Seurat::GetAssayData(
        object = object,
        assay  = assay,
        layer  = "data"
      )
    } else {
      # Seurat v4: use slot
      data <- Seurat::GetAssayData(
        object = object,
        assay  = assay,
        slot   = "data"
      )
    }
    
    ## data sanity check
    if (min(data) < 0) {
      stop("The data matrix contains negative values.")
    }
    
    ## meta
    if (is.null(meta)) {
      message("The `meta.data` slot in the Seurat object is used as cell meta information.\n")
      meta <- object@meta.data
    }
    
    meta$ident <- Seurat::Idents(object)
    
    if (is.null(input.ident)) {
      input.ident <- "ident"
    }
  }
  
  ## -------------------------------
  ## 3. Meta handling
  ## -------------------------------
  if (!is.null(meta)) {
    
    if (inherits(x = meta, what = c("matrix", "Matrix"))) {
      meta <- as.data.frame(meta)
    }
    
    if (!is.data.frame(meta)) {
      stop("The input `meta` should be a data frame.")
    }
    
    if (!identical(rownames(meta), colnames(data))) {
      warning("The cell barcodes in `meta` do not match expression matrix. ",
              "Resetting rownames of `meta` to match data.\n")
      rownames(meta) <- colnames(data)
    }
    
  } else {
    meta <- data.frame(meta)
  }
  
  ## -------------------------------
  ## 4. Create PlantCellChat object
  ## -------------------------------
  pcc_obj <- methods::new(
    Class = "PlantCellChat",
    data  = data,
    meta  = meta
  )
  
  ## -------------------------------
  ## 5. Spatial coordinates (optional)
  ## -------------------------------
  if (!is.null(spatial.coord)) {
    
    if (is.matrix(spatial.coord)) {
      pcc_obj@spatial$cell.embeddings <- spatial.coord
    } else if (is.data.frame(spatial.coord)) {
      pcc_obj@spatial$cell.embeddings <- as.matrix(spatial.coord)
    } else {
      stop("The `spatial.coord` must be either a matrix or a data frame.")
    }
  }
  
  ## -------------------------------
  ## 6. Set cell identity
  ## -------------------------------
  if (!is.null(meta) && nrow(meta) > 0) {
    
    if (!(input.ident %in% colnames(meta))) {
      stop("The `input.ident` is not a column name in the `meta`, ",
           "which will be used for cell grouping.")
    }
    
    suppressMessages(
      pcc_obj <- SetPccIdent(pcc_obj, input.ident)
    )
    
    message(
      "The cell groups used for PlantCellChat analysis are:\n",
      paste(levels(pcc_obj@idents), collapse = ","),
      "\n"
    )
  }
  
  ## -------------------------------
  ## 7. Done
  ## -------------------------------
  message(
    "------------The PlantCellChat Object has been successfully created.------------ ",
    Sys.time()
  )
  
  return(pcc_obj)
}

