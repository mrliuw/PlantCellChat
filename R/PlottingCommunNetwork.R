#' @title Plot Cell Communication Network
#'
#' @description
#' This function visualizes cell communication networks based on ligand-receptor pairs or signaling molecules.
#'     It provides options to plot either paracrine or autocrine communication patterns in circular or chord plot types.
#'
#'
#' @param pcc_obj A PlantCellChat object.
#' @param ligand.type A character string specifying the type of ligand, either "lrpairs" (ligand-receptor pairs) or "signal".
#' @param comm.pattern A character string specifying the communication pattern. Should be "paracrine" or "autocrine" if `ligand.type` is "lrpairs".
#' @param key.lrpairs A character string specifying the key ligand-receptor pair to filter the network. Only applicable if `ligand.type` is "lrpairs".
#' @param key.signal A character string specifying the key signaling molecule to filter the network. Required if `ligand.type` is "signal".
#' @param key.cluster A character string specifying specific clusters to include in the network. Only applicable when `comm.pattern` is "paracrine".
#' @param edge.size Numeric, controlling the thickness of edges based on communication probabilities.
#' @param edge.alpha Numeric, transparency level for edges, ranging from 0 to 1.
#' @param arrow.size Numeric, controlling the size of arrowheads on directed edges.
#' @param arrow.width Numeric, controlling the width of arrowheads on directed edges.
#' @param label.cex Numeric, scaling factor for vertex label font size.
#' @param label.distance Numeric, distance for vertex labels from the nodes in the chord diagram.
#' @param title.cex Numeric, scaling factor for title font size.
#' @param show.title Logical, indicating whether to display the plot title.
#' @param input.color A character vector of custom colors for nodes. Length should match the number of clusters. If NULL, default colors are used.
#' @param plot.type A character string specifying the type of plot to generate: "circle" for a circular plot or "chord" for a chord diagram.
#'
#' @return A network plot of the cell communication graph.
#' @export
#'
PlottingCommunNetwork <- function(pcc_obj,
                                  ligand.type,
                                  comm.pattern = NULL,
                                  key.lrpairs = NULL,
                                  key.signal = NULL,
                                  key.cluster = NULL,
                                  edge.size = 8,
                                  edge.alpha = 0.5,
                                  arrow.size = 0.2,
                                  arrow.width = 0.6,
                                  label.cex = 1,
                                  label.distance = -0.5,
                                  title.cex = 1.2,
                                  show.title = T,
                                  input.color = NULL,
                                  plot.type = "circle",
                                  normalize = T) {
  check.lrpairs <- pcc_obj@lrpairs$lr.signaling %>%
    dplyr::filter(Interaction_type == "Interacting protein") %>%
    dplyr::pull(Interaction_name)

  check.signal <- as.character(na.omit(unique(pcc_obj@lrpairs$lr.signaling$Signal)))

  default_colors <- c(
    "#B07836", "#B41BBC", "#16AACC", "#DF6182", "#FCC910",
    "#D9BC97", "#C72026", "#0F499C", "#3E83A7", "#7C96E2",
    "#4DAC5C", "#A1CC44", "#502B85", "#8DC791", "#EE4C1E",
    "#F88E1B", "#790E2B", "#E672C2", "#0C5936", "#54e5cd"
  )

  if (length(default_colors) < length(levels(PccIdents(pcc_obj)))) {
    stop("The length of `default_colors` is less than the number of clusters. Please use `input.color` parameter to provide a longer color vector.")
  }

  if (!(ligand.type %in% c("lrpairs", "signal"))) {
    stop("The `ligand.type` parameter must be either 'lrpairs' or 'signal'. Please provide a valid input.")
  }

  if (!plot.type %in% c("circle", "chord")) {
    stop("The `plot.type` parameter must be either 'circle' or 'chord'. Please provide a valid input.")
  }

  if (!is.null(input.color) && length(input.color) != length(levels(PccIdents(pcc_obj)))) {
    stop("The `input.color` parameter must have the same length as the number of clusters. Please provide a valid input color vector.")
  }

  if (ligand.type == "lrpairs") {
    if (is.null(comm.pattern)) {
      stop("The `comm.pattern` parameter must be specified when `ligand.type` = 'lrpairs'. Please provide either 'paracrine' or 'autocrine'.")
    }
    if (!comm.pattern %in% c("paracrine", "autocrine")) {
      stop("The `comm.pattern` parameter must be either 'paracrine' or 'autocrine' when `ligand.type` = 'lrpairs'.")
    }
    if (!is.null(key.lrpairs) && !is.null(key.signal)) {
      stop("Both `key.lrpairs` and `key.signal` cannot be specified at the same time when `ligand.type` = 'lrpairs'. Please provide only one of them.")
    }
    if (!is.null(key.lrpairs) && !(key.lrpairs %in% check.lrpairs)) {
      stop(paste("The `key.lrpairs`", key.lrpairs, "is not found in the data. Please provide a valid lrpairs."))
    }
    if (!is.null(key.signal) && !(key.signal %in% check.signal)) {
      stop(paste("The `key.signal`", key.signal, "is not found in the data. Please provide a valid signal."))
    }

    if (!is.null(key.lrpairs)) {
      df <- pcc_obj@net$df %>%
        dplyr::filter(Interaction_name == key.lrpairs)
    }
    if (!is.null(key.signal)) {
      df <- pcc_obj@net$df %>%
        dplyr::filter(Signal == key.signal)
    } else {
      df <- pcc_obj@net$df
    }

    idx <- df %>%
      dplyr::group_by(Target, Source) %>%
      dplyr::summarise(Prob = sum(Prob), .groups = 'drop')

    if (normalize == T) {
      # max-min
      min_prob <- min(idx$Prob)
      max_prob <- max(idx$Prob)
      idx$Prob <- (idx$Prob - min_prob) / (max_prob - min_prob)
    }

    # nodes <- unique(c(unique(edges$Source),unique(edges$Target)))
    nodes <- levels(PccIdents(pcc_obj))

    if (comm.pattern == "paracrine") {
      edges <- idx[idx$Source != idx$Target, c("Source", "Target", "Prob")] %>%
        dplyr::mutate(Source = factor(Source, levels = nodes),
                      Target = factor(Target, levels = nodes))

    }
    if (comm.pattern == "autocrine") {
      if (!is.null(key.cluster)) {
        stop("The `key.cluster` parameter should not be specified when `comm.pattern` = 'autocrine'.")
      }
      edges <- idx[idx$Source == idx$Target, c("Source", "Target", "Prob")] %>%
        dplyr::mutate(Source = factor(Source, levels = nodes),
                      Target = factor(Target, levels = nodes))
    }

    if (!is.null(key.cluster)) {
      edges <- edges %>%
        dplyr::filter(Source %in% key.cluster)
    }

    graph <- graph_from_data_frame(edges,directed = T,vertices = nodes)

    V(graph)$color <- c(input.color %||% default_colors)[1:length(nodes)]
    ncells <- table(PccIdents(pcc_obj))
    min_size <- 10
    max_size <- 20
    node_size <- (ncells - min(ncells)) / (max(ncells) - min(ncells)) * (max_size - min_size) + min_size
    V(graph)$size <- node_size[match(V(graph)$name, names(node_size))]
    V(graph)$label.color <- "black"
    V(graph)$label.cex <- label.cex
    V(graph)$frame.color <- NA

    edge_colors <- V(graph)$color[match(edges$Source, V(graph)$name)]

    E(graph)$color <- adjustcolor(edge_colors, alpha.f = edge.alpha)
    E(graph)$width <- ifelse(edges$Prob == 0, 0.000000001, edges$Prob * edge.size)
    E(graph)$curved <- 0.3

    layout_coords <- layout.circle(graph)
    angles <- -atan2(layout_coords[,2], layout_coords[,1]) + pi
    E(graph)$loop.angle <- ifelse(edges$Source == edges$Target,
                                  angles[as.numeric(as.factor(edges$Source))], 0)

    if (plot.type == "circle") {
      plot(graph,
           layout = layout.circle(graph),
           edge.width = E(graph)$width,
           edge.color = E(graph)$color,
           edge.curved = E(graph)$curved,
           vertex.size = V(graph)$size,
           vertex.label = V(graph)$name,
           vertex.label.cex = V(graph)$label.cex,
           edge.arrow.size = arrow.size,
           edge.arrow.width = arrow.width,
           arrow.mode = 2
      )
    }

    if (plot.type == "chord") {
      if (comm.pattern == "paracrine") {
        chord_matrix <- stats::xtabs(Prob ~ Source + Target, data = edges)

        # for (i in 1:nrow(chord_matrix)) {
        #   chord_matrix[i, ] <- chord_matrix[i, ] * ncells[i]/sum(ncells)
        #   chord_matrix[, i] <- chord_matrix[, i] * ncells[i]/sum(ncells)
        # }

        circos.clear()
        chordDiagram(chord_matrix,
                     grid.col = setNames(c(input.color %||% default_colors)[1:length(nodes)], nodes),
                     grid.border = "NA",
                     order = nodes,
                     annotationTrack = c("grid"),
                     directional = 1,
                     direction.type = c("arrows", "diffHeight"),
                     link.arr.length = 0.01,
                     link.arr.type = "big.arrow")
      }
      if (comm.pattern == "autocrine") {
        chord_matrix <- stats::xtabs(Prob ~ Source + Target, data = edges)

        chordDiagram(chord_matrix,
                     grid.col = setNames(c(input.color %||% default_colors)[1:length(nodes)], nodes),
                     grid.border = "NA",
                     order = nodes,
                     annotationTrack = c("grid"),
                     directional = 1,
                     direction.type = c("arrows", "diffHeight"),
                     link.arr.length = 0.2,
                     link.arr.type = "triangle")
      }
      suppressMessages({
        circos.track(track.index = 1,
                     panel.fun = function(x, y) {
                       xlim = get.cell.meta.data("xlim")
                       ylim = get.cell.meta.data("ylim")
                       sector.name = get.cell.meta.data("sector.index")
                       circos.text(mean(xlim),
                                   label.distance + ylim[2],
                                   sector.name,
                                   facing = "reverse.clockwise",
                                   niceFacing = TRUE,
                                   cex = label.cex,
                                   font = 1)
                     },
                     bg.border = NA)
      })
    }

    if (show.title == TRUE) {
      if (is.null(key.lrpairs) && is.null(key.signal)) {
        if (comm.pattern == "paracrine") {
          title(main = paste(key.signal, "signal Paracrine Communication Network"),
                cex.main = title.cex, font.main = 1)
        }
        if (comm.pattern == "autocrine") {
          title(main = paste(key.signal, "signal Autocrine Communication Network"),
                cex.main = title.cex, font.main = 1)
        }
      }

      if (!is.null(key.lrpairs)) {
        if (comm.pattern == "paracrine") {
          title(main = paste0(key.lrpairs, "-mediated Paracrine Communication Network"),
                cex.main = title.cex, font.main = 1)
        }
        if (comm.pattern == "autocrine") {
          title(main = paste0(key.lrpairs, "-mediated Autocrine Communication Network"),
                cex.main = title.cex, font.main = 1)
        }
      }

      if (!is.null(key.signal)) {
        if (comm.pattern == "paracrine") {
          title(main = paste0(key.signal, "-mediated Paracrine Communication Network"),
                cex.main = title.cex, font.main = 1)
        }
        if (comm.pattern == "autocrine") {
          title(main = paste0(key.signal, "-mediated Autocrine Communication Network"),
                cex.main = title.cex, font.main = 1)
        }
      }
    }
  }

  if (ligand.type == "signal") {
    if (!is.null(comm.pattern)) {
      stop("The `comm.pattern` parameter should not be specified when `ligand.type` = 'signal'.")
    }
    if (!is.null(key.lrpairs)) {
      stop("The `key.lrpairs` parameter should not be specified when `ligand.type` = 'signal'.")
    }
    if (is.null(key.signal)) {
      stop("The `key.signal` parameter must be specified when `ligand.type` = 'signal'. Please provide a valid signal.")
    }
    if (!is.null(key.signal) && !(key.signal %in% check.signal)) {
      stop(paste("The `key.signal`", key.signal, "is not found in the data. Please provide a valid signal."))
    }
    if (!is.null(key.cluster)) {
      stop("The `key.cluster` parameter should not be specified when `ligand.type` = 'signal'.")
    }

    if (!is.null(key.signal)) {
      df <- pcc_obj@netSignal$df %>%
        dplyr::filter(Signal == key.signal)

      idx <- df %>%
        dplyr::group_by(Target) %>%
        dplyr::summarise(Prob = sum(Prob), .groups = 'drop') %>%
        dplyr::mutate(Source = key.signal) %>%
        dplyr::select(Source,Target,Prob)
    }

    if (normalize == T) {
      # max-min
      min_prob <- min(idx$Prob)
      max_prob <- max(idx$Prob)
      idx$Prob <- (idx$Prob - min_prob) / (max_prob - min_prob)
    }

    nodes <- c(unique(idx$Source),levels(PccIdents(pcc_obj)))
    edges <- idx[,c("Source","Target","Prob")] %>%
      mutate(Source = factor(Source, levels = nodes),
             Target = factor(Target, levels = nodes))

    graph <- graph_from_data_frame(edges,directed = T,vertices = nodes)

    V(graph)$color <- c("#DECCDF", input.color %||% default_colors)[1:length(nodes)]
    ncells <- table(PccIdents(pcc_obj))
    min_size <- 10
    max_size <- 20
    node_size <- (ncells - min(ncells)) / (max(ncells) - min(ncells)) * (max_size - min_size) + min_size
    V(graph)$size <- ifelse(V(graph)$name == key.signal, max_size,
                            node_size[match(V(graph)$name, names(node_size))])
    V(graph)$label.color <- "black"
    V(graph)$label.cex <- label.cex
    V(graph)$frame.color <- NA

    edge_colors <- V(graph)$color[match(edges$Target, V(graph)$name)]
    E(graph)$color <- adjustcolor(edge_colors, alpha.f = edge.alpha)
    E(graph)$width <- ifelse(edges$Prob == 0, 0.000000001, edges$Prob * edge.size)
    E(graph)$curved <- 0.3

    if (plot.type == "circle") {
      plot(graph,
           layout = layout_as_star(graph),
           edge.width = E(graph)$width,
           edge.color = E(graph)$color,
           edge.curved = E(graph)$curved,
           vertex.size = V(graph)$size,
           vertex.label = V(graph)$name,
           vertex.label.cex = V(graph)$label.cex,
           edge.arrow.size = arrow.size,
           edge.arrow.width = arrow.width,
           arrow.mode = 2)
    }

    if (plot.type == "chord") {
      edges$Prob <- ifelse(edges$Prob == 0, 1e-10, edges$Prob)

      chord_matrix <- matrix(edges %>%
                               arrange(Target) %>%
                               pull(Prob), nrow = 1, byrow = TRUE)
      rownames(chord_matrix) <- nodes[1]
      colnames(chord_matrix) <- nodes[-1]

      chord_matrix <- t(chord_matrix)

      circos.clear()
      chordDiagram(
        chord_matrix,
        grid.col = setNames(c("#DECCDF", input.color %||% default_colors)[1:length(nodes)], nodes),
        grid.border = NA,
        order = nodes,
        annotationTrack = c("grid"),
        directional = -1,
        direction.type = c("diffHeight")
      )

      suppressMessages({
        circos.track(track.index = 1,
                     panel.fun = function(x, y) {
                       xlim = get.cell.meta.data("xlim")
                       ylim = get.cell.meta.data("ylim")
                       sector.name = get.cell.meta.data("sector.index")
                       circos.text(mean(xlim),
                                   label.distance + ylim[2],
                                   sector.name,
                                   facing = "reverse.clockwise",
                                   niceFacing = TRUE,
                                   cex = label.cex,
                                   font = 1)
                     },
                     bg.border = NA)
      })
    }
    if (show.title == T) {
      title(paste(key.signal, "signaling molecule Communication Network"),
            cex.main = title.cex, font.main = 1)
    }
  }
}


