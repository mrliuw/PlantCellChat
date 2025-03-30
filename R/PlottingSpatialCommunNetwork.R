#' @title Plot Spatial Communication Network or Barplot
#'
#' @description
#' This function visualizes a spatial communication network based on a given `pcc_obj` (a spatial object) containing communication probabilities and related information. The network can be plotted in either a "dimplot" format (graph-based) or a "barplot" format (quantitative summary). It supports filtering based on various user-provided parameters such as specific signals, interactions, sources, targets, and others.
#'
#' @param pcc_obj A PlantCellChat object.
#' @param key.signal A specific signal of interest (optional). If provided, filters the network based on this signal.
#' @param key.lrpairs A specific list of ligand-receptor pairs (optional). If provided, filters the network based on these pairs.
#' @param key.source A specific source cell type or node (optional). Filters the network to show communication starting from this source.
#' @param key.target A specific target cell type or node (optional). Filters the network to show communication ending at this target.
#' @param input.color A vector of colors to use for plotting (optional). Defaults to a predefined color set.
#' @param input.gradient A vector of colors to use for gradient fill in the circle plot (optional).
#' @param top.n The top N edges to display, sorted by communication strength (optional).
#' @param edge.size The size of the edges representing communication strength.
#' @param node.size The size of the nodes (cell types or spatial points) in the network.
#' @param arrow.size The size of the arrows on directed edges (only applicable in graph-based plots).
#' @param arrow.width The width of the arrows on directed edges (only applicable in graph-based plots).
#' @param edge.alpha The transparency of the edges (values between 0 and 1).
#' @param fig.ratio The aspect ratio of the figure, typically 1 for square plots.
#' @param label.cex The size of the labels (e.g., axis labels, legends).
#' @param title.cex The size of the title in the plot.
#' @param label.distance The distance of the labels from the circle center.
#' @param plot.type The type of plot to generate: "dimplot", "barplot", or "circle".
#' @param show.title Whether to display the title of the plot. Defaults to `TRUE`.
#' @param coord.flip Whether to flip the coordinates for the bar plot (optional).
#'
#' @return A spatial communication network plot, either as a graph or a bar plot, depending on the chosen `plot.type`.
#'
#' @export
PlottingSpatialCommunNetwork <- function(pcc_obj,
                                         key.signal = NULL,
                                         key.lrpairs = NULL,
                                         key.source = NULL,
                                         key.target = NULL,
                                         input.color = NULL,
                                         input.gradient = c("#fb7676","#7681fb"),
                                         top.n = NULL,
                                         edge.size = 2,
                                         node.size = 3,
                                         arrow.size = 0.3,
                                         arrow.width = 0.8,
                                         edge.alpha = 0.5,
                                         fig.ratio = 1,
                                         label.cex = 10,
                                         title.cex = 11,
                                         label.distance = 0,
                                         plot.type = "dimplot",
                                         show.title = T,
                                         coord.flip = F) {
  coord <- pcc_obj@spatial$cell.embeddings
  celltype <- levels(PccIdents(pccob))

  info <- data.frame(
    Cell = rownames(coord),
    Celltype = as.character(PccIdents(pcc_obj))
  )

  df <- pcc_obj@spatial$sp.net$df

  if (!is.null(key.lrpairs) && !is.null(key.signal)) {
    stop("Both `key.lrpairs` and `key.signal` cannot be specified at the same time. Please provide only one of them.")
  }

  if (!is.null(key.signal)) {
    if (!key.signal %in% df$Signal) {
      stop("The specified `key.signal` parameter does not exist in the data.")
    }
    idx <- df %>%
      dplyr::filter(Signal == key.signal) %>%
      dplyr::group_by(Signal, Source, Target) %>%
      dplyr::summarise(Prob = sum(Prob), .groups = "drop") %>%
      dplyr::arrange(desc(Prob))
  } else if (!is.null(key.lrpairs)) {
    if (!key.lrpairs %in% df$Interaction_name) {
      stop("The specified `key.lrpairs` parameter does not exist in the data.")
    }
    idx <- df %>%
      dplyr::filter(Interaction_name %in% key.lrpairs) %>%
      dplyr::group_by(Interaction_name, Source, Target) %>%
      dplyr::summarise(Prob = sum(Prob), .groups = "drop") %>%
      dplyr::arrange(desc(Prob))
  } else {
    idx <- df %>%
      dplyr::group_by(Source, Target) %>%
      dplyr::summarise(Prob = sum(Prob), .groups = "drop") %>%
      dplyr::arrange(desc(Prob))
  }

  if (!is.null(key.source) & !is.null(key.target)) {
    if (!key.source %in% df$Source | !key.target %in% df$Target) {
      stop("The specified `key.source` or `key.target` parameter does not exist in the data.")
    }
    idx <- idx %>% dplyr::filter(Source %in% key.source & Target %in% key.target) %>%
      dplyr::arrange(desc(Prob))
  } else if (!is.null(key.source)) {
    if (!key.source %in% df$Source) {
      stop("The specified `key.source` parameter does not exist in the data.")
    }
    idx <- idx %>% dplyr::filter(Source %in% key.source) %>%
      dplyr::arrange(desc(Prob))
  } else if (!is.null(key.target)) {
    if (!key.target %in% df$Target) {
      stop("The specified `key.target` parameter does not exist in the data.")
    }
    idx <- idx %>% dplyr::filter(Target %in% key.target) %>%
      dplyr::arrange(desc(Prob))
  }

  if (!is.null(top.n)) {
    idx <- idx %>%
      dplyr::top_n(top.n, Prob)
  }

  default_colors <-   c('#ffed6f', '#1f7ab4', '#fdcde6', '#cdecc6', '#fd8d3c',
                        '#66c1a7', '#e4191c', '#a6cee5', '#bc80bd', '#a55628',
                        '#fcd1a2', "#7C96E2", "#E672C2", "#4DAC5C", "#16AACC")

  colors <- if (is.null(input.color)) default_colors else input.color

  if (plot.type == "dimplot") {
    graph <- make_empty_graph(n = nrow(coord))
    V(graph)$name <- rownames(coord)

    edges <- idx %>%
      select(Source, Target) %>%
      as.matrix() %>%
      t()

    graph <- add_edges(graph, edges)

    # vertexs
    cell_colors <- setNames(colors[seq_along(celltype)], celltype)
    vertex_colors <- cell_colors[as.character(PccIdents(pcc_obj))]

    V(graph)$color <- vertex_colors
    V(graph)$frame.color <- NA

    vertexs <- unique(c(idx$Source, idx$Target))
    V(graph)$size <- ifelse(V(graph)$name %in% vertexs,
                            node.size + 7,
                            node.size)

    # edges
    edge_colors <- V(graph)$color[match(idx$Source, V(graph)$name)]
    E(graph)$color <- adjustcolor(edge_colors, alpha.f = edge.alpha)
    E(graph)$width <- idx$Prob * edge.size
    E(graph)$curved <- 0.2

    plot(graph,
         layout = as.matrix(coord),
         edge.width = E(graph)$width,
         edge.color = E(graph)$color,
         edge.arrow.size = arrow.size,
         edge.arrow.width = arrow.width,
         edge.curved = E(graph)$curved,
         vertex.label = NA,
         vertex.size = V(graph)$size,
         asp = fig.ratio,
         arrow.mode = 1)

    if (show.title == TRUE) {
      if (!is.null(key.signal)) {
        title(main = paste(key.signal, "-mediated\nSpatial Communication Network"),
              cex.main = title.cex - 10, font.main = 1)
      } else if (!is.null(key.lrpairs)) {
        title(main = paste0(key.lrpairs, "-mediated\nSpatial Communication Network"),
              cex.main = title.cex - 10, font.main = 1)
      } else {
        title(main = "Spatial Communication Network",
              cex.main = title.cex - 10, font.main = 1)
      }
    }
  }

  if (plot.type == "barplot") {
    cell_colors <- setNames(colors[seq_along(celltype)], celltype)

    if (coord.flip == F) {
      plot <- idx %>%
        mutate(Pairs = paste(Source, " - ", Target)) %>%
        left_join(info, by = c("Source" = "Cell")) %>%
        rename("Source_cell" = Celltype) %>%
        left_join(info, by = c("Target" = "Cell")) %>%
        rename("Target_cell" = Celltype) %>%
        mutate(Source_color = cell_colors[Source_cell],
               Target_color = cell_colors[Target_cell]) %>%
        ggplot(aes(x = Prob, y = reorder(Pairs, Prob), fill = Source_cell)) +
        geom_bar(stat = "identity",width = 0.7) +
        scale_fill_manual(values = cell_colors) +
        theme_bw() +
        theme(
          axis.text.y = element_text(size = label.cex, color = "black"),
          axis.text.x = element_text(size = label.cex, color = "black"),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = title.cex, color = "black"),
          plot.title = element_text(hjust = 0.5, size = title.cex, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
        labs(x = "Communication Strength",
             y = "",
             fill = "Sender cell",
             title = if(show.title == TRUE) {
               paste("Top",top.n,"Interaction Cells")
             })
    } else {
      plot <- idx %>%
        mutate(Pairs = paste(Source, " - ", Target)) %>%
        left_join(info, by = c("Source" = "Cell")) %>%
        mutate(Color = cell_colors[Celltype]) %>%
        ggplot(aes(x = reorder(Pairs, -Prob), y = Prob, fill = Celltype)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = cell_colors) +
        scale_y_continuous(expand = c(0, 0),
                           limits = c(0,max(idx$Prob) + max(idx$Prob)/10)) +
        theme_bw() +
        theme(
          axis.text.y = element_text(size = label.cex, color = "black"),
          axis.text.x = element_text(angle = 80, hjust = 1,size = label.cex, color = "black"),
          axis.title.y = element_text(size = title.cex, color = "black"),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = title.cex, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
        labs(x = "",
             y = "Communication Strength",
             fill = "Sender cell",
             title = if(show.title == TRUE) {
               paste("Top",top.n,"Interaction Cells")
             })
    }
  }

  if (plot.type == "circle") {
    idx$id <- seq(nrow(idx),1)
    idx$angle <- 90 + 180/nrow(idx) + (360/nrow(idx))*(nrow(idx)-idx$id)
    idx$angle <- ifelse(idx$id > nrow(idx)/2, idx$angle+180, idx$angle)

    plot <- idx %>%
      mutate(Pairs = paste(Source, " - ", Target)) %>%
      ggplot(aes(x = reorder(Pairs, Prob), y = Prob, fill = Prob)) +
      geom_bar(stat = "identity", width = 0.9) +
      geom_text(aes(x = id, y = Prob + max(Prob) + label.distance, label = Pairs, angle = angle), size = label.cex - 7) +
      scale_fill_gradient(low = input.gradient[1], high = input.gradient[2]) +
      coord_polar(theta = "x") +
      ylim(-3,max(idx$Prob) + 8) +
      theme_bw() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = title.cex, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
      labs(fill = "Communication Strength",
           title = if(show.title == TRUE) {
             paste("Top",top.n,"Interaction Cells")
           })
  }
  print(plot)
}


