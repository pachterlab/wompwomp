#' alluvialmatch: Cluster-matching alluvial plots
#'
#' Main plotting function and helpers for bipartite-matching-based alluvial diagrams
#' @docType package
#' @name alluvialmatch
#'
#' @importFrom dplyr mutate group_by summarise arrange desc ungroup slice n pull filter
#' @importFrom tibble is_tibble
#' @importFrom utils read.csv
#' @importFrom rlang sym .data
#' @importFrom magrittr %>%
#' @importFrom data.table :=

utils::globalVariables(c(
    ".data", ":=", "%>%", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "total", "cum_y", "best_cluster_agreement"
))


#' Determine sum of products of overlapping edges
#'
#' @param crossing_edges A CSV path or list as outputted from determine_crossing_edges.
#' @param minimum_edge_weight Optional positive integer that represents the minimum edge weight to count an edge in the calculation.
#'
#' @return A data frame.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' greedy_wolf(df)
#' crossing_edges <- determine_crossing_edges(df, column1 = "method1", column2 = "method2")
#' objective <- determine_weighted_layer_free_objective(crossing_edges)
#' }
#'
#' @export
determine_weighted_layer_free_objective <- function(crossing_edges, minimum_edge_weight = 0) {
    # Read the file (tab-separated or CSV)
    if (is.character(crossing_edges) && length(crossing_edges) == 1 && file.exists(crossing_edges)) {
        # df_in <- read.table(crossing_edges, sep="\t", header = TRUE, stringsAsFactors = FALSE)
        df_in <- read.csv(crossing_edges, stringsAsFactors = FALSE)

        if (ncol(df_in) != 6) {
            stop("File must contain exactly 6 columns: left1, right1, weight1, left2, right2, weight2")
        }

        crossing_edges <- apply(df_in, 1, function(row) {
            edge1 <- as.character(row[1:3])
            edge2 <- as.character(row[4:6])
            list(edge1, edge2)
        })

        crossing_edges <- unname(as.list(crossing_edges))

        # Case 2: Already a list of edge pairs
    } else if (is.list(crossing_edges)) {
        # do nothing
    } else
        stop("Input must be either a file path or a list.")

    total_weighted_crossings <- 0
    for (pair in crossing_edges) {
        w1 <- as.numeric(pair[[1]][3])
        w2 <- as.numeric(pair[[2]][3])

        if (is.na(w1) || is.na(w2) || w1 < minimum_edge_weight || w2 < minimum_edge_weight) {
            next
        }

        total_weighted_crossings <- total_weighted_crossings + w1 * w2
    }

    # Correct for double-counting
    total_weighted_crossings <- total_weighted_crossings / 2

    return(total_weighted_crossings)
}


# list(
#    list(c(l1, r1, e1), c(l2, r2, e2)),
#    ...
# )
#' Determine overlapping edges of graph run through data_sort
#'
#' @param df A data frame, tibble, or CSV file path. Must contain at least two columns, each representing a clustering/grouping of the same entities (rows).
#' @param column1 Character. Name of the first column to plot. Optional if \code{df} has exactly two columns, or if \code{df} has exactly three columns including \code{column_weights}.
#' @param column2 Character. Name of the second column to plot. Optional if \code{df} has exactly two columns, or if \code{df} has exactly three columns including \code{column_weights}
#' @param column_weights Optional numeric vector of weights (same length as number of rows in \code{df}) to weight each row differently when calculating flows.
#' @param minimum_edge_weight Optional positive integer that represents the minimum edge weight to count an edge in the calculation.
#' @param output_df_path Character. File path to save the dataframe (e.g., "df.csv"). If \code{NULL}, the dataframe is not saved.
#' @param return_weighted_layer_free_objective Bool Whether to return a list of overlapping edges (FALSE) or the sum of products of overlapping edges (TRUE)
#'
#' @return A data frame.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' df <- data_sort(df)
#' crossing_edges <- determine_crossing_edges(df, column1 = "col1_int", column2 = "col2_int")
#' }
#'
#' @export
determine_crossing_edges <- function(df, column1 = NULL, column2 = NULL, column_weights = "value", minimum_edge_weight = 0, output_df_path = NULL, map_dict = NULL, map_dict_1 = NULL, map_dict_2 = NULL, fixed_column = NULL, return_weighted_layer_free_objective = FALSE) {
    if (is.character(df)) {
        if (grepl("\\.rds$", df, ignore.case = TRUE)) {
            df <- readRDS(df)
        } else if (grepl("\\.csv$", df, ignore.case = TRUE)) {
            df <- read.csv(df)
        } else {
            stop("Input path must be a .csv or .rds file.")
        }
    } else if (!is.data.frame(df)) {
        stop("Input must be a data frame or a file path to a .csv or .rds file.")
    }

    if (!(column_weights %in% colnames(df))) {
        stop(sprintf("column_weights '%s' is not a column in the dataframe.", column_weights))
    }

    cluster_cols <- setdiff(colnames(df), column_weights)

    if (length(cluster_cols) <= 1) {
        stop(sprintf("Dataframe has %d columns. It must have at least three columns.", ncol(df) + 1))
    } else if (length(cluster_cols) == 2) {
        if (is.null(column1) && is.null(column2)) {
            column1 <- cluster_cols[1]
            column2 <- cluster_cols[2]
        } else if (is.null(column1)) {
            column1 <- setdiff(cluster_cols, column2)
        } else if (is.null(column2)) {
            column2 <- setdiff(cluster_cols, column1)
        }
    } else {
        if (is.null(column1) && is.null(column2)) {
            stop(sprintf("Dataframe has more than three columns. Please specify column1 and column2"))
        }
    }

    if (!(column1 %in% colnames(df))) {
        stop(sprintf("column1 '%s' is not a column in the dataframe.", column1))
    }
    if (!(column2 %in% colnames(df))) {
        stop(sprintf("column2 '%s' is not a column in the dataframe.", column2))
    }

    # # set to factors if not already
    # if (!is.factor(df[[column1]])) df[[column1]] <- factor(df[[column1]])
    # if (!is.factor(df[[column2]])) df[[column2]] <- factor(df[[column2]])

    # Assign fixed coordinates for bipartite layout
    df <- df %>%
        mutate(x1 = as.integer(!!sym(column1)),
               x2 = as.integer(!!sym(column2)),
               y1 = 0,
               y2 = 1)

    if (!is.null(map_dict)) {
        if (fixed_column == 1 || fixed_column == column1) {
            df$x2 <- map_dict[as.character(df$x2)]
        } else if (fixed_column == 2 || fixed_column == column2) {
            df$x1 <- map_dict[as.character(df$x1)]
        } else {
            stop("`fixed_column` must be 1, 2, column1, or column2")
        }
    }

    # Initialize result list and seen pair tracker
    crossing_edges <- list()
    n <- nrow(df)

    # Compare each pair of edges
    for (i in 1:n) {
        for (j in 1:n) {
            if (i != j) {
                if ((df$x1[i] < df$x1[j] && df$x2[i] > df$x2[j]) | (df$x1[i] > df$x1[j] && df$x2[i] < df$x2[j])) {
                    edge1 <- c(
                        as.character(df[[column1]][i]),
                        as.character(df[[column2]][i]),
                        as.character(df[[column_weights]][i])
                    )
                    edge2 <- c(
                        as.character(df[[column1]][j]),
                        as.character(df[[column2]][j]),
                        as.character(df[[column_weights]][j])
                    )

                    w1 <- as.numeric(edge1[3])
                    w2 <- as.numeric(edge2[3])
                    if (w1 < minimum_edge_weight || w2 < minimum_edge_weight) {
                        next
                    }

                    crossing_edges <- append(crossing_edges, list(list(edge1, edge2)))
                }
            }
        }
    }

    if (is.character(output_df_path) && grepl("\\.csv$", output_df_path, ignore.case = TRUE)) {
        df_out <- do.call(rbind, lapply(crossing_edges, function(pair) {
            c(pair[[1]], pair[[2]])
        }))
        colnames(df_out) <- c("left1", "right1", "weight1", "left2", "right2", "weight2")

        df_out <- as.data.frame(df_out, stringsAsFactors = FALSE)
        write.csv(df_out, file = output_df_path, row.names = FALSE, quote = FALSE)
    }

    if (return_weighted_layer_free_objective) {
        WLF_objective <- determine_weighted_layer_free_objective(crossing_edges)
        return(WLF_objective)
    }

    return(crossing_edges)
}


plot_alluvial_internal_multicol <- function(clus_df_gather,graphing_columns,
                                            sorting_algorithm = NULL,fixed_column=NULL,
                                            color_list = NULL, color_boxes = TRUE,
                                            color_bands = FALSE, color_band_list = NULL,
                                            color_band_column=NULL, color_band_boundary=FALSE,
                                            alluvial_alpha = 0.5, match_colors = TRUE, output_plot_path = NULL,
                                            include_labels_in_boxes = FALSE, include_axis_titles = FALSE,
                                            include_group_sizes = FALSE
) {
    if (!is.null(color_list)){
        ditto_colors <- color_list
    } else{
        ditto_colors <- default_colors
    }

    # Extract colors for each factor, assuming ditto_colors is long enough
    if (match_colors) {
        remaining_colors <- ditto_colors
        first <- TRUE
        final_colors <- c()
        n <- 1
        for (col_group in graphing_columns) {
            num_levels <- length(levels(clus_df_gather[[col_group]]))
            if (first) {
                old_colors <- remaining_colors[1:num_levels]
                final_colors <- c(final_colors, rev(old_colors))
            } else {
                temp_colors <- find_group2_colors(clus_df_gather, remaining_colors,
                                                  group_1_name = paste0('col',n-1,'_int'), group_2_name = paste0('col',n,'_int'))
                remaining_colors <- remaining_colors[1:length(old_colors)]
                old_colors <- temp_colors
                final_colors <- c(final_colors, rev(old_colors))

            }
            n <- n+1
        }
    } else {
        remaining_colors <- ditto_colors
        final_colors <- c()
        for (col_group in graphing_columns) {
            num_levels <- length(levels(clus_df_gather[[col_group]]))
            old_colors <- remaining_colors[1:num_levels]
            final_colors <- c(final_colors, rev(old_colors))
            remaining_colors <- remaining_colors[num_levels:length(remaining_colors)]
        }
    }

    remaining_colors <- ditto_colors[!(ditto_colors %in% final_colors)]

    # remove duplicate dims
    temp_df <- clus_df_gather[1:as.integer(dim(clus_df_gather)[1]/2),1:dim(clus_df_gather)[2]]

    # uncomment to attempt mapping
    p <- ggplot(data = temp_df, aes(y = value),
    )
    for (x in seq_along(graphing_columns)) {
        p$mapping[[paste0('axis',x)]] = sym(paste0('col', x,'_int'))
    }

    if (color_bands) {
        if (!is.null(color_band_column)) {
            if (is.null(color_band_list)) {
                color_band_list <- final_colors
            }
            if (color_band_boundary){
                p <- p +
                    geom_alluvium(aes(fill = !!sym(color_band_column), color=!!sym(color_band_column)),
                                  alpha = alluvial_alpha) +
                    scale_fill_manual(values = color_band_list) + scale_color_manual(values = color_band_list)+
                    labs(fill = NULL)+guides(fill='none')

            } else{
                p <- p +
                    geom_alluvium(aes(fill = !!sym(color_band_column)), alpha = alluvial_alpha) +
                    scale_fill_manual(values = color_band_list) +
                    labs(fill = NULL)+guides(fill='none')
            }
        } else {
            colors_group1 <- rev(final_colors[1:length(levels(clus_df_gather[['col1_int']]))])
            if (color_band_boundary){
                p <- p +
                    geom_alluvium(aes(fill = !!sym('col1_int'), color = !!sym('col1_int')), alpha = alluvial_alpha) +
                    scale_fill_manual(values = colors_group1) + scale_color_manual(values = colors_group1)+
                    labs(fill = NULL)+guides(fill='none')
            } else{
                p <- p +
                    geom_alluvium(aes(fill = !!sym('col1_int')), alpha = alluvial_alpha) +
                    scale_fill_manual(values = colors_group1) +
                    labs(fill = NULL)+guides(fill='none')
            }
        }
    } else {
        if (color_band_boundary){
            p <- p + geom_alluvium(color='grey2',alpha = alluvial_alpha)
        } else{
            p <- p + geom_alluvium(alpha = alluvial_alpha)
        }
    }

    if (color_boxes) {
        p <- p + geom_stratum(fill = final_colors)
    } else {
        p <- p + geom_stratum()
    }

    if (!(include_labels_in_boxes==FALSE)) {
        final_label_names <- c()
        for (col_int in seq_along(graphing_columns)) {
            int_name <- paste0('col', col_int, '_int')
            group_name <- graphing_columns[[col_int]]

            curr_label <- as.character(unique(clus_df_gather[order(clus_df_gather[[int_name]]),][[group_name]]))

            final_label_names <- c(final_label_names, rev(curr_label))
        }
        p <- p +
            geom_text(stat = StatStratum, aes(label = after_stat(final_label_names)))
    }

    top_y = 0
    for (test_x in unique(clus_df_gather$x)) {
        curr_y <- clus_df_gather %>%
            filter(x == test_x) %>%
            group_by(y) %>%
            summarise(total = sum(value), .groups = "drop") %>%
            arrange(desc(total)) %>%
            mutate(cum_y = cumsum(total)) %>%
            pull(cum_y) %>%
            max()
        top_y <- max(curr_y/2, top_y/2)
    }# top_y1 and top_y2 are probably the same

    if (include_axis_titles) {
        # Offset to place labels a bit above
        offset <- 1.1 * top_y
        x<-1
        for (col_group in graphing_columns) {
            p <- p +
                annotate("text", x = x, y = top_y + offset, label = col_group, size = 5, hjust = 0.5)
            x <- x+1
        }
    }

    if (include_group_sizes) {
        offset_below <- top_y * 0.075
        x<-1
        for (col_group in graphing_columns) {
            p <- p +
                annotate("text", x = x, y = -offset_below, label = length(levels(clus_df_gather[[col_group]])), hjust = 0.5, size = 5) # Adjust x, y for Scanpy
            x <- x+1
        }
    }

    p <- p +
        theme_void() +
        theme(
            text = element_text(family = "sans"),
            legend.text = element_text(size = rel(axis_text_size))
        )

    p <- p + theme(legend.position = "none")  # to hide legend

    if (!is.null(output_plot_path)) {
        ggsave(output_plot_path, plot = p, dpi = 300, bg = "white")
    }

    return(p)
}


