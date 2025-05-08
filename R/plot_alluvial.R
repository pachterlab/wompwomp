#' alluvialmatch: Cluster-matching alluvial plots
#'
#' Main plotting function and helpers for bipartite-matching-based alluvial diagrams
#' @docType package
#' @name alluvialmatch
#'
#' @importFrom dplyr mutate select group_by summarise arrange desc ungroup slice n pull filter
#' @importFrom ggplot2 ggplot aes geom_text scale_fill_manual labs after_stat annotate theme_void theme element_text rel ggsave
#' @importFrom ggalluvial geom_alluvium geom_stratum stat_stratum
#' @importFrom ggforce gather_set_data
#' @importFrom igraph max_bipartite_match V graph_from_data_frame
#' @importFrom tibble is_tibble
#' @importFrom utils read.csv
#' @importFrom stats setNames
#' @importFrom rlang sym .data
#' @importFrom magrittr %>%
#' @importFrom data.table :=

utils::globalVariables(c(
    ".data", ":=", "%>%", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "stratum", "total", "cum_y", "best_cluster_agreement"
))

# library(dplyr)
# library(ggplot2)
# library(ggalluvial)
# library(ggforce)
# library(igraph)
# library(tibble)

axis_text_size <- 1.7
axis_numbering_size <- 1.4

default_colors <- c(
    "#D55E00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#E69F00", "#CC79A7", "#666666", "#AD7700", "#1C91D4", "#007756", "#D5C711", "#005685",
    "#A04700", "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71", "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C", "#FFCB57", "#9AD2F2",
    "#2CFFC6", "#F6EF8E", "#38B7FF", "#FF9B4D", "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45", "#AA9F0D", "#00446B", "#803800", "#8D3666",
    "#3D3D3D", 'black', 'white', 'red', 'blue'
)

if (!exists("group1_color")) {
    group1_color <- "#D55E00"
}

if (!exists("group2_color")) {
    group2_color <- "#56B4E9"
}

StatStratum <- ggalluvial::StatStratum  # avoid the error Can't find stat called "stratum"

reorder_clusters_descending <- function(clusters) {
    # Count the size of each cluster
    cluster_sizes <- table(clusters)

    # Sort the clusters by size and get the ordered names
    ordered_cluster_names <- names(sort(cluster_sizes, decreasing = TRUE))

    # Create a mapping from old to new cluster numbers
    cluster_mapping <- setNames(seq_along(ordered_cluster_names), ordered_cluster_names)

    # Apply the mapping to renumber the clusters
    renumbered_clusters <- cluster_mapping[as.character(clusters)]

    # Convert the factor to numeric
    renumbered_clusters_factor <- factor(as.numeric(as.character(renumbered_clusters)))

    return(renumbered_clusters_factor)
}


increment_if_zeros <- function(clus_df_gather, column) {
    clus_df_gather <- clus_df_gather %>% mutate(group_numeric = as.numeric(as.character(.data[[column]])))

    if (any(clus_df_gather$group_numeric == 0, na.rm = TRUE)) {
        clus_df_gather$group_numeric <- clus_df_gather$group_numeric + 1
        clus_df_gather <- clus_df_gather %>% mutate(!!column := factor(group_numeric))
    }

    clus_df_gather <- clus_df_gather %>% select(-group_numeric)

    return(clus_df_gather)
}

sort_clusters_by_agreement <- function(clus_df_gather, stable_column = "A", reordered_column = "B") {
    for (n in 1:2) {
        clus_df_gather$y <- -2  # tmp
        reordered_column_original_clusters_name <- paste0(reordered_column, "_original_clusters")

        clus_df_gather <- increment_if_zeros(clus_df_gather, stable_column)
        clus_df_gather <- increment_if_zeros(clus_df_gather, reordered_column)
        clus_df_gather <- increment_if_zeros(clus_df_gather, "col2_int")

        # Initialize variables
        half_rows <- nrow(clus_df_gather) / 2

        subset_data <- clus_df_gather %>%
            ungroup() %>%
            dplyr::slice((half_rows + 1):nrow(clus_df_gather))

        subset_data <- subset_data %>% mutate(
            !!reordered_column_original_clusters_name := as.numeric(as.character(.data[[reordered_column]])),
            !!reordered_column := -as.numeric(as.character(.data[[reordered_column]])),
            y := -as.numeric(as.character(y)),
            best_cluster_agreement := .data[[reordered_column]]
        )

        # Loop over each unique cluster number in Group 2
        for (cluster_number in sort(unique(subset_data[[reordered_column_original_clusters_name]]))) {
            # Subset the data for the current cluster number
            subset_data2 <- subset_data[subset_data[[reordered_column_original_clusters_name]] == cluster_number, ]

            # Find the row with the largest overlap (value)
            largest_overlap_row <- subset_data2[which.max(subset_data2$value), ]

            # Check if the corresponding Group 1 number is available
            new_cluster_number <- as.numeric(as.character(largest_overlap_row[[stable_column]]))

            best_cluster_number <- new_cluster_number

            subset_data$best_cluster_agreement[subset_data[[reordered_column_original_clusters_name]] == cluster_number] <- best_cluster_number

            any(subset_data[[reordered_column]] == best_cluster_number)

            while (any(subset_data[[reordered_column]] == new_cluster_number)) {
                new_cluster_number <- new_cluster_number + 1
                if (!any(subset_data$best_cluster_agreement[subset_data[[reordered_column]] == new_cluster_number] <= best_cluster_number)) {
                    subset_data[[reordered_column]] <- ifelse(subset_data[[reordered_column]] >= new_cluster_number, subset_data[[reordered_column]] + 1, subset_data[[reordered_column]])
                    subset_data$y <- ifelse(subset_data$y >= new_cluster_number, subset_data$y + 1, subset_data$y)
                    break
                }
            }

            # Assign the new cluster number
            subset_data[[reordered_column]][subset_data[[reordered_column_original_clusters_name]] == cluster_number] <- new_cluster_number

            for (i in (new_cluster_number + 1):(max(as.numeric(as.character(clus_df_gather$y))))) {
                if (!any(subset_data[[reordered_column]] == i)) {
                    subset_data[[reordered_column]] <- ifelse(subset_data[[reordered_column]] > i, subset_data[[reordered_column]] - 1, subset_data[[reordered_column]])
                    subset_data$y <- ifelse(subset_data$y > i, subset_data$y - 1, subset_data$y)
                    break
                }
            }
        }

        mapping <- setNames(
            seq_along(sort(unique(subset_data[[reordered_column]]))),
            sort(unique(subset_data[[reordered_column]]))
        )

        subset_data <- subset_data %>% mutate(!!reordered_column := mapping[as.character(.data[[reordered_column]])])

        clus_df_gather[[reordered_column]][(1:half_rows)] <- subset_data[[reordered_column]]
        clus_df_gather[[reordered_column]][((half_rows + 1):nrow(clus_df_gather))] <- subset_data[[reordered_column]]
        clus_df_gather$y[((half_rows + 1):nrow(clus_df_gather))] <- subset_data[[reordered_column]]

        sorted_levels <- sort(as.numeric(levels(clus_df_gather[[reordered_column]])))
        sorted_levels <- as.character(sorted_levels)
        clus_df_gather[[reordered_column]] <- factor(clus_df_gather[[reordered_column]], levels = sorted_levels)

        sorted_levels <- sort(as.numeric(levels(clus_df_gather$y)))
        sorted_levels <- as.character(sorted_levels)
        clus_df_gather$y <- factor(clus_df_gather$y, levels = sorted_levels)

        # set clus_df_gather$y to col1_int for rows where x=1, and clus_df_gather$y to col2_int for rows where x=2
        clus_df_gather$y <- ifelse(clus_df_gather$x == 1,
                                   clus_df_gather[[stable_column]],
                                   clus_df_gather[[reordered_column]])

    }

    return(clus_df_gather)
}


find_group2_colors <- function(clus_df_gather, group1_name, group2_name, ditto_colors) {
    group1_name = 'col1_int'
    group2_name = 'col2_int'
    clus_df_filtered <- clus_df_gather[, c(group1_name, group2_name, "value")]

    clus_df_filtered[[group1_name]] <- paste0("G1_", clus_df_filtered[[group1_name]])
    clus_df_filtered[[group2_name]] <- paste0("G2_", clus_df_filtered[[group2_name]])

    g <- igraph::graph_from_data_frame(d = clus_df_filtered, directed = FALSE)
    igraph::V(g)$type <- ifelse(igraph::V(g)$name %in% clus_df_filtered[[group1_name]], TRUE, FALSE)

    matching <- igraph::max_bipartite_match(g, weights = clus_df_filtered$value)

    keys <- names(matching$matching)
    number_group1_clusters <- length(sub("^G1_", "", keys[grep("^G1_", keys)]))
    number_group2_clusters <- length(sub("^G2_", "", keys[grep("^G2_", keys)]))

    # Extract and filter the matching pairs
    non_na_pairs <- matching$matching[!is.na(matching$matching)]

    # Filter out pairs where S is matched to C (excluding C matched to S or NA)
    g1_to_g2_pairs <- non_na_pairs[grep("^G1_", names(non_na_pairs))]


    # Extract numeric indices from the filtered pairs
    g1_indices <- as.numeric(sub("G1_", "", names(g1_to_g2_pairs)))
    g2_indices <- as.numeric(sub("G2_", "", g1_to_g2_pairs))

    # Initialize the new colors vector
    group2_colors <- vector("character", number_group2_clusters)

    # Assign colors based on the matching
    for (i in seq_along(g1_indices)) {
        # Check if the C index is within the bounds of ditto_colors
        if (g2_indices[i] <= length(ditto_colors)) {
            group2_colors[g2_indices[i]] <- ditto_colors[g1_indices[i]]
        }
    }

    remaining_colors <- ditto_colors[number_group1_clusters+1:length(ditto_colors)]

    # logic: take parents color if fraction > .5
    for (i in which(group2_colors %in% c(''))){
        current_g2 <-  paste0("G2_", i)

        test <- clus_df_filtered[, c(group1_name, group2_name, "value")]
        test2 <- test[test[[group2_name]] == current_g2,]

        potential_g1 = test2[[group1_name]]

        if (length(unique(potential_g1)) > 1) {
            test2[['fraction']] <- test2['value']/sum(test2['value'])
            if (sum(test2[['fraction']] > .1) > 3) {
                group2_colors[i] <- remaining_colors[1]
                remaining_colors <- remaining_colors[2:length(remaining_colors)]
            } else if (max(test2[['fraction']]) > .5) {
                test3 <- test2[test2[['fraction']]==max(test2[['fraction']]),]
                potential_g1 = test3[[group1_name]]
                group2_colors[i] = ditto_colors[as.numeric(sub("G1_", "", unique(potential_g1)))]
            } else {
                group2_colors[i] <- remaining_colors[1]
                remaining_colors <- remaining_colors[2:length(remaining_colors)]
            }
        } else{
            group2_colors[i] = ditto_colors[as.numeric(sub("G1_", "", unique(potential_g1)))]
        }


    }
    #group2_colors[(group2_colors == "")] <- remaining_colors[1:sum(group2_colors == "")]


    return (group2_colors)
}

plot_alluvial_internal <- function(clus_df_gather, group1_name = "A", group2_name = "B", group1_name_mapping = "A", group2_name_mapping = "B", color_list = 'DEFAULT', color_boxes = TRUE, color_bands = FALSE, alluvial_alpha = 0.5, match_colors = TRUE, output_path = NULL, include_labels_in_boxes = FALSE, include_axis_titles = FALSE, include_group_sizes = FALSE, show_group_2_box_labels_in_ascending = FALSE) {
    if (!is.null(color_list)){
        ditto_colors <- color_list
    } else{
        ditto_colors <- default_colors
    }
    num_levels_group1 <- length(levels(clus_df_gather[[group1_name]]))
    num_levels_group2 <- length(levels(clus_df_gather[[group2_name]]))

    if (show_group_2_box_labels_in_ascending) {
        group2_name_mapping <- group2_name_mapping
    }

    # Extract colors for each factor, assuming ditto_colors is long enough
    colors_group1 <- ditto_colors[1:num_levels_group1]

    if (match_colors) {
        colors_group2 <- find_group2_colors(clus_df_gather, group1_name, group2_name, ditto_colors)
    } else {
        colors_group2 <- ditto_colors[(num_levels_group1+1):(num_levels_group1+num_levels_group2)]
    }

    colors_group1_reverse <- rev(colors_group1)
    colors_group2_reverse <- rev(colors_group2)

    # Combine the colors
    combined_colors <- c(colors_group1, colors_group2)
    combined_colors_reverse <- c(colors_group1_reverse, colors_group2_reverse)

    # uncomment to attempt mapping
    p <- ggplot(data = clus_df_gather, aes(axis1 = !!sym(group1_name_mapping), axis2 = !!sym(group2_name_mapping), y = value))
    # p <- ggplot(data = clus_df_gather, aes(axis1 = !!sym(group1_name), axis2 = !!sym(group2_name), y = value))

    if (color_bands) {
        if (num_levels_group2 > num_levels_group1) {
            p <- p +
                geom_alluvium(aes(fill = !!sym(group2_name)), alpha = alluvial_alpha) +
                scale_fill_manual(values = colors_group2) +
                labs(fill = NULL)
        } else {
            p <- p +
                geom_alluvium(aes(fill = !!sym(group1_name)), alpha = alluvial_alpha) +
                scale_fill_manual(values = colors_group1) +
                labs(fill = NULL)
        }
    } else {
        p <- p + geom_alluvium()
    }

    if (color_boxes) {
        p <- p + geom_stratum(fill = combined_colors_reverse)
    } else {
        p <- p + geom_stratum()
    }

    if (!(include_labels_in_boxes==FALSE)) {
        if (match_colors) {p <- p +
            geom_text(stat = StatStratum, aes(label = after_stat(stratum)), size = 3, color = "black")
        } else{
            combined_names <- c(levels(clus_df_gather[[group1_name]]), levels(clus_df_gather[[group2_name]]))
            label_names <- c()
            for (i in seq_along(combined_colors)){
                current_color <- combined_colors[i]
                if (!(current_color %in% combined_colors[0:(i-1)])){
                    label_names <- c(label_names, combined_names[i])
                } else{
                    label_names <- c(label_names,label_names[which(combined_colors == current_color)[1]])
                }
            }
            A_label_names <- rev(label_names[1:num_levels_group1])
            B_label_names <- rev(label_names[(1+num_levels_group1):length(label_names)])

            final_label_names <- c(A_label_names, B_label_names)
            p <- p +
                geom_text(stat = StatStratum, aes(label = after_stat(final_label_names)))
        }
    }

    top_y1 <- clus_df_gather %>%
        filter(x == 1) %>%
        group_by(y) %>%
        summarise(total = sum(value), .groups = "drop") %>%
        arrange(desc(total)) %>%
        mutate(cum_y = cumsum(total)) %>%
        pull(cum_y) %>%
        max()

    top_y2 <- clus_df_gather %>%
        filter(x == 2) %>%
        group_by(y) %>%
        summarise(total = sum(value), .groups = "drop") %>%
        arrange(desc(total)) %>%
        mutate(cum_y = cumsum(total)) %>%
        pull(cum_y) %>%
        max()

    top_y <- max(top_y1, top_y2)  # top_y1 and top_y2 are probably the same

    if (include_axis_titles) {
        # Offset to place labels a bit above
        offset <- 1.1 * top_y

        p <- p +
            annotate("text", x = 1, y = top_y + offset, label = group1_name, size = 5, hjust = 0.5) +
            annotate("text", x = 2, y = top_y + offset, label = group2_name, size = 5, hjust = 0.5)
    }

    if (include_group_sizes) {
        offset_below <- top_y * 0.075

        p <- p +
            annotate("text", x = 1, y = -offset_below, label = num_levels_group1, hjust = 0.5, size = 5) + # Adjust x, y for Seurat
            annotate("text", x = 2, y = -offset_below, label = num_levels_group2, hjust = 0.5, size = 5) # Adjust x, y for Scanpy
    }

    p <- p +
        # geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
        theme_void() +
        theme(
            text = element_text(family = "sans"),
            legend.text = element_text(size = rel(axis_text_size))
        )

    p <- p + theme(legend.position = "none")  # to hide legend

    if (!is.null(output_path)) {
        ggsave(output_path, plot = p, dpi = 300, bg = "white")
    }

    return(p)
}


get_alluvial_df <- function(df) {
    # Convert numeric clustering columns to ordered factors
    df <- df |>
        dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
        dplyr::group_by_all() |>
        dplyr::count(name = "value")
    ggforce::gather_set_data(df, 1:2)
}


#' Generate an Alluvial Plot with Minimal Cluster Cross-over
#'
#' Creates a two-axis alluvial plot to visualize the relationship between two categorical groupings (e.g., cluster assignments from different methods),
#' minimizing crossover and aligning clusters based on agreement.
#'
#' @param df A data frame, tibble, or CSV file path. Must contain at least two columns, each representing a clustering/grouping of the same entities (rows).
#' @param column1 Character. Name of the first column to plot. Optional if \code{df} has exactly two columns, or if \code{df} has exactly three columns including \code{column_weights}.
#' @param column2 Character. Name of the second column to plot. Optional if \code{df} has exactly two columns, or if \code{df} has exactly three columns including \code{column_weights}.
#' @param show_group_2_box_labels_in_ascending Logical. If \code{TRUE}, forces the group 2 axis to be displayed in ascending label order rather than greedy matching.
#' @param color_boxes Logical. Whether to color the rectangular strata boxes representing groups.
#' @param color_bands Logical. Whether to color the alluvial bands connecting the groups.
#' @param match_colors Logical. If \code{TRUE}, assigns consistent colors between column1 and column2 where matched.
#' @param alluvial_alpha Numeric between 0 and 1. Transparency level for the alluvial bands.
#' @param include_labels_in_boxes Logical. Whether to include text labels inside the rectangular group boxes.
#' @param include_axis_titles Logical. Whether to display axis titles for column1 and column2.
#' @param include_group_sizes Logical. If \code{TRUE}, includes group sizes in the labels (e.g., "Group A (42)").
#' @param column_weights Optional numeric vector of weights (same length as number of rows in \code{df}) to weight each row differently when calculating flows.
#' @param output_path Character. File path to save the plot (e.g., "plot.png"). If \code{NULL}, the plot is not saved.
#' @param color_list Optional named list or vector of colors to override default group colors.
#' @param return_greedy_wolf Logical. If \code{TRUE}, returns the underlying bipartite match data frame instead of a ggplot object.
#'
#' @return A \code{ggplot2} object representing the alluvial plot.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' plot_alluvial(df)
#' }
#'
#' @export
plot_alluvial <- function(df, column1 = NULL, column2 = NULL, show_group_2_box_labels_in_ascending = FALSE, color_boxes = TRUE, color_bands = TRUE, match_colors = TRUE, alluvial_alpha = 0.5, include_labels_in_boxes = TRUE, include_axis_titles = TRUE, include_group_sizes = TRUE, column_weights = NULL, output_path = NULL, color_list = NULL, return_greedy_wolf = FALSE) {
    if (is.character(df) && grepl("\\.csv$", df)) {
        df <- read.csv(df)  # load in CSV as dataframe
    } else if (tibble::is_tibble(df)) {
        df <- as.data.frame(df)  # convert tibble to dataframe
    } else if (!is.data.frame(df)) {
        stop("Input must be a data frame, tibble, or CSV file path.")
    }

    if (!is.null(column_weights)) {
        if (!(column_weights %in% colnames(df))) {
            stop(sprintf("column_weights '%s' is not a column in the dataframe.", column_weights))
        }
        df <- tidyr::uncount(df, weights = !!rlang::sym(column_weights))
    }

    if (ncol(df) <= 1) {
        stop(sprintf("Dataframe has %d columns. It must have at least two columns.", ncol(df)))
    } else if (ncol(df) == 2) {
        if (is.null(column1) && is.null(column2)) {
            column1 <- colnames(df)[1]
            column2 <- colnames(df)[2]
        } else if (is.null(column1)) {
            column1 <- setdiff(colnames(df), column2)
        } else if (is.null(column2)) {
            column2 <- setdiff(colnames(df), column1)
        }
    } else if (ncol(df) > 2) {
        if (is.null(column1) || is.null(column2)) {
            stop("Dataframe has more than two columns. Please specify column1 and column2 for the plot.")
        }
        warning(sprintf(
            "Dataframe has %d columns. Only columns %s and %s will be used.",
            ncol(df), column1, column2
        ))
    }

    if (!(column1 %in% colnames(df))) {
        stop(sprintf("column1 '%s' is not a column in the dataframe.", column1))
    }
    if (!(column2 %in% colnames(df))) {
        stop(sprintf("column2 '%s' is not a column in the dataframe.", column2))
    }

    df[['col1_int']] <- as.integer(as.factor(df[[column1]]))
    df[['col2_int']] <- as.integer(as.factor(df[[column2]]))

    #factorize input columns
    df[[column1]] <- as.factor(df[[column1]])
    df[[column2]] <- as.factor(df[[column2]])

    clus_df_gather <- get_alluvial_df(df)

    # clus_df_gather <- clus_df_gather %>% mutate(
    #     group1_column_original_clusters := as.character(.data[[column1]]),
    #     group2_column_original_clusters := as.character(.data[[column2]])
    # )

    # clus_df_gather <- sort_clusters_by_agreement(clus_df_gather, stable_column = column1, reordered_column = column2)
    clus_df_gather <- sort_clusters_by_agreement(clus_df_gather, stable_column = 'col1_int', reordered_column = 'col2_int')
    clus_df_gather[[column2]] <- clus_df_gather$col2_int

    if (return_greedy_wolf) {
        clus_df_gather <- clus_df_gather %>%
            ungroup() %>%
            slice(1:(n() %/% 2)) %>%              # keep first half of rows
            select(-id, -x, -y, -col1_int, -col2_int)         # drop columns
        return(clus_df_gather)
    }

    # clus_df_gather <- merge(
    #     clus_df_gather,
    #     adata_obs[, c(col1_int, column1)],
    #     by.x = column1,
    #     by.y = col1_int,
    #     all.x = TRUE
    # )
    #
    # clus_df_gather <- merge(
    #     clus_df_gather,
    #     adata_obs[, c(col2_int, column2)],
    #     by.x = column2,
    #     by.y = col2_int,
    #     all.x = TRUE
    # )
    #
    # clus_df_gather <- clus_df_gather %>%
    #     mutate(
    #         column1_mapping = factor(group1_column_original_clusters, levels = unique(group1_column_original_clusters[order(column1)])),
    #         column2_mapping = factor(group2_column_original_clusters, levels = unique(group2_column_original_clusters[order(column2)]))
    #     )

    alluvial_plot <- plot_alluvial_internal(clus_df_gather, group1_name = column1, group2_name = column2, group1_name_mapping = column1, group2_name_mapping = column2, color_list = color_list, color_boxes = color_boxes, color_bands = color_bands, match_colors = match_colors, alluvial_alpha = alluvial_alpha, include_labels_in_boxes = include_labels_in_boxes, include_axis_titles = include_axis_titles, include_group_sizes = include_group_sizes, show_group_2_box_labels_in_ascending = show_group_2_box_labels_in_ascending, output_path = output_path)

    return(alluvial_plot)
}


#' Perform a greedy algorithm heuristic for the Weighted One Layer Free problem
#'
#' Perform a greedy algorithm heuristic for the Weighted One Layer Free problem
#'
#' @param df A data frame, tibble, or CSV file path. Must contain at least two columns, each representing a clustering/grouping of the same entities (rows).
#' @param column1 Character. Name of the first column to plot. Optional if \code{df} has exactly two columns, or if \code{df} has exactly three columns including \code{column_weights}.
#' @param column2 Character. Name of the second column to plot. Optional if \code{df} has exactly two columns, or if \code{df} has exactly three columns including \code{column_weights}
#' @param column_weights Optional numeric vector of weights (same length as number of rows in \code{df}) to weight each row differently when calculating flows.
#'
#' @return A data frame.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' greedy_wolf(df)
#' }
#'
#' @export
greedy_wolf <- function(df, column1 = NULL, column2 = NULL, column_weights = NULL) {
    clus_df_gather <-  plot_alluvial(df = df, column1 = column1, column2 = column2, column_weights = column_weights, return_greedy_wolf = TRUE)
    return(clus_df_gather)
}
