#' alluvialmatch: Cluster-matching alluvial plots
#'
#' Main plotting function and helpers for bipartite-matching-based alluvial diagrams
#' @docType package
#' @name alluvialmatch
#'
#' @importFrom dplyr mutate select group_by summarise arrange desc ungroup slice n pull filter
#' @importFrom ggplot2 ggplot aes geom_text scale_fill_manual labs after_stat annotate theme_void theme element_text rel ggsave guides
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

StatStratum <- ggalluvial::StatStratum  # avoid the error Can't find stat called "stratum"

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
    "#3D3D3D"
)

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


find_group2_colors <- function(clus_df_gather, ditto_colors,
                               group1_name = 'col1_int', group2_name = 'col2_int') {
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
    return (group2_colors)
}

plot_alluvial_internal_2col <- function(clus_df_gather,
                                        sorting_algorithm = NULL,
                                        group1_name = "A", group2_name = "B", fixed_column=NULL,
                                        color_list = NULL, color_boxes = TRUE,
                                        color_bands = FALSE, color_band_list = NULL,
                                        color_band_column=NULL, color_band_boundary=FALSE,
                                        alluvial_alpha = 0.5, match_colors = TRUE, output_plot_path = NULL,
                                        include_labels_in_boxes = FALSE, include_axis_titles = FALSE,
                                        include_group_sizes = FALSE
) {
    if (length(setdiff(c("id", "x", "y"), colnames(clus_df_gather))) > 0) {
        clus_df_gather <- ggforce::gather_set_data(clus_df_gather, 1:2)
    }

    if (!is.null(color_list)){
        ditto_colors <- color_list
    } else{
        ditto_colors <- default_colors
    }
    num_levels_group1 <- length(levels(clus_df_gather[[group1_name]]))
    num_levels_group2 <- length(levels(clus_df_gather[[group2_name]]))


    # Extract colors for each factor, assuming ditto_colors is long enough
    if (match_colors) {
        if (is.null(fixed_column) | fixed_column == 'col1_int') {
            colors_group1 <- ditto_colors[1:num_levels_group1]
            colors_group2 <- find_group2_colors(clus_df_gather, ditto_colors)

        } else {
            colors_group2 <- ditto_colors[1:num_levels_group2]
            colors_group1 <- find_group2_colors(clus_df_gather, ditto_colors,
                                                group1_name = 'col2_int', group2_name = 'col1_int')
        }
    } else {
        colors_group1 <- ditto_colors[1:num_levels_group1]
        colors_group2 <- ditto_colors[(num_levels_group1+1):(num_levels_group1+num_levels_group2)]
    }

    colors_group1_reverse <- rev(colors_group1)
    colors_group2_reverse <- rev(colors_group2)

    # Combine the colors
    combined_colors_reverse <- c(colors_group1_reverse, colors_group2_reverse)
    remaining_colors <- ditto_colors[!(ditto_colors %in% combined_colors_reverse)]

    temp_df <- clus_df_gather[1:as.integer(dim(clus_df_gather)[1]/2),1:dim(clus_df_gather)[2]]
    # uncomment to attempt mapping
    p <- ggplot(data = temp_df, aes(axis1 = !!sym('col1_int'),
                                    axis2 = !!sym('col2_int'), y = value),
    )

    if (color_bands) {
        if (!is.null(color_band_column)) {
            if (is.null(color_band_list)) {
                color_band_list <- remaining_colors
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
        } else if (fixed_column == 'col1_int') {
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
        } else {
            if (color_band_boundary){
                p <- p +
                    geom_alluvium(aes(fill = !!sym('col2_int'), color=!!sym('col2_int')), alpha = alluvial_alpha) +
                    scale_fill_manual(values = colors_group2) + scale_color_manual(values = colors_group2)+
                    labs(fill = NULL)
            } else{
                p <- p +
                    geom_alluvium(aes(fill = !!sym('col2_int')), alpha = alluvial_alpha) +
                    scale_fill_manual(values = colors_group2) +
                    labs(fill = NULL)
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
        p <- p + geom_stratum(fill = combined_colors_reverse,
        )
    } else {
        p <- p + geom_stratum()
    }

    if (!(include_labels_in_boxes==FALSE)) {
        A_label_names <- as.character(unique(clus_df_gather[order(clus_df_gather$col1_int),][[group1_name]]))
        B_label_names <- as.character(unique(clus_df_gather[order(clus_df_gather$col2_int),][[group2_name]]))

        final_label_names <- c(rev(A_label_names), rev(B_label_names))
        p <- p +
            geom_text(stat = StatStratum, aes(label = after_stat(final_label_names)))
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

    top_y <- max(top_y1/2, top_y2/2)  # top_y1 and top_y2 are probably the same

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

plot_alluvial_internal_multicol <- function(clus_df_gather,graphing_columns,
                                            sorting_algorithm = NULL,fixed_column=NULL,
                                            color_list = NULL, color_boxes = TRUE,
                                            color_bands = FALSE, color_band_list = NULL,
                                            color_band_column=NULL, color_band_boundary=FALSE,
                                            alluvial_alpha = 0.5, match_colors = TRUE, output_plot_path = NULL,
                                            include_labels_in_boxes = FALSE, include_axis_titles = FALSE,
                                            include_group_sizes = FALSE
) {
    if (length(setdiff(c("id", "x", "y"), colnames(clus_df_gather))) > 0) {
        clus_df_gather <- ggforce::gather_set_data(clus_df_gather, 1:2)
    }

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


add_int_columns <- function(df, graphing_columns) {
    n<-1
    for (col in graphing_columns) {
        col_int_name <- paste0('col', n, '_int')
        n <- n+1

        #factorize input columns
        df[[col]] <- factor(as.character(df[[col]]), levels = sort(unique(as.character(df[[col]])), method = "radix"))

        if (!(col_int_name %in% colnames(df))) {
            # make columns integer for sorting
            df[[col_int_name]] <- as.integer(factor(df[[col]], levels = sort(unique(df[[col]]), method = "radix")))
        }
    }
    return(df)
}

get_alluvial_df <- function(df, do_gather_set_data = FALSE) {
    # Convert numeric clustering columns to ordered factors
    df <- df |>
        dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
        dplyr::group_by_all() |>
        dplyr::count(name = "value")
    if (do_gather_set_data) {
        df <- ggforce::gather_set_data(df, 1:2)
    }
    return(df)
}

load_in_df <- function(df, graphing_columns = NULL, column_weights = NULL) {
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

    for (col in graphing_columns) {
        if (!(col %in% colnames(df))) {
            stop(sprintf("column '%s' is not a column in the dataframe.", col))
        }
        # convert to factor
        if (!is.factor(df[[col]])) {
            df[[col]] <- as.factor(df[[col]])
        }
    }

    return(df)
}

#' Preprocess data
#'
#' Preprocess data (load in, add integer columns, and group as needed)
#'
#' @param df A data frame, tibble, or CSV file path. Must contain at least two columns, each representing a clustering/grouping of the same entities (rows).
#' @param graphing_columns Optional Character. List of columns to use. Incompatible with column1 and column2. Optional if \code{df} has exactly two columns, or if \code{df} has exactly three columns including \code{column_weights}.
#' @param column_weights Optional numeric vector of weights (same length as number of rows in \code{df}) to weight each row differently when calculating flows.
#'
#' @return A data frame where each row provides the weight in the 'value' column that connects a given combination of values in graphing_columns.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- data_preprocess(df, graphing_columns = c('method1', 'method2'))
#' }
#'
#' @export
data_preprocess <- function(df, graphing_columns = NULL, column_weights = NULL, load_df = TRUE, do_gather_set_data = FALSE) {
    if (load_df) {
        df <- load_in_df(df = df, graphing_columns = graphing_columns, column_weights = column_weights)
    }
    df <- add_int_columns(df, graphing_columns = graphing_columns)

    if (is.null(column_weights) || !(column_weights %in% colnames(df))) {
        clus_df_gather <- get_alluvial_df(df, do_gather_set_data = do_gather_set_data)
    } else {
        clus_df_gather <- df
    }

    return(clus_df_gather)
}



sort_neighbornet <- function(clus_df_gather, graphing_columns = NULL, optimize_column_order = TRUE) {
    cycle <- run_neighbornet(clus_df_gather, graphing_columns=graphing_columns)
    res <- determine_optimal_cycle_start(clus_df_gather, cycle, graphing_columns=graphing_columns, optimize_column_order = optimize_column_order)
    clus_df_gather_neighbornet <- res$clus_df_gather
    return(clus_df_gather_neighbornet)
}

sort_greedy_wolf <- function(clus_df_gather, graphing_columns = NULL, column1 = NULL, column2 = NULL, fixed_column = 1, random_initializations = 1, set_seed = 42, column_weights = "value", sorting_algorithm = "greedy_WBLF") {
    if (length(graphing_columns) != 2) {
        stop(sprintf("graphing_columns must be of length 2 for greedy_wblf/greedy_wolf"))
    }

    if (is.null(fixed_column)) {
        fixed_column <- column1
    } else if ((is.integer(fixed_column) | (is.double(fixed_column)))) {
        if (fixed_column > length(colnames(clus_df_gather))){
            stop(sprintf("fixed_column index '%s' is not a column in the dataframe.", fixed_column))
        } else{
            fixed_column <- colnames(clus_df_gather)[fixed_column]
        }
    } else if (!(fixed_column %in% colnames(clus_df_gather))) {
        stop(sprintf("fixed_column '%s' is not a column in the dataframe.", fixed_column))
    }

    if (isTRUE(fixed_column == column1)) {
        fixed_column <- "col1_int"
        reordered_column <- "col2_int"
    } else if (isTRUE(fixed_column == column2)) {
        fixed_column <- "col2_int"
        reordered_column <- "col1_int"
    }
    crossing_edges_objective_minimum <- Inf
    set.seed(set_seed)


    if (length(setdiff(c("id", "x", "y"), colnames(clus_df_gather))) > 0) {
        clus_df_gather <- ggforce::gather_set_data(clus_df_gather, 1:2)
    }

    for (i in seq_len(random_initializations)) {
        if (sorting_algorithm == 'greedy_WBLF') {
            # randomize clus_df_gather order
            for (column_num in c('col1_int', 'col2_int')){
                clus_df_gather_tmp[[column_num]] = as.factor(clus_df_gather_tmp[[column_num]])
                clus_df_gather_tmp[[column_num]] = factor(clus_df_gather_tmp[[column_num]], levels=sample(levels(clus_df_gather_tmp[[column_num]])))
                # clus_df_gather_tmp[[column_num]] = as.integer(clus_df_gather_tmp[[column_num]])
            }
            # WBLF
            clus_df_gather_tmp <- sort_clusters_by_agreement(clus_df_gather_tmp, stable_column = 'col1_int',
                                                             reordered_column = 'col2_int')
            clus_df_gather_tmp <- sort_clusters_by_agreement(clus_df_gather_tmp, stable_column = 'col2_int',
                                                             reordered_column = 'col1_int')
        } else if (sorting_algorithm == 'greedy_WOLF') {
            # randomize clus_df_gather order
            column_num = reordered_column
            clus_df_gather_tmp[[column_num]] = as.factor(clus_df_gather_tmp[[column_num]])
            clus_df_gather_tmp[[column_num]] = factor(clus_df_gather_tmp[[column_num]], levels=sample(levels(clus_df_gather_tmp[[column_num]])))
            # clus_df_gather_tmp[[column_num]] = as.integer(clus_df_gather_tmp[[column_num]])
            # WOLF

            # browser()
            clus_df_gather_tmp <- sort_clusters_by_agreement(clus_df_gather_tmp, stable_column = fixed_column,
                                                             reordered_column = reordered_column)
        }
        if (random_initializations > 1) {
            crossing_edges_objective <- determine_crossing_edges(clus_df_gather_tmp, column1=column1, column2=column2,
                                                                 column_weights = column_weights, minimum_edge_weight = 0,
                                                                 output_df_path = NULL, return_weighted_layer_free_objective = TRUE)
            if (crossing_edges_objective < crossing_edges_objective_minimum) {
                crossing_edges_objective_minimum <- crossing_edges_objective
                clus_df_gather_best <- clus_df_gather_tmp
            }
        } else {
            clus_df_gather_best <- clus_df_gather_tmp
        }
    }

    if (length(setdiff(c("id", "x", "y"), colnames(clus_df_gather_best))) > 0) {
        clus_df_gather_best <- clus_df_gather_best %>%
            ungroup() %>%
            slice(1:(n() %/% 2)) %>%  # keep first half of rows
            select(-id, -x, -y)       # drop columns
    }

    return(clus_df_gather_best)
}

#' Sorts a dataframe (e.g., the output of data_preprocess)
#'
#' Maps the integers in integer columns corresponding to graphing_columns such that the alluvial plot will be sorted correctly.
#'
#' @param df A data frame, tibble, or CSV file path. Must contain at least two columns, each representing a clustering/grouping of the same entities (rows).
#'
#' @return A \code{ggplot2} object representing the alluvial plot.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' data_sort(df)
#' }
#'
#' @export
data_sort <- function(df, graphing_columns = NULL, column1 = NULL, column2 = NULL, column_weights = NULL, sorting_algorithm = "neighbornet", optimize_column_order = TRUE, fixed_column = 1, random_initializations = 1, set_seed = 42, output_df_path = NULL, load_df = TRUE) {
    #* Type Checking Start
    valid_algorithms <- c("neighbornet", "greedy_WOLF", "greedy_WBLF", "None")
    if (!(sorting_algorithm %in% valid_algorithms)) {
        stop(sprintf("Invalid sorting_algorithm: '%s'. Must be one of: %s",
                     sorting_algorithm, paste(valid_algorithms, collapse = ", ")))
    }

    if ((sorting_algorithm == 'None') & (random_initializations > 1)) {
        warning("random_initializations > 1 but sorting algorithm is None. Setting random_initializations to 1.")
        random_initializations=1
    }

    if ((is.null(fixed_column)) && (sorting_algorithm=='greedy_WOLF')) {
        stop(sprintf("Column to fix for One-Sided matching is not specified.", fixed_column))
    }

    if (!is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        stop("Specify either graphing_columns or column1/column2, not both.")
    }

    if (length(graphing_columns) == 2) {
        column1 <- graphing_columns[1]
        column2 <- graphing_columns[2]
    }

    if (is.null(graphing_columns)) {
        graphing_columns <- c(column1, column2)
    }
    #* Type Checking End

    # Preprocess (i.e., add int columns and do the grouping)
    clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns, column_weights = column_weights, load_df = load_df, do_gather_set_data = FALSE)

    if (sorting_algorithm == "neighbornet") {
        clus_df_gather_sorted <- sort_neighbornet(clus_df_gather = clus_df_gather, graphing_columns = graphing_columns, column_weights = column_weights, optimize_column_order = optimize_column_order)
    } else if (sorting_algorithm == "greedy_WBLF" || sorting_algorithm == "greedy_WOLF") {
        clus_df_gather_sorted <- sort_greedy_wolf(clus_df_gather = clus_df_gather, graphing_columns = graphing_columns, column1 = column1, column2 = column2, column_weights = column_weights, fixed_column = fixed_column, random_initializations = random_initializations, set_seed = set_seed, sorting_algorithm = sorting_algorithm)
    } else if (sorting_algorithm == "None") {
        clus_df_gather_sorted <- clus_df_gather
    } else {
        stop(sprintf("Invalid sorting_algorithm: '%s'. Must be one of: %s", sorting_algorithm, paste(valid_algorithms, collapse = ", ")))
    }

    # Save if desired
    if ((is.character(output_df_path) && grepl("\\.csv$", output_df_path, ignore.case = TRUE))) {
        write.csv(clus_df_gather_to_save, file = output_df_path, row.names = FALSE)
    }

    return(clus_df_gather_sorted)
}

#' Generate an Alluvial Plot with Minimal Cluster Cross-over
#'
#' Creates a two-axis alluvial plot to visualize the relationship between two categorical groupings (e.g., cluster assignments from different methods),
#' minimizing crossover and aligning clusters based on agreement.
#'
#' @param df A data frame, tibble, or CSV file path. Must contain at least two columns, each representing a clustering/grouping of the same entities (rows).
#' @param graphing_columns Optional Character. List of columns to use. Incompatible with column1 and column2. Optional if \code{df} has exactly two columns, or if \code{df} has exactly three columns including \code{column_weights}.
#' @param column1 Optional Character. Name of the first column to plot. Incompatible with graphing_columns. Optional if \code{df} has exactly two columns, or if \code{df} has exactly three columns including \code{column_weights}.
#' @param column2 Optional Character. Name of the second column to plot. Incompatible with graphing_columns. Optional if \code{df} has exactly two columns, or if \code{df} has exactly three columns including \code{column_weights}.
#' @param fixed_column Optional Character or Integer. Name of the column to fix, if desiring a one-layer free algorithm. If NULL, then implement both layers free.
#' @param random_initializations Optional Integer. Number of random initializations of the WLF heuristic to perform.
#' @param set_seed Optional Integer. Seed for random initializations of the WLF heuristic to perform.
#' @param color_band_column Optional Character. Which column to use for coloring bands.
#' @param color_band_boundary Logical. Whether or not to color boundaries between bands
#' @param sorting_algorithm Character. Must be neighbornet, greedy_WBLF, greedy_WOLF, or None.
#' @param optimize_column_order Logical. Only used for neighbornet.
#' @param color_boxes Logical. Whether to color the rectangular strata boxes representing groups.
#' @param color_bands Logical. Whether to color the alluvial bands connecting the groups.
#' @param match_colors Logical. If \code{TRUE}, assigns consistent colors between column1 and column2 where matched.
#' @param alluvial_alpha Numeric between 0 and 1. Transparency level for the alluvial bands.
#' @param include_labels_in_boxes Logical. Whether to include text labels inside the rectangular group boxes.
#' @param include_axis_titles Logical. Whether to display axis titles for column1 and column2.
#' @param include_group_sizes Logical. If \code{TRUE}, includes group sizes in the labels (e.g., "Group A (42)").
#' @param column_weights Optional numeric vector of weights (same length as number of rows in \code{df}) to weight each row differently when calculating flows.
#' @param output_plot_path Character. File path to save the plot (e.g., "plot.png"). If \code{NULL}, the plot is not saved.
#' @param output_df_path Character. File path to save the dataframe (e.g., "df.csv"). If \code{NULL}, the dataframe is not saved.
#' @param color_list Optional named list or vector of colors to override default group colors.
#' @param color_band_list Optional named list or vector of colors to override default band colors.
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
plot_alluvial <- function(df, graphing_columns = NULL, column1 = NULL, column2 = NULL, fixed_column = 1,
                          sorting_algorithm = 'neighbornet',random_initializations = 1, color_list = NULL,
                          color_boxes = TRUE,
                          color_bands = FALSE, color_band_list = NULL,
                          color_band_column=NULL, color_band_boundary=FALSE,
                          match_colors = TRUE, alluvial_alpha = 0.5,
                          include_labels_in_boxes = TRUE, include_axis_titles = TRUE, include_group_sizes = TRUE,
                          column_weights = NULL, output_plot_path = NULL, output_df_path = NULL,
                          set_seed=42, optimize_column_order=TRUE) {
    #* Type Checking Start
    if (!is.null(graphing_columns) && any(!graphing_columns %in% colnames(df))) {
        stop("Some graphing_columns are not present in the dataframe.")
    }

    # ensure someone doesn't specify both graphing_columns and column1/2
    if (!is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        stop("Specify either graphing_columns or column1/column2, not both.")
    }

    df <- load_in_df(df = df, graphing_columns = graphing_columns, column_weights = column_weights)

    if (ncol(df) < 2) {
        stop("Dataframe must have at least 2 columns when column_weights is NULL.")
    } else if (ncol(df) > 2) {
        if (is.null(graphing_columns) && is.null(column1) && is.null(column2)) {
            stop("graphing_columns must be specified when dataframe has more than 2 columns and column_weights is NULL.")
        }
    } else {  # length 2
        if (is.null(column1) && !is.null(column2)) {
            column1 <- setdiff(colnames(df), column2)
        } else if (is.null(column2) && !is.null(column1)) {
            column2 <- setdiff(colnames(df), column1)
        } else if (is.null(column1) && is.null(column2)) {
            column1 <- colnames(df)[1]
            column2 <- colnames(df)[2]
        }
    }

    # if someone specifies column1/2, then use it
    if (length(graphing_columns) == 2) {
        column1 <- graphing_columns[1]
        column2 <- graphing_columns[2]
    }

    if (is.null(graphing_columns)) {
        graphing_columns <- c(column1, column2)
    }

    if (is.null(fixed_column)) {
        fixed_column <- column1
    } else if ((is.integer(fixed_column) | (is.double(fixed_column)))) {
        if (fixed_column > length(colnames(df))){
            stop(sprintf("fixed_column index '%s' is not a column in the dataframe.", fixed_column))
        } else{
            fixed_column <- colnames(df)[fixed_column]
        }
    } else if (!(fixed_column %in% colnames(df))) {
        stop(sprintf("fixed_column '%s' is not a column in the dataframe.", fixed_column))
    }

    if (!is.null(color_band_column)) {
        color_bands <- TRUE
    }
    #* Type Checking End

    # Sort
    clus_df_gather <- data_sort(df = df, graphing_columns = graphing_columns, column_weights = column_weights, sorting_algorithm = sorting_algorithm, optimize_column_order = optimize_column_order, fixed_column = fixed_column, random_initializations = random_initializations, set_seed = set_seed, output_df_path = output_df_path, load_df = FALSE)

    # Plot
    # alluvial_plot <- plot_alluvial_internal()
    if (length(graphing_columns) == 2) {
        alluvial_plot <- plot_alluvial_internal_2col(clus_df_gather, group1_name = column1, group2_name = column2, fixed_column = fixed_column,
                                                     color_list = color_list, color_boxes = color_boxes,
                                                     color_bands = color_bands, color_band_list = color_band_list,
                                                     color_band_column=color_band_column, color_band_boundary=color_band_boundary,
                                                     match_colors = match_colors, alluvial_alpha = alluvial_alpha,
                                                     include_labels_in_boxes = include_labels_in_boxes, include_axis_titles = include_axis_titles,
                                                     include_group_sizes = include_group_sizes,
                                                     output_plot_path = output_plot_path)
    } else {
        alluvial_plot <- plot_alluvial_internal_multicol(clus_df_gather, graphing_columns=graphing_columns, fixed_column = fixed_column,
                                                         color_list = color_list, color_boxes = color_boxes,
                                                         color_bands = color_bands, color_band_list = color_band_list,
                                                         color_band_column=color_band_column, color_band_boundary=color_band_boundary,
                                                         match_colors = match_colors, alluvial_alpha = alluvial_alpha,
                                                         include_labels_in_boxes = include_labels_in_boxes, include_axis_titles = include_axis_titles,
                                                         include_group_sizes = include_group_sizes,
                                                         output_plot_path = output_plot_path)
    }

    return(alluvial_plot)
}


