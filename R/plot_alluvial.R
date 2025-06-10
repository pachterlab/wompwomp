#' alluvialmatch: Cluster-matching alluvial plots
#'
#' Main plotting function and helpers for bipartite-matching-based alluvial diagrams
#' @docType package
#' @name alluvialmatch
#'
#' @importFrom dplyr mutate select group_by summarise arrange desc ungroup slice n pull filter bind_rows across matches all_of
#' @importFrom ggplot2 ggplot aes geom_text scale_fill_manual labs after_stat annotate theme_void theme element_text rel ggsave guides scale_color_manual
#' @importFrom ggalluvial geom_alluvium geom_stratum stat_stratum stat_alluvium
#' @importFrom ggfittext geom_fit_text
#' @importFrom ggforce gather_set_data
#' @importFrom igraph max_bipartite_match V graph_from_data_frame
#' @importFrom tibble is_tibble
#' @importFrom utils read.csv write.csv combn
#' @importFrom stats setNames
#' @importFrom rlang sym .data
#' @importFrom magrittr %>%
#' @importFrom data.table :=

utils::globalVariables(c(
    ".data", ":=", "%>%", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "stratum", "total", "cum_y", "best_cluster_agreement", "neighbor_net", "alluvium", "pos", "count"
))

StatStratum <- ggalluvial::StatStratum  # avoid the error Can't find stat called "stratum"

neighbornet_script_path <- system.file("scripts", "run_neighbornet.py", package = "yourpackagename")
if (neighbornet_script_path == "") {
    # Fallback to development location
    neighbornet_script_path <- file.path(here::here("inst", "scripts", "run_neighbornet.py"))
    # neighbornet_script_path <- file.path(here::here("scripts", "run_neighbornet.py"))
}
stopifnot(file.exists(neighbornet_script_path))

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

determine_column_order <- function(clus_df_gather_neighbornet, graphing_columns, column_weights = "value", verbose = FALSE) {
    # this doesn't strictly need its own condition (2 choose 2 is 1 anyways), but does avoid a little overhead
    if (length(graphing_columns) == 2) {
        return(graphing_columns)
    }

    column_dist_matrix <- matrix(1e6, nrow = length(graphing_columns), ncol = length(graphing_columns),
                                 dimnames = list(graphing_columns, graphing_columns))

    pairs <- combn(graphing_columns, 2)
    for (i in 1:ncol(pairs)) {
        column1 <- pairs[1, i]
        column2 <- pairs[2, i]

        # Step 1: Get their positions
        integer1 <- which(names(clus_df_gather_neighbornet) == column1)
        integer2 <- which(names(clus_df_gather_neighbornet) == column2)

        # Step 2: Construct corresponding col{integer}_int names
        col1_int <- paste0("col", integer1, "_int")
        col2_int <- paste0("col", integer2, "_int")

        # Step 3: Reorder the data frame
        cols_to_keep <- c(column1, column2, col1_int, col2_int, column_weights)
        clus_df_gather_neighbornet_tmp <- clus_df_gather_neighbornet[, cols_to_keep]

        # Step 4: Rename columns col{integer1}_int → col1_int, col{integer2}_int → col2_int
        names(clus_df_gather_neighbornet_tmp)[match(c(col1_int, col2_int), names(clus_df_gather_neighbornet_tmp))] <- c("col1_int", "col2_int")
        graphing_columns_tmp <- c(column1, column2)

        neighbornet_objective <- determine_crossing_edges(
            clus_df_gather_neighbornet_tmp,
            graphing_columns = graphing_columns_tmp,
            column_weights = column_weights,
            return_weighted_layer_free_objective = TRUE,
            load_df = FALSE,
            preprocess_data = FALSE
        )
        neighbornet_objective <- log1p(neighbornet_objective)  # log1p to avoid issue of log(0)

        column_dist_matrix[column1, column2] <- neighbornet_objective
        column_dist_matrix[column2, column1] <- neighbornet_objective
    }
    # Prepare data in R
    labels <- graphing_columns  # assuming this is a character vector
    column_dist_matrix <- column_dist_matrix
    mat_list <- split(column_dist_matrix, row(column_dist_matrix))  # convert R matrix to list of row-vectors
    if (verbose) message("Running neighbornet for column order")
    reticulate::source_python(neighbornet_script_path)
    result <- neighbor_net(labels, column_dist_matrix)
    cycle <- result[[1]]
    cycle_mapped <- labels[cycle]

    # determine the optimal starting point for cycle
    adj_distances <- sapply(seq_len(length(cycle_mapped)), function(i) {
        from <- cycle_mapped[i]
        to <- cycle_mapped[(i %% length(cycle_mapped)) + 1]  # wraps around
        column_dist_matrix[from, to]
    })
    max_index <- which.max(adj_distances)

    cycle_mapped_optimal_start <- rotate_left(cycle_mapped, max_index)
    return(cycle_mapped_optimal_start)
}

run_neighbornet <- function(df, graphing_columns = NULL, column1 = NULL, column2 = NULL, column_weights = "value", verbose = FALSE) {
    # ensure someone doesn't specify both graphing_columns and column1/2
    if (!is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        stop("Specify either graphing_columns or column1/column2, not both.")
    }
    # if someone specifies column1/2, then use it
    if (is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        graphing_columns <- c(column1, column2)
    }

    # map from string to int if needed
    if (is.null(column_weights) || !(column_weights %in% colnames(df))) {
        clus_df_gather <- get_alluvial_df(df)
    } else {
        clus_df_gather <- df
    }

    # Add prefixes to distinguish node types

    # prefix is "tissue~~"
    for (col in graphing_columns) {
        clus_df_gather[[col]] <- paste0(col, "~~", clus_df_gather[[col]])
    }

    # prefix is "column1_"
    # for (i in seq_along(graphing_columns)) {
    #     col <- graphing_columns[i]
    #     clus_df_gather[[col]] <- paste0("column", i, "_", clus_df_gather[[col]])
    # }

    # Get all node names
    all_nodes <- sort(unique(unlist(clus_df_gather[graphing_columns])))

    # Compute full distance matrix based on -log(edge weight)
    # Initialize distance matrix
    full_dist_matrix <- matrix(1e6, nrow = length(all_nodes), ncol = length(all_nodes),
                               dimnames = list(all_nodes, all_nodes))

    # Get all 2-column combinations
    pairwise_groupings <- combn(graphing_columns, 2, simplify = FALSE)

    # For each combination, group and summarize
    summarized_results <- purrr::map(pairwise_groupings, function(cols) {
        clus_df_gather %>%
            group_by(across(all_of(cols))) %>%
            summarise(total_value = sum(value), .groups = "drop") %>%
            mutate(grouping = paste(cols, collapse = "+"))
    })

    # Combine into a single data frame
    final_result <- bind_rows(summarized_results)

    for (i in seq_len(nrow(final_result))) {
        grouping_str <- final_result$grouping[i]
        parts <- strsplit(grouping_str, "\\+")[[1]]
        column1_tmp <- parts[1]
        column2_tmp <- parts[2]

        n1 <- as.character(final_result[[column1_tmp]][i])
        n2 <- as.character(final_result[[column2_tmp]][i])
        w <- final_result$total_value[i]

        if (w > 0) {
            full_dist_matrix[n1, n2] <- -log(w)
            full_dist_matrix[n2, n1] <- -log(w)  # symmetric since graph is undirected
        }
    }

    # make sure all numbers are positive for neighbornet
    min_val_abs <- abs(min(full_dist_matrix))
    full_dist_matrix <- full_dist_matrix + (min_val_abs + 1)

    # Prepare data in R
    labels <- all_nodes  # assuming this is a character vector
    mat <- full_dist_matrix
    mat[is.infinite(mat)] <- 1e6
    mat[is.na(mat)] <- 1e6
    mat_list <- split(mat, row(mat))  # convert R matrix to list of row-vectors

    # Call Python function
    # result <- nn_mod$neighbor_net(labels, mat_list)
    if (verbose) message("Running neighbornet for stratum order")
    reticulate::source_python(neighbornet_script_path)
    result <- neighbor_net(labels, mat)
    cycle <- result[[1]]
    splits <- result[[2]]

    cycle_mapped <- labels[cycle]

    return(cycle_mapped)
}

rotate_left <- function(vec, k = 1) {
    n <- length(vec)
    k <- k %% n
    if (k == 0) return(vec)
    c(vec[(k + 1):n], vec[1:k])
}

get_graph_groups <- function(cycle) {
    groups <- list()

    for (node in cycle) {
        prefix <- sub("~~.*", "", node)  # Extract everything before the first `~~`

        if (!prefix %in% names(groups)) {
            groups[[prefix]] <- c()
        }

        groups[[prefix]] <- c(groups[[prefix]], node)
    }

    return(groups)
}

# swap_columns_in_clus_df_gather <- function(clus_df_gather, graphing_columns_int) {
#   # Find all *_int columns in the dataframe
#   original_int_cols <- grep("^col[0-9]+_int$", names(clus_df_gather), value = TRUE)
#
#   # Create new names in the desired order
#   new_names <- paste0("col", seq_along(graphing_columns_int), "_int")
#
#   # Map from current column name -> new column name
#   names_map <- setNames(new_names, graphing_columns_int)
#
#   # Rename the columns accordingly
#   matched_cols <- names(clus_df_gather) %in% names_map
#   names(clus_df_gather)[matched_cols] <- names_map[names(clus_df_gather)[matched_cols]]
#
#   return(clus_df_gather)
# }

# swap_graphing_column_order_based_on_graphing_column_int_order <- function(graphing_columns, graphing_columns_int) {
#     # Get the index of each graphing_columns_int entry (e.g., "col2_int" → 2)
#     int_positions <- as.integer(gsub("col([0-9]+)_int", "\\1", graphing_columns_int))
#
#     # Create an empty character vector of the correct length
#     reordered_graphing_columns <- character(length(graphing_columns))
#
#     # Place each graphing column at its new position
#     reordered_graphing_columns[int_positions] <- graphing_columns
#
#     return(reordered_graphing_columns)
# }

# # example:
# graphing_columns          <- c("tissue", "sex", "cluster")
# graphing_columns_int      <- c("col2_int", "col3_int", "col1_int")
# output: c("sex", "cluster", "tissue")
swap_graphing_column_order_based_on_graphing_column_int_order <- function(graphing_columns, graphing_columns_int) {
    # Extract suffixes to determine new order
    suffixes <- as.integer(gsub("col([0-9]+)_int", "\\1", graphing_columns_int))

    # Match suffix to index in original graphing_columns (which are col1 = graphing_columns[1], etc.)
    return(graphing_columns[suffixes])
}



determine_optimal_cycle_start <- function(df, cycle, graphing_columns = NULL, column1 = NULL, column2 = NULL, column_weights = "value", optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE, verbose = FALSE) {
    # ensure someone doesn't specify both graphing_columns and column1/2
    if (!is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        stop("Specify either graphing_columns or column1/column2, not both.")
    }
    # if someone specifies column1/2, then use it
    if (is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        graphing_columns <- c(column1, column2)
    }

    #factorize input columns
    for (col in graphing_columns) {
        df[[col]] <- as.factor(as.character(df[[col]]))
    }

    neighbornet_objective_minimum <- Inf
    p_best_neighbornet <- NULL
    cycle_best <- NULL
    clus_df_gather_best <- NULL
    graphing_columns_best <- NULL
    objective_matrix_vector <- c()

    n <- length(cycle)
    for (i in 0:(n - 1)) {
        if (i == 0) {
            if (verbose) message(sprintf("Starting iteration 1"))
        } else if (i == 1) {
            if (verbose) message(sprintf("Starting subsequent iterations (should go much faster than iteration 1)"))
        }

        cycle_shifted <- rotate_left(cycle, i)
        graphs_list <- get_graph_groups(cycle_shifted)

        # remove prefix (column1_, etc)
        graphs_list_stripped <- lapply(graphs_list, function(x) {
            sub("^.*?~~", "", x)
        })

        if (is.null(column_weights) || !(column_weights %in% colnames(df))) {
            clus_df_gather_neighbornet <- get_alluvial_df(df)
        } else {
            clus_df_gather_neighbornet <- df
        }

        graphing_columns_int <- c()
        for (j in seq_along(graphing_columns)) {
            col_name <- graphing_columns[j]
            int_col_name <- paste0("col", j, "_int")
            graph <- graphs_list_stripped[[col_name]]

            # Assign the new integer-mapped column
            clus_df_gather_neighbornet[[int_col_name]] <- match(clus_df_gather_neighbornet[[col_name]], graph)

            # Collect the new column name
            graphing_columns_int <- c(graphing_columns_int, int_col_name)
        }

        if (optimize_column_order) {
            # optimize order either on the first iteration if optimize_column_order_per_cycle is FALSE, or each time if optimize_column_order_per_cycle is TRUE
            if ((optimize_column_order_per_cycle) || (i == 0)) {
                verbose_tmp <- if (i == 0) verbose else FALSE  # only have the option for verbose on first iteration
                graphing_columns_tmp <- determine_column_order(clus_df_gather_neighbornet, graphing_columns = graphing_columns, column_weights = column_weights, verbose = verbose_tmp)
            } else {
                graphing_columns_tmp <- graphing_columns
            }
        } else {
            graphing_columns_tmp <- graphing_columns
        }

        # graphing_columns_tmp <- swap_graphing_column_order_based_on_graphing_column_int_order(graphing_columns=graphing_columns, graphing_columns_int=graphing_columns_int_sorted)
        clus_df_gather_neighbornet_tmp <- reorder_and_rename_columns(clus_df_gather_neighbornet, graphing_columns_tmp)

        if ((optimize_column_order_per_cycle) || (i == 0)) {
            neighbornet_objective_output <- determine_crossing_edges(
                clus_df_gather_neighbornet_tmp,
                graphing_columns = graphing_columns_tmp,
                include_output_objective_matrix_vector = TRUE,
                return_weighted_layer_free_objective = FALSE
            )
        } else {
            swapped_node <- cycle_shifted[length(cycle_shifted)]
            parts <- strsplit(swapped_node, "~~", fixed = TRUE)[[1]]
            layer_name <- parts[1]
            node_name <- parts[2]

            # determine the int column that matches to this layer_name - because I haven't reordered any columns in clus_df_gather_neighbornet, I should use the order as determined here
            int_column_int <- match(layer_name, graphing_columns_tmp)
            int_column <- paste0("col", int_column_int, "_int")

            # find the value in clus_df_gather_neighbornet[[int_column]] that maps to node_name of clus_df_gather_neighbornet[[layer_name]]
            matched <- clus_df_gather_neighbornet_tmp[[int_column]][clus_df_gather_neighbornet_tmp[[layer_name]] == node_name]
            stratum_int_name <- matched[1]
            stratum_column_and_value_to_keep <- setNames(list(stratum_int_name), as.character(int_column_int))

            neighbornet_objective_output <- determine_crossing_edges(
                clus_df_gather_neighbornet_tmp,
                graphing_columns = graphing_columns_tmp,
                stratum_column_and_value_to_keep = stratum_column_and_value_to_keep,
                input_objective = neighbornet_objective,
                input_objective_matrix_vector = objective_matrix_vector,
                include_output_objective_matrix_vector = TRUE,
                return_weighted_layer_free_objective = FALSE
            )
        }
        objective_matrix_vector <- neighbornet_objective_output$objective_matrix_vector
        neighbornet_objective <- neighbornet_objective_output$output_objective

        # print(neighbornet_objective)
        if (verbose) message(sprintf("neighbornet_objective for iteration %s = %s", i+1, neighbornet_objective))
        if (neighbornet_objective < neighbornet_objective_minimum) {

            neighbornet_objective_minimum <- neighbornet_objective
            cycle_best <- cycle_shifted
            individual_graphs <- graphs_list_stripped
            # p_best_neighbornet <- p_neighbornet
            graphing_columns_best <- graphing_columns_tmp
            clus_df_gather_best <- clus_df_gather_neighbornet_tmp
        }
    }

    # clus_df_gather_best <- reorder_and_rename_columns(clus_df_gather_best, graphing_columns_best)  # done earlier now

    # make factors
    for (j in seq_along(graphing_columns_best)) {
        int_col_name <- paste0("col", j, "_int")
        clus_df_gather_best[[int_col_name]] <- factor(clus_df_gather_best[[int_col_name]])
    }

    return(list(cycle = cycle_best, individual_graphs = individual_graphs, neighbornet_objective = neighbornet_objective_minimum, clus_df_gather = clus_df_gather_best, graphing_columns = graphing_columns_best))
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


find_group2_colors <- function(clus_df_gather, ditto_colors,unused_colors, current_g1_colors,
                               group1_name = 'col1_int', group2_name = 'col2_int', group2_colors=NULL) {
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
    if (is.null(group2_colors)) {
        group2_colors <- vector("character", number_group2_clusters)
    }

    # Assign colors based on the matching
    for (i in g1_indices) {
        g2_indice <-  as.integer(sub("G2_", "", g1_to_g2_pairs[paste0('G1_',i)]))
        group2_colors[g2_indice] <- current_g1_colors[names(current_g1_colors) == i]
    }

    remaining_colors <- unused_colors

    num_remaining_group2 <- length(which(group2_colors %in% c('')))

    if (num_remaining_group2==0){
        return (group2_colors)
    } else{
        if ((length(remaining_colors)==0) | (length(remaining_colors) < num_remaining_group2)){
            smaller_clus_df_filtered <- clus_df_gather[!(clus_df_filtered[[group2_name]] %in% g1_to_g2_pairs),]
            smaller_clus_df_filtered[[group2_name]] <- as.integer(droplevels(as.factor(smaller_clus_df_filtered[[group2_name]])))
            sub_group2_colors <- find_group2_colors(smaller_clus_df_filtered, ditto_colors,remaining_colors,current_g1_colors,
                                                    group1_name, group2_name)
            for (i in which(group2_colors %in% c(''))) {
                group2_colors[i] <- sub_group2_colors[1]
                sub_group2_colors <- sub_group2_colors[2:length(sub_group2_colors)]
            }
        } else {
            for (i in which(group2_colors %in% c(''))) {
                group2_colors[i] <- remaining_colors[1]
                remaining_colors <- remaining_colors[2:length(remaining_colors)]
            }
            }
        return(group2_colors)
    }
}

plot_alluvial_internal <- function(clus_df_gather,graphing_columns,
                                            color_list = NULL, color_boxes = TRUE,
                                            color_bands = FALSE, color_band_list = NULL,
                                            color_band_column=NULL, color_band_boundary=FALSE,
                                            alluvial_alpha = 0.5, match_colors = TRUE, output_plot_path = NULL,
                                            include_labels_in_boxes = FALSE, include_axis_titles = FALSE,
                                            include_group_sizes = FALSE, verbose = FALSE,
                                   box_width = 1/3, text_width = 1/4, min_text = 4, 
                                   save_height = 6, save_width = 6
) {
    clus_df_gather <- ggforce::gather_set_data(clus_df_gather, 1:2)
    clus_df_gather <- clus_df_gather[clus_df_gather$x == 1,]

    if (!is.null(color_list)){
        ditto_colors <- color_list
    } else{
        ditto_colors <- default_colors
    }

    # Extract colors for each factor, assuming ditto_colors is long enough
    if (match_colors) {
        unused_colors <- ditto_colors
        first <- TRUE
        final_colors <- c()
        n <- 1
        for (col_group in graphing_columns) {
            num_levels <- length(levels(clus_df_gather[[col_group]]))
            if (first) {
                old_colors <- unused_colors[1:num_levels]
                names(old_colors) <- 1:num_levels
                unused_colors <- unused_colors[!(unused_colors %in% old_colors)]
                final_colors <- c(final_colors, rev(old_colors))
                first <- FALSE
            } else {
                temp_colors <- find_group2_colors(clus_df_gather, ditto_colors,unused_colors,old_colors,
                                                  group1_name = paste0('col',n-1,'_int'), group2_name = paste0('col',n,'_int'))
                unused_colors <- unused_colors[!(unused_colors %in% old_colors)]
                old_colors <- temp_colors
                names(old_colors) <- 1:num_levels
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
    temp_df <- clus_df_gather#[1:as.integer(dim(clus_df_gather)[1]/2),1:dim(clus_df_gather)[2]]

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
        p <- p + geom_stratum(width = box_width, fill = final_colors)
    } else {
        p <- p + geom_stratum(width = box_width)
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
            #geom_text(stat = StatStratum, aes(label = after_stat(final_label_names)))+
            geom_fit_text(reflow = TRUE,stat = "stratum", width = text_width, min.size = min_text, 
                          aes(label = after_stat(final_label_names)))
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
        top_y <- max(curr_y, top_y)
    }# top_y1 and top_y2 are probably the same

    if (include_axis_titles) {
        # Offset to place labels a bit above
        text_size <- 5
        #x<-1
        #for (col_group in graphing_columns) {
        p <- p +
            annotate("text", x = seq(1, length(graphing_columns)), y = rep(2*text_size+top_y, length(graphing_columns)), 
                     label = graphing_columns, size = text_size, hjust = 0.5)
        #    x <- x+1
        #}
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
        if (verbose) message(sprintf("Saving plot to=%s", output_plot_path))
        ggsave(output_plot_path, plot = p, 
               height = save_height, width = save_width, dpi = 300, bg = "white")
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

    if (!(is.null(graphing_columns))) {
        for (col in graphing_columns) {
            if (!(col %in% colnames(df))) {
                stop(sprintf("column '%s' is not a column in the dataframe.", col))
            }
            # convert to factor
            if (!is.factor(df[[col]])) {
                df[[col]] <- as.factor(df[[col]])
            }
        }
    }

    return(df)
}

# reorders int columns to match graphing_columns - eg if df has tissue, cluster, col1_int (for tissue), col2_int (for cluster) and graphing_columns = (cluster, tissue), then the output df will have cluster, tissue, col1_int (for cluster), col2_int (for tissue)
reorder_and_rename_columns <- function(df, graphing_columns) {
    # Find the order in df of the columns listed in graphing_columns
    original_graphing_columns <- intersect(colnames(df), graphing_columns)

    # Get original colX_int names based on original order
    original_int_cols <- paste0("col", seq_along(original_graphing_columns), "_int")

    # Target int column names based on desired new graphing order
    new_int_cols <- paste0("col", seq_along(graphing_columns), "_int")

    # Fix: old col name → new col name
    old_int_cols <- original_int_cols[match(graphing_columns, original_graphing_columns)]
    old_to_new_int_names <- setNames(new_int_cols, old_int_cols)

    # Rename df
    names(df)[names(df) %in% names(old_to_new_int_names)] <-
        old_to_new_int_names[names(df)[names(df) %in% names(old_to_new_int_names)]]

    # Final column order
    everything_else <- setdiff(names(df), c(graphing_columns, new_int_cols))
    df <- df[, c(graphing_columns, new_int_cols, everything_else)]

    return(df)
}






#' Preprocess data
#'
#' Preprocess data (load in, add integer columns, reorder columns to match graphing_columns, and group as needed)
#'
#' @param df A data frame, tibble, or CSV file path. Must be in one of two formats:
#' (1) column_weights == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two graphing_columns).
#' (2) column_weights != NULL: Each row represents a combination of groupings, each column from \code{graphing_columns} represents a grouping, and the column \code{column_weights} represents the number of entities in that combination of groupings. Must contain at least three columns (two \code{graphing_columns}, one \code{column_weights}).
#' @param graphing_columns Character vector. Vector of column names from \code{df} to be used in graphing (i.e., alluvial plotting).
#' @param column_weights Optional character. Column name from \code{df} that contains the weights of each combination of groupings if \code{df} is in format (2) (see above).
#' @param output_df_path Optional character. Output path for the output data frame, in CSV format. If \code{NULL}, then will not be saved.
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param load_df Internal flag; not recommended to modify.
#' @param do_gather_set_data Internal flag; not recommended to modify.
#'
#' @return A data frame where each row represents a combination of groupings, each column from \code{graphing_columns} represents a grouping, and the column \code{column_weights} ('value' if \code{column_weights} == NULL) represents the number of entities in that combination of groupings. For each column in \code{graphing_columns}, there will be an additional column \code{col1_int}, \code{col2_int}, etc. where each column corresponds to a position mapping of groupings in the respective entry of \code{graphing_columns} - for example, \code{col1_int} corresponds to \code{graphing_columns[1]}, \code{col2_int} corresponds to \code{graphing_columns[2]}, etc.
#'
#' @examples
#' # Example 1: df format 1
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- data_preprocess(df, graphing_columns = c('method1', 'method2'))
#'
#' # Example 2: df format 2
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- df |>
#'     dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
#'     dplyr::group_by_all() |>
#'     dplyr::count(name = "value")
#' clus_df_gather <- data_preprocess(clus_df_gather, graphing_columns = c('method1', 'method2'), column_weights = 'value')
#'
#' @export
data_preprocess <- function(df, graphing_columns, column_weights = NULL, output_df_path = NULL, verbose = FALSE, load_df = TRUE, do_gather_set_data = FALSE) {
    if (load_df) {
        df <- load_in_df(df = df, graphing_columns = graphing_columns, column_weights = column_weights)
    }

    # Fill in NA values with "Missing"
    df[is.na(df)] <- "Missing"

    df <- add_int_columns(df, graphing_columns = graphing_columns)

    # sort columns according to graphing_columns
    # df <- df %>% dplyr::relocate(all_of(graphing_columns))  # put graphing_columns in front
    if (!all(intersect(colnames(df), graphing_columns) == graphing_columns)) {
        df <- reorder_and_rename_columns(df, graphing_columns)
    }

    if (is.null(column_weights) || !(column_weights %in% colnames(df))) {
        clus_df_gather <- get_alluvial_df(df, do_gather_set_data = do_gather_set_data)
    } else {
        clus_df_gather <- df
    }

    if ((is.character(output_df_path) && grepl("\\.csv$", output_df_path, ignore.case = TRUE))) {
        if (verbose) message(sprintf("Saving sorted dataframe to=%s", output_df_path))
        write.csv(clus_df_gather, file = output_df_path, row.names = FALSE)
    }

    return(clus_df_gather)
}



sort_neighbornet <- function(clus_df_gather, graphing_columns = NULL, column_weights = "value", optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE, verbose = FALSE) {
    if (!reticulate::py_module_available("splitspy")) {
        stop("Python module 'splitspy' is not available, which is required for the neighbornet algorithm (default). Please run alluvialmatch::setup_python_env().")
    }
    if (verbose) message("Running neighbornet")
    cycle <- run_neighbornet(clus_df_gather, graphing_columns=graphing_columns, column_weights=column_weights, verbose=verbose)
    if (verbose) message("Determining optimal cycle start")
    res <- determine_optimal_cycle_start(clus_df_gather, cycle, graphing_columns=graphing_columns, column_weights=column_weights, optimize_column_order = optimize_column_order, optimize_column_order_per_cycle=optimize_column_order_per_cycle, verbose=verbose)
    clus_df_gather_neighbornet <- res$clus_df_gather
    # graphing_columns_neighbornet <- res$graphing_columns
    if (verbose) message(sprintf("crossing edges objective = %s", res$neighbornet_objective))
    return(clus_df_gather_neighbornet)
}

sort_greedy_wolf <- function(clus_df_gather, graphing_columns = NULL, column1 = NULL, column2 = NULL, fixed_column = NULL, random_initializations = 1, set_seed = 42, column_weights = "value", sorting_algorithm = "greedy_WBLF", verbose = FALSE) {
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
    } else {
        stop(sprintf("fixed_column '%s' is not recognized.", fixed_column))
    }
    crossing_edges_objective_minimum <- Inf
    set.seed(set_seed)


    length_clus_df_gather_original <- nrow(clus_df_gather)
    if (length(setdiff(c("id", "x", "y"), colnames(clus_df_gather))) > 0) {
        clus_df_gather <- ggforce::gather_set_data(clus_df_gather, 1:2)
    }

    for (i in seq_len(random_initializations)) {
        clus_df_gather_tmp <- clus_df_gather
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

            clus_df_gather_tmp <- sort_clusters_by_agreement(clus_df_gather_tmp, stable_column = fixed_column,
                                                             reordered_column = reordered_column)
        }
        if (random_initializations > 1) {
            crossing_edges_objective <- determine_crossing_edges(clus_df_gather_tmp, column1=column1, column2=column2,
                                                                 column_weights = column_weights, load_df = FALSE, preprocess_data = FALSE,
                                                                 output_df_path = NULL, return_weighted_layer_free_objective = TRUE)
            if (crossing_edges_objective < crossing_edges_objective_minimum) {
                crossing_edges_objective_minimum <- crossing_edges_objective
                clus_df_gather_best <- clus_df_gather_tmp
            }
            if (verbose) message(sprintf("Complete with iteration=%d", i))

        } else {
            clus_df_gather_best <- clus_df_gather_tmp
        }
    }

    # checks if id, x, y in clus_df_gather_best, as well as if clus_df_gather_best has more rows than original
    if ((length(setdiff(c("id", "x", "y"), colnames(clus_df_gather_best))) == 0) && nrow(clus_df_gather_best) > length_clus_df_gather_original) {
        clus_df_gather_best <- clus_df_gather_best %>%
            ungroup() %>%
            slice(1:(n() %/% 2)) %>%  # keep first half of rows
            select(-id, -x, -y)       # drop columns
    }

    return(clus_df_gather_best)
}


#' Sorts a dataframe.
#'
#' Sorts a dataframe with the algorithm specified by \code{sorting_algorithm}.
#'
#' @param df A data frame, tibble, or CSV file path. Must be in one of two formats:
#' (1) column_weights == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two graphing_columns).
#' (2) column_weights != NULL: Each row represents a combination of groupings, each column from \code{graphing_columns} represents a grouping, and the column \code{column_weights} represents the number of entities in that combination of groupings. Must contain at least three columns (two \code{graphing_columns}, one \code{column_weights}).
#' @param graphing_columns Optional character vector. Vector of column names from \code{df} to be used in graphing (i.e., alluvial plotting). Mutually exclusive with \code{column1} and \code{column2}.
#' @param column1 Optional character. Can be used along with \code{column2} in place of \code{graphing_columns} if working with two columns only. Mutually exclusive with \code{graphing_columns}.
#' @param column2 Optional character. Can be used along with \code{column1} in place of \code{graphing_columns} if working with two columns only. Mutually exclusive with \code{graphing_columns}.
#' @param column_weights Optional character. Column name from \code{df} that contains the weights of each combination of groupings if \code{df} is in format (2) (see above).
#' @param sorting_algorithm Character. Algorithm with which to sort the values in the dataframe. Can choose from {'neighbornet', 'greedy_WOLF', 'greedy_WBLF', 'None'}. 'neighbornet' performs sorting with NeighborNet (Bryant and Moulton, 2004). 'greedy_WOLF' implements a custom greedy algorithm where one layer is fixed, and the other layer is sorted such that each node is positioned as close to its largest parent from the fixed side as possible in a greedy fashion. 'greedy_WBLF' implements the 'greedy_WOLF' algorithm described previously twice, treating each column as fixed in one iteration and free in the other iteration. 'greedy_WOLF' and 'greedy_WBLF' are only valid when \code{graphing_columns} has exactly two entries.
#' @param optimize_column_order Logical. If TRUE, will optimize the order of \code{graphing_columns} to minimize edge overlap. Only applies when \code{sorting_algorithm == 'neighbornet'} and \code{length(graphing_columns) > 2}.
#' @param optimize_column_order_per_cycle Logical. If TRUE, will optimize the order of \code{graphing_columns} to minimize edge overlap upon each cycle. If FALSE, will optimize the order of \code{graphing_columns} to minimize edge overlap on the beginning cycle only. Only applies when \code{sorting_algorithm == 'neighbornet'} and \code{length(graphing_columns) > 2}.
#' @param fixed_column Character or Integer. Name or position of the column in \code{graphing_columns} to keep fixed during sorting. Only applies when \code{sorting_algorithm == 'greedy_WOLF'}.
#' @param random_initializations Integer. Number of random initializations for the positions of each grouping in \code{graphing_columns}. Only applies when \code{sorting_algorithm == 'greedy_WOLF' or sorting_algorithm == 'greedy_WBLF'}.
#' @param set_seed Integer. Random seed for the \code{random_initializations} parameter. Only applies when \code{sorting_algorithm == 'greedy_WOLF' or sorting_algorithm == 'greedy_WBLF'}.
#' @param output_df_path Optional character. Output path for the output data frame, in CSV format. If \code{NULL}, then will not be saved.
#' @param preprocess_data Logical. If TRUE, will preprocess the data with the \code{data_preprocess} function.
#' @param return_updated_graphing_columns Logical. If FALSE, will only return the updated data frame. If TRUE, will return both the updated data frame and the updated graphing_columns parameter in the order in which the columns should be graphed.
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param load_df Internal flag; not recommended to modify.
#'
#' @return
#' If return_updated_graphing_columns == FALSE (default): A data frame where each row represents a combination of groupings, each column from \code{graphing_columns} represents a grouping, and the column \code{column_weights} ('value' if \code{column_weights} == NULL) represents the number of entities in that combination of groupings. For each column in \code{graphing_columns}, there will be an additional column \code{col1_int}, \code{col2_int}, etc. where each column corresponds to a position mapping of groupings in the respective entry of \code{graphing_columns} - for example, \code{col1_int} corresponds to \code{graphing_columns[1]}, \code{col2_int} corresponds to \code{graphing_columns[2]}, etc. The position mappings in these columns, as well as the order of the columns (if \code{optimize_column_order} is TRUE), will be sorted according to \code{sorting_algorithm}.
#' If return_updated_graphing_columns == TRUE: A list of the data frame described above and the sorted \code{graphing_columns}, in the keys 'clus_df_gather' and 'graphing_columns', respectively.
#'
#' @examples
#' # Example 1: df format 1
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- data_sort(df, graphing_columns = c('method1', 'method2'))
#'
#' # Example 2: df format 2
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- df |>
#'     dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
#'     dplyr::group_by_all() |>
#'     dplyr::count(name = "value")
#' clus_df_gather <- data_sort(clus_df_gather, graphing_columns = c('method1', 'method2'), column_weights = 'value')
#'
#' @export
data_sort <- function(df, graphing_columns = NULL, column1 = NULL, column2 = NULL, column_weights = NULL, sorting_algorithm = "neighbornet", optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE, fixed_column = NULL, random_initializations = 1, set_seed = 42, output_df_path = NULL, preprocess_data = TRUE, return_updated_graphing_columns = FALSE, verbose = FALSE, load_df = TRUE) {
    #* Type Checking Start
    valid_algorithms <- c("neighbornet", "greedy_WOLF", "greedy_WBLF", "None")
    if (!(sorting_algorithm %in% valid_algorithms)) {
        stop(sprintf("Invalid sorting_algorithm: '%s'. Must be one of: %s",
                     sorting_algorithm, paste(valid_algorithms, collapse = ", ")))
    }

    if (sorting_algorithm == 'neighbornet') {
        for (col in graphing_columns) {
            if (grepl("~~", col)) {
                stop(sprintf("No entry of graphing_columns can contain '~~' when sorting_algorithm == neighbornet. Issue with column '%s'.", col))
            }
        }
    }

    if ((sorting_algorithm == 'None') & (random_initializations > 1)) {
        warning("random_initializations > 1 but sorting algorithm is None. Setting random_initializations to 1.")
        random_initializations=1
    }

    if (!is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        stop("Specify either graphing_columns or column1/column2, not both.")
    }

    if (!is.null(graphing_columns) && any(!graphing_columns %in% colnames(df))) {
        stop("Some graphing_columns are not present in the dataframe.")
    }

    if ((is.character(df)) && (load_df)) {
        if (verbose) message("Loading in data")
        df <- load_in_df(df = df, graphing_columns = graphing_columns, column_weights = column_weights)
        load_df <- FALSE
    }

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

    if (length(graphing_columns) == 2) {
        column1 <- graphing_columns[1]
        column2 <- graphing_columns[2]
    }

    if (is.null(graphing_columns)) {
        graphing_columns <- c(column1, column2)
    }

    if ((is.null(fixed_column)) && (sorting_algorithm=='greedy_WOLF')) {
        if (verbose) message(sprintf("Using column %s as fixed_column for greedy_WOLF by default", column1))
        fixed_column <- column1
    }
    #* Type Checking End


    # Preprocess (i.e., add int columns and do the grouping)
    if (preprocess_data) {
        if (verbose) message("Preprocessing data before sorting")
        clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns, column_weights = column_weights, load_df = load_df, do_gather_set_data = FALSE)
        if (is.null(column_weights)) {
            column_weights <- "value"  # is set during data_preprocess
        }
    } else {
        clus_df_gather <- df
    }

    if (sorting_algorithm == "neighbornet") {
        clus_df_gather_sorted <- sort_neighbornet(clus_df_gather = clus_df_gather, graphing_columns = graphing_columns, column_weights = column_weights, optimize_column_order = optimize_column_order, optimize_column_order_per_cycle = optimize_column_order_per_cycle, verbose = verbose)
    } else if (sorting_algorithm == "greedy_WBLF" || sorting_algorithm == "greedy_WOLF") {
        clus_df_gather_sorted <- sort_greedy_wolf(clus_df_gather = clus_df_gather, graphing_columns = graphing_columns, column1 = column1, column2 = column2, column_weights = column_weights, fixed_column = fixed_column, random_initializations = random_initializations, set_seed = set_seed, sorting_algorithm = sorting_algorithm, verbose = verbose)
    } else if (sorting_algorithm == "None") {
        clus_df_gather_sorted <- clus_df_gather
    } else {
        stop(sprintf("Invalid sorting_algorithm: '%s'. Must be one of: %s", sorting_algorithm, paste(valid_algorithms, collapse = ", ")))
    }

    # print objective - don't do for neighbornet because I did it right before
    if ((verbose) && (sorting_algorithm != "neighbornet")) {
        objective <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns=graphing_columns,
                                 column_weights = column_weights, load_df = FALSE, preprocess_data = TRUE, return_weighted_layer_free_objective = TRUE)
        message(sprintf("crossing edges objective = %s", objective))
    }

    if (verbose) message("Complete with sorting")

    # Save if desired
    if ((is.character(output_df_path) && grepl("\\.csv$", output_df_path, ignore.case = TRUE))) {
        if (verbose) message(sprintf("Saving sorted dataframe to=%s", output_df_path))
        write.csv(clus_df_gather_sorted, file = output_df_path, row.names = FALSE)
    }

    if (return_updated_graphing_columns) {
        graphing_columns_sorted <- graphing_columns[order(match(graphing_columns, names(clus_df_gather_sorted)))]  # reorder graphing_columns to match any changed order in clus_df_gather_sorted
        return(list(clus_df_gather = clus_df_gather_sorted, graphing_columns = graphing_columns_sorted))
    } else {
        return(clus_df_gather_sorted)
    }
}

#' Generate an Alluvial Plot with Minimal Cluster Cross-over
#'
#' Creates a two-axis alluvial plot to visualize the relationship between two categorical groupings (e.g., cluster assignments from different methods),
#' minimizing crossover and aligning clusters based on agreement.
#'
#' @param df A data frame, tibble, or CSV file path. Must be in one of two formats:
#' (1) column_weights == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two graphing_columns).
#' (2) column_weights != NULL: Each row represents a combination of groupings, each column from \code{graphing_columns} represents a grouping, and the column \code{column_weights} represents the number of entities in that combination of groupings. Must contain at least three columns (two \code{graphing_columns}, one \code{column_weights}).
#' @param graphing_columns Optional character vector. Vector of column names from \code{df} to be used in graphing (i.e., alluvial plotting). Mutually exclusive with \code{column1} and \code{column2}.
#' @param column1 Optional character. Can be used along with \code{column2} in place of \code{graphing_columns} if working with two columns only. Mutually exclusive with \code{graphing_columns}.
#' @param column2 Optional character. Can be used along with \code{column1} in place of \code{graphing_columns} if working with two columns only. Mutually exclusive with \code{graphing_columns}.
#' @param column_weights Optional character. Column name from \code{df} that contains the weights of each combination of groupings if \code{df} is in format (2) (see above).
#' @param sorting_algorithm Character. Algorithm with which to sort the values in the dataframe. Can choose from {'neighbornet', 'greedy_WOLF', 'greedy_WBLF', 'None'}. 'neighbornet' performs sorting with NeighborNet (Bryant and Moulton, 2004). 'greedy_WOLF' implements a custom greedy algorithm where one layer is fixed, and the other layer is sorted such that each node is positioned as close to its largest parent from the fixed side as possible in a greedy fashion. 'greedy_WBLF' implements the 'greedy_WOLF' algorithm described previously twice, treating each column as fixed in one iteration and free in the other iteration. 'greedy_WOLF' and 'greedy_WBLF' are only valid when \code{graphing_columns} has exactly two entries.
#' @param optimize_column_order Logical. If TRUE, will optimize the order of \code{graphing_columns} to minimize edge overlap. Only applies when \code{sorting_algorithm == 'neighbornet'} and \code{length(graphing_columns) > 2}.
#' @param optimize_column_order_per_cycle Logical. If TRUE, will optimize the order of \code{graphing_columns} to minimize edge overlap upon each cycle. If FALSE, will optimize the order of \code{graphing_columns} to minimize edge overlap on the beginning cycle only. Only applies when \code{sorting_algorithm == 'neighbornet'} and \code{length(graphing_columns) > 2}.
#' @param fixed_column Character or Integer. Name or position of the column in \code{graphing_columns} to keep fixed during sorting. Only applies when \code{sorting_algorithm == 'greedy_WOLF'}.
#' @param random_initializations Integer. Number of random initializations for the positions of each grouping in \code{graphing_columns}. Only applies when \code{sorting_algorithm == 'greedy_WOLF' or sorting_algorithm == 'greedy_WBLF'}.
#' @param set_seed Integer. Random seed for the \code{random_initializations} parameter. Only applies when \code{sorting_algorithm == 'greedy_WOLF' or sorting_algorithm == 'greedy_WBLF'}.
#' @param color_boxes Logical. Whether to color the strata/boxes (representing groups).
#' @param color_bands Logical. Whether to color the alluvia/edges (connecting the strata).
#' @param color_list Optional named list or vector of colors to override default group colors.
#' @param color_band_list Optional named list or vector of colors to override default band colors.
#' @param color_band_column Optional Character. Which column to use for coloring bands.
#' @param color_band_boundary Logical. Whether or not to color boundaries between bands
#' @param match_colors Logical. If \code{TRUE}, assigns consistent colors between column1 and column2 where matched.
#' @param alluvial_alpha Numeric between 0 and 1. Transparency level for the alluvial bands.
#' @param include_labels_in_boxes Logical. Whether to include text labels inside the rectangular group boxes.
#' @param include_axis_titles Logical. Whether to display axis titles for column1 and column2.
#' @param include_group_sizes Logical. If \code{TRUE}, includes group sizes in the labels (e.g., "Group A (42)").
#' @param output_plot_path Character. File path to save the plot (e.g., "plot.png"). If \code{NULL}, then will not be saved.
#' @param output_df_path Optional character. Output path for the output data frame, in CSV format. If \code{NULL}, then will not be saved.
#' @param preprocess_data Logical. If TRUE, will preprocess the data with the \code{data_preprocess} function.
#' @param verbose Logical. If TRUE, will display messages during the function.
#'
#' @return A \code{ggplot2} object representing the alluvial plot.
#'
#' @examples
#' # Example 1: df format 1
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' p <- plot_alluvial(df, graphing_columns = c('method1', 'method2'))
#'
#' # Example 2: df format 2
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- df |>
#'     dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
#'     dplyr::group_by_all() |>
#'     dplyr::count(name = "value")
#' p <- plot_alluvial(clus_df_gather, graphing_columns = c('method1', 'method2'), column_weights = 'value')
#'
#' @export
plot_alluvial <- function(df, graphing_columns = NULL, column1 = NULL, column2 = NULL, 
                          column_weights = NULL, sorting_algorithm = 'neighbornet', optimize_column_order=TRUE, 
                          optimize_column_order_per_cycle=FALSE, fixed_column = NULL, random_initializations = 1, 
                          set_seed=42, color_boxes = TRUE, color_bands = FALSE, color_list = NULL, 
                          color_band_list = NULL, color_band_column=NULL, color_band_boundary=FALSE, 
                          match_colors = TRUE, alluvial_alpha = 0.5, include_labels_in_boxes = TRUE, 
                          include_axis_titles = TRUE, include_group_sizes = TRUE, output_plot_path = NULL, 
                          output_df_path = NULL, preprocess_data=TRUE, verbose=FALSE,
                          box_width = 1/3, text_width = 1/4, min_text = 4, 
                          save_height = 6, save_width = 6) {
    #* Type Checking Start
    # ensure someone doesn't specify both graphing_columns and column1/2
    if (!is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        stop("Specify either graphing_columns or column1/column2, not both.")
    }

    if (verbose) message("Loading in data")
    df <- load_in_df(df = df, graphing_columns = graphing_columns, column_weights = column_weights)

    if (!is.null(graphing_columns) && any(!graphing_columns %in% colnames(df))) {
        stop("Some graphing_columns are not present in the dataframe.")
    }

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

    # Preprocess
    if (preprocess_data) {
        if (verbose) message("Preprocessing data before sorting")
        clus_df_gather_unsorted <- data_preprocess(df = df, graphing_columns = graphing_columns, column_weights = column_weights, load_df = FALSE, do_gather_set_data = FALSE)
        if (is.null(column_weights)) {
            column_weights <- "value"  # is set during data_preprocess
        }
    } else {
        clus_df_gather_unsorted <- df
    }

    # Sort
    if (verbose) message(sprintf("Sorting data with sorting_algorithm=%s", sorting_algorithm))
    data_sort_output <- data_sort(df = clus_df_gather_unsorted, graphing_columns = graphing_columns, column_weights = column_weights, sorting_algorithm = sorting_algorithm, optimize_column_order = optimize_column_order, optimize_column_order_per_cycle = optimize_column_order_per_cycle, fixed_column = fixed_column, random_initializations = random_initializations, set_seed = set_seed, output_df_path = output_df_path, return_updated_graphing_columns = TRUE, preprocess_data = FALSE, load_df = FALSE, verbose = verbose)
    clus_df_gather <- data_sort_output$clus_df_gather
    graphing_columns <- data_sort_output$graphing_columns

    # Plot
    if (verbose) message("Plotting data")
    alluvial_plot <- plot_alluvial_internal(clus_df_gather, graphing_columns=graphing_columns,
                                                         color_list = color_list, color_boxes = color_boxes,
                                                         color_bands = color_bands, color_band_list = color_band_list,
                                                         color_band_column=color_band_column, color_band_boundary=color_band_boundary,
                                                         match_colors = match_colors, alluvial_alpha = alluvial_alpha,
                                                         include_labels_in_boxes = include_labels_in_boxes, include_axis_titles = include_axis_titles,
                                                         include_group_sizes = include_group_sizes,
                                                         output_plot_path = output_plot_path, verbose = verbose,
                                            box_width = box_width, text_width = text_width, min_text = min_text, 
                                            save_height = save_height, save_width = save_width)

    return(alluvial_plot)
}
