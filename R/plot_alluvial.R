#' wompwomp: Cluster-matching alluvial plots
#'
#' Main plotting function and helpers for bipartite-matching-based alluvial diagrams
#' @docType package
#' @name wompwomp
#'
#' @importFrom dplyr mutate select group_by summarise arrange desc ungroup slice n pull filter bind_rows across matches all_of add_count distinct
#' @importFrom ggplot2 ggplot aes geom_text scale_fill_manual labs after_stat annotate theme_void theme element_text rel ggsave guides scale_color_manual scale_x_continuous element_blank
#' @importFrom ggalluvial geom_alluvium geom_stratum stat_stratum stat_alluvium
#' @importFrom ggfittext geom_fit_text
#' @importFrom ggforce gather_set_data
#' @importFrom purrr map
#' @importFrom igraph max_bipartite_match V graph_from_data_frame cluster_louvain cluster_leiden E
#' @importFrom tibble is_tibble
#' @importFrom utils read.csv write.csv combn
#' @importFrom stats setNames
#' @importFrom rlang sym .data
#' @importFrom magrittr %>%
#' @importFrom data.table :=

utils::globalVariables(c(
    ".data", ":=", "%>%", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "stratum", "total", "cum_y", "best_cluster_agreement", "neighbor_net", "alluvium", "pos", "count", "group1", "group2", "value", "group1_size", "group2_size", "weight", "parent", "group_name"
))

StatStratum <- ggalluvial::StatStratum # avoid the error Can't find stat called "stratum" - and make sure to do stat = StatStratum instead of stat = "stratum"

neighbornet_script_path <- system.file("scripts", "run_neighbornet.py")
if (neighbornet_script_path == "") {
    # Fallback to development location
    neighbornet_script_path <- file.path(here::here("inst", "scripts", "run_neighbornet.py"))
}
stopifnot(file.exists(neighbornet_script_path))

# reticulate::source_python(neighbornet_script_path)  # Error: Unable to access object (object is from previous session and is now invalid)

# library(dplyr)
# library(ggplot2)
# library(ggalluvial)
# library(ggforce)
# library(igraph)
# library(tibble)

default_colors <- c(
    "#D55E00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#E69F00", "#CC79A7", "#666666", "#AD7700", "#1C91D4", "#007756", "#D5C711", "#005685",
    "#A04700", "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71", "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C", "#FFCB57", "#9AD2F2",
    "#2CFFC6", "#F6EF8E", "#38B7FF", "#FF9B4D", "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45", "#AA9F0D", "#00446B", "#803800", "#8D3666",
    "#3D3D3D"
)

compute_alluvial_statistics <- function(clus_df_gather, graphing_columns, column_weights = "value") {
    message(sprintf("Alluvial statistics: n = number of elements; m = number of graphing columns; a = number of alluvia/edges; k_i = number of blocks in layer i (where i goes from 1:m); K_sum = number of blocks across all layers; K_prod = product of blocks across all layers"))
    message(sprintf("n = %s", sum(clus_df_gather[[column_weights]], na.rm = TRUE)))
    message(sprintf("m = %s", length(graphing_columns)))
    message(sprintf("a = %s", nrow(clus_df_gather)))

    K_sum <- 0
    K_prod <- 1
    for (i in 1:length(graphing_columns)) {
        colname <- paste0("col", i, "_int")
        k_i <- length(unique(clus_df_gather[[colname]]))
        message(sprintf("k_%s = %s", i, k_i))

        K_sum <- K_sum + k_i
        K_prod <- K_prod * k_i
    }

    message(sprintf("K_sum = %s", K_sum))
    message(sprintf("K_prod = %s", K_prod))
}

determine_column_order <- function(clus_df_gather_neighbornet, graphing_columns, column_weights = "value", matrix_initialization_value_column_order = 1e6, weight_scalar_column_order = 1, column_sorting_metric = "edge_crossing", column_sorting_algorithm = "tsp", verbose = FALSE, environment = "wompwomp_env", use_conda = TRUE) {
    # this doesn't strictly need its own condition (2 choose 2 is 1 anyways), but does avoid a little overhead
    if (length(graphing_columns) == 2) {
        return(graphing_columns)
    }

    if (column_sorting_metric == "ARI") {
        if (!requireNamespace("mclust", quietly = TRUE)) {
            stop("The 'mclust' package is required to compute Adjusted Rand Index (ARI) with column_sorting_metric == 'ARI'. Please install it with install.packages('mclust').")
        }
    }

    column_dist_matrix <- matrix(matrix_initialization_value_column_order,
        nrow = length(graphing_columns), ncol = length(graphing_columns),
        dimnames = list(graphing_columns, graphing_columns)
    )

    pairs <- combn(graphing_columns, 2)
    if (verbose) message("Computing objectives for each pair of columns in order to determine column order")
    for (i in 1:ncol(pairs)) {
        if (verbose) message(sprintf("Computing objective for column pairs %s / %s", i, ncol(pairs)))
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

        if (column_sorting_metric == "ARI") {
            expanded_df <- clus_df_gather_neighbornet_tmp[rep(seq_len(nrow(clus_df_gather_neighbornet_tmp)), clus_df_gather_neighbornet_tmp[[column_weights]]), ]
            neighbornet_objective <- mclust::adjustedRandIndex(expanded_df$col1_int, expanded_df$col2_int)
            neighbornet_objective <- -neighbornet_objective + 1 # convert from [-1,1] to [0,2], and flip the sign (so that 1 becomes smallest ie perfect cluster agreement --> smallest distance)
            neighbornet_objective <- weight_scalar_column_order * 50 * neighbornet_objective # built-in scalar
        } else if (column_sorting_metric == "edge_crossing") {
            neighbornet_objective <- determine_crossing_edges(
                clus_df_gather_neighbornet_tmp,
                graphing_columns = graphing_columns_tmp,
                column_weights = column_weights,
                return_weighted_layer_free_objective = TRUE,
                load_df = FALSE,
                # verbose = verbose,
                preprocess_data = FALSE, environment = environment, use_conda = use_conda
            )
            neighbornet_objective <- weight_scalar_column_order * log1p(neighbornet_objective) # log1p to avoid issue of log(0)
        } else {
            stop(sprintf("column_sorting_metric '%s' is not a valid option.", column_sorting_metric))
        }

        column_dist_matrix[column1, column2] <- neighbornet_objective
        column_dist_matrix[column2, column1] <- neighbornet_objective
    }
    # Prepare data in R
    labels <- graphing_columns # assuming this is a character vector
    # # Call Python function
    # mat_list <- split(column_dist_matrix, row(column_dist_matrix))  # convert R matrix to list of row-vectors
    # result <- nn_mod$neighbor_net(labels, mat_list)
    if (verbose) message(sprintf("Running '%s' for column order", column_sorting_algorithm))

    if (column_sorting_algorithm == "tsp") {
        tsp_instance <- TSP::TSP(column_dist_matrix)
        tour <- TSP::solve_TSP(tsp_instance)
        cycle <- as.integer(tour)
    } else if (column_sorting_algorithm == "neighbornet") {
        reticulate::source_python(neighbornet_script_path)
        result <- neighbor_net(labels, column_dist_matrix) # from python
        cycle <- result[[1]]
    } else {
        stop(sprintf("column_sorting_algorithm '%s' is not a valid option.", column_sorting_algorithm))
    }

    cycle_mapped <- labels[cycle]

    # determine the optimal starting point for cycle
    adj_distances <- sapply(seq_len(length(cycle_mapped)), function(i) {
        from <- cycle_mapped[i]
        to <- cycle_mapped[(i %% length(cycle_mapped)) + 1] # wraps around
        column_dist_matrix[from, to]
    })
    max_index <- which.max(adj_distances)

    cycle_mapped_optimal_start <- rotate_left(cycle_mapped, max_index)
    if (verbose) message("Done with neighbornet for column order")
    return(cycle_mapped_optimal_start)
}

run_neighbornet <- function(df, graphing_columns = NULL, column1 = NULL, column2 = NULL, column_weights = "value", matrix_initialization_value = 1e6, same_side_matrix_initialization_value = 1e6, weight_scalar = 5e5, sorting_algorithm = "neighbornet", verbose = FALSE) {
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
        clus_df_gather <- get_alluvial_df(df, column_weights = column_weights)
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
    full_dist_matrix <- matrix(matrix_initialization_value,
        nrow = length(all_nodes), ncol = length(all_nodes),
        dimnames = list(all_nodes, all_nodes)
    )

    # Same-side initialization
    if (same_side_matrix_initialization_value != matrix_initialization_value) {
        for (col in graphing_columns) {
            prefix <- paste0(col, "~~")
            node_indices <- which(startsWith(all_nodes, prefix))
            full_dist_matrix[node_indices, node_indices] <- same_side_matrix_initialization_value
        }
    }

    # Get all 2-column combinations
    pairwise_groupings <- combn(graphing_columns, 2, simplify = FALSE)

    # For each combination, group and summarize
    summarized_results <- purrr::map(pairwise_groupings, function(cols) {
        clus_df_gather %>%
            group_by(across(all_of(cols))) %>%
            summarise(total_value = sum(!!sym(column_weights)), .groups = "drop") %>%
            mutate(grouping = paste(cols, collapse = "+"))
    })

    # summarized_results <- purrr::map(pairwise_groupings, function(cols) {
    #     clus_df_gather %>%
    #         group_by(across(all_of(cols))) %>%
    #         summarise(total_value = sum(value), .groups = "drop") %>%
    #         mutate(grouping = paste(cols, collapse = "+"))
    # })

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
            full_dist_matrix[n1, n2] <- weight_scalar * -log(w)
            full_dist_matrix[n2, n1] <- weight_scalar * -log(w) # symmetric since graph is undirected
        }
    }

    # make sure all numbers are positive for neighbornet
    min_val_abs <- abs(min(full_dist_matrix))
    full_dist_matrix <- full_dist_matrix + (min_val_abs + 1)

    # Prepare data in R
    labels <- all_nodes # assuming this is a character vector
    mat <- full_dist_matrix
    mat[is.infinite(mat)] <- 1e6
    mat[is.na(mat)] <- 1e6

    # # Call Python function
    # mat_list <- split(mat, row(mat))  # convert R matrix to list of row-vectors
    # result <- nn_mod$neighbor_net(labels, mat_list)
    if (verbose) message(sprintf("Running '%s' for stratum order", sorting_algorithm))

    if (sorting_algorithm == "tsp") {
        tsp_instance <- TSP::TSP(mat)
        tour <- TSP::solve_TSP(tsp_instance)
        cycle <- as.integer(tour)
    } else if (sorting_algorithm == "neighbornet") {
        reticulate::source_python(neighbornet_script_path)
        result <- neighbor_net(labels, mat) # from python
        cycle <- result[[1]]
        # splits <- result[[2]]
    }

    cycle_mapped <- labels[cycle]

    return(cycle_mapped)
}

rotate_left <- function(vec, k = 1) {
    n <- length(vec)
    k <- k %% n
    if (k == 0) {
        return(vec)
    }
    c(vec[(k + 1):n], vec[1:k])
}

get_graph_groups <- function(cycle) {
    groups <- list()

    for (node in cycle) {
        prefix <- sub("~~.*", "", node) # Extract everything before the first `~~`

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



determine_optimal_cycle_start <- function(df, cycle, graphing_columns = NULL, column1 = NULL, column2 = NULL, column_weights = "value", optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE, matrix_initialization_value_column_order = 1e6, weight_scalar_column_order = 1, column_sorting_metric = "edge_crossing", column_sorting_algorithm = "tsp", cycle_start_positions = NULL, verbose = FALSE, make_intermediate_neighbornet_plots = FALSE, environment = "wompwomp_env", use_conda = TRUE) {
    # ensure someone doesn't specify both graphing_columns and column1/2
    if (!is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        stop("Specify either graphing_columns or column1/column2, not both.")
    }
    # if someone specifies column1/2, then use it
    if (is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        graphing_columns <- c(column1, column2)
    }

    # # Commented out because I'm not sold on ARI
    # if (optimize_column_order_per_cycle && (column_sorting_metric == "edge_crossing")) {
    #     if (verbose) message("column_sorting_metric == 'edge_crossing' and optimize_column_order_per_cycle is TRUE. This might be a bit slow. Consider setting column_sorting_metric == 'ARI' and/or optimize_column_order_per_cycle to FALSE.")
    # }

    # factorize input columns
    for (col in graphing_columns) {
        df[[col]] <- as.factor(as.character(df[[col]]))
    }

    neighbornet_objective_minimum <- Inf
    # p_best_neighbornet <- NULL
    cycle_best <- NULL
    clus_df_gather_best <- NULL
    graphing_columns_best <- NULL
    objective_matrix_vector <- c()

    n <- length(cycle)

    graphing_columns_tmp <- graphing_columns
    for (i in 0:(n - 1)) {
        if ((!is.null(cycle_start_positions)) && !((i + 1) %in% cycle_start_positions)) {
            next
        }
        if (verbose) message(sprintf("Starting iteration %s / %s", i + 1, n))
        # if (i == 0) {
        #     if (verbose) message(sprintf("Starting iteration 1"))
        # } else if (i == 1) {
        #     if (verbose) message(sprintf("Starting subsequent iterations (should go much faster than iteration 1 if optimize_column_order is FALSE and/or optimize_column_order_per_cycle is FALSE)"))
        # }

        cycle_shifted <- rotate_left(cycle, i)
        graphs_list <- get_graph_groups(cycle_shifted)

        # remove prefix (column1_, etc)
        graphs_list_stripped <- lapply(graphs_list, function(x) {
            sub("^.*?~~", "", x)
        })

        if (is.null(column_weights) || !(column_weights %in% colnames(df))) {
            clus_df_gather_neighbornet <- get_alluvial_df(df, column_weights = column_weights)
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

        graphing_columns_tmp_previous_iteration <- graphing_columns_tmp
        if (optimize_column_order) {
            # optimize order either on the first iteration if optimize_column_order_per_cycle is FALSE, or each time if optimize_column_order_per_cycle is TRUE
            if ((optimize_column_order_per_cycle) || (i == 0)) {
                # verbose_tmp <- verbose
                verbose_tmp <- if (i == 0) verbose else FALSE # only have the option for verbose on first iteration
                graphing_columns_tmp <- determine_column_order(clus_df_gather_neighbornet, graphing_columns = graphing_columns, column_weights = column_weights, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_sorting_metric = column_sorting_metric, column_sorting_algorithm = column_sorting_algorithm, verbose = verbose_tmp, environment = environment, use_conda = use_conda)
            }
        }

        clus_df_gather_neighbornet_tmp <- reorder_and_rename_columns(clus_df_gather_neighbornet, graphing_columns_tmp)

        neighbornet_objective <- determine_crossing_edges(
            clus_df_gather_neighbornet_tmp,
            graphing_columns = graphing_columns_tmp,
            column_weights = column_weights,
            # verbose = verbose,
            return_weighted_layer_free_objective = TRUE, environment = environment, use_conda = use_conda
        )

        # # if ((optimize_column_order_per_cycle) || (i == 0)) {
        # # # recalculate from scratch if i == 0 (first iteration, no matrices made yet) OR if graphing_columns_tmp != graphing_columns_tmp_previous_iteration (optimize_column_order_per_cycle is TRUE and we have a new column order, and therefore our matrices are no help)
        # if ((!all(graphing_columns_tmp == graphing_columns_tmp_previous_iteration)) || (i == 0)) {
        #     # include_output_objective_matrix_vector <- !optimize_column_order_per_cycle  # if optimize_column_order_per_cycle is TRUE, then no need to return this big matrix
        #     neighbornet_objective_output <- determine_crossing_edges(
        #         clus_df_gather_neighbornet_tmp,
        #         graphing_columns = graphing_columns_tmp,
        #         column_weights = column_weights,
        #         verbose = verbose,
        #         include_output_objective_matrix_vector = TRUE,
        #         return_weighted_layer_free_objective = FALSE,
        #         environment = environment, use_conda = use_conda
        #     )
        # } else {
        #     swapped_node <- cycle_shifted[length(cycle_shifted)]
        #     parts <- strsplit(swapped_node, "~~", fixed = TRUE)[[1]]
        #     layer_name <- parts[1]
        #     node_name <- parts[2]
        #
        #     # determine the int column that matches to this layer_name - because I haven't reordered any columns in clus_df_gather_neighbornet, I should use the order as determined here
        #     int_column_int <- match(layer_name, graphing_columns_tmp)
        #     int_column <- paste0("col", int_column_int, "_int")
        #
        #     # find the value in clus_df_gather_neighbornet[[int_column]] that maps to node_name of clus_df_gather_neighbornet[[layer_name]]
        #     matched <- clus_df_gather_neighbornet_tmp[[int_column]][clus_df_gather_neighbornet_tmp[[layer_name]] == node_name]
        #     stratum_int_name <- matched[1]
        #     stratum_column_and_value_to_keep <- setNames(list(stratum_int_name), as.character(int_column_int))
        #
        #     neighbornet_objective_output <- determine_crossing_edges(
        #         clus_df_gather_neighbornet_tmp,
        #         graphing_columns = graphing_columns_tmp,
        #         column_weights = column_weights,
        #         stratum_column_and_value_to_keep = stratum_column_and_value_to_keep,
        #         input_objective = neighbornet_objective,
        #         input_objective_matrix_vector = objective_matrix_vector,
        #         verbose = verbose,
        #         include_output_objective_matrix_vector = TRUE,
        #         return_weighted_layer_free_objective = FALSE,
        #         environment = environment, use_conda = use_conda
        #     )
        # }
        #
        # objective_matrix_vector <- neighbornet_objective_output$objective_matrix_vector
        # neighbornet_objective <- neighbornet_objective_output$output_objective

        # to save each plot
        if (make_intermediate_neighbornet_plots) {
            output_plot_path <- sprintf("neighbornet_intermediate_%s.png", i + 1)
            for (j in seq_along(graphing_columns_tmp)) {
                int_col_name <- paste0("col", j, "_int")
                clus_df_gather_neighbornet_tmp[[int_col_name]] <- factor(clus_df_gather_neighbornet_tmp[[int_col_name]])
            }
            plot_alluvial(clus_df_gather_neighbornet_tmp, graphing_columns = graphing_columns, column_weights = column_weights, sorting_algorithm = "none", preprocess_data = FALSE, color_boxes = FALSE, color_bands = TRUE, output_plot_path = output_plot_path)
        }

        # print(neighbornet_objective)
        if (verbose) message(sprintf("neighbornet_objective for iteration %s = %s", i + 1, neighbornet_objective))
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
        clus_df_gather$y <- -2 # tmp
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
            clus_df_gather[[reordered_column]]
        )
    }

    return(clus_df_gather)
}

find_colors_advanced <- function(clus_df_gather, ditto_colors, graphing_columns, coloring_algorithm_advanced_option = "leiden", resolution = 1) {
    column_int_names <- c()
    for (col_int in seq_along(graphing_columns)) {
        int_name <- paste0("col", col_int, "_int")
        column_int_names <- c(column_int_names, int_name)
    }
    clus_df_ungrouped <- clus_df_gather[, c(column_int_names, "value")]

    first <- TRUE
    compared <- c()
    for (group1_name in column_int_names) {
        for (group2_name in column_int_names) {
            if (!(group1_name == group2_name)) {
                comp1 <- paste0(group1_name, group2_name)
                comp2 <- paste0(group2_name, group1_name)
                if (!(comp1 %in% compared | comp2 %in% compared)) {
                    if (first) {
                        clus_df_filtered <- clus_df_ungrouped[, c(group1_name, group2_name, "value")]
                        clus_df_filtered <- clus_df_filtered %>%
                            add_count(!!rlang::sym(group1_name), !!rlang::sym(group2_name), wt = value) %>%
                            select(!!rlang::sym(group1_name), !!rlang::sym(group2_name), n)
                        clus_df_filtered <- distinct(clus_df_filtered)
                        colnames(clus_df_filtered) <- c("group1", "group2", "value")

                        clus_df_filtered <- clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(group1_size = sum(value))
                        clus_df_filtered <- clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(group2_size = sum(value))

                        clus_df_filtered <- clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(weight = value) # /group2_size)

                        clus_df_filtered$group1 <- sub("^", paste0(group1_name, "_"), clus_df_filtered[["group1"]])
                        clus_df_filtered$group2 <- sub("^", paste0(group2_name, "_"), clus_df_filtered[["group2"]])

                        first <- FALSE
                    } else {
                        temp_clus_df_filtered <- clus_df_ungrouped[, c(group1_name, group2_name, "value")]
                        temp_clus_df_filtered <- temp_clus_df_filtered %>%
                            add_count(!!rlang::sym(group1_name), !!rlang::sym(group2_name), wt = value) %>%
                            select(!!rlang::sym(group1_name), !!rlang::sym(group2_name), n)
                        temp_clus_df_filtered <- distinct(temp_clus_df_filtered)
                        colnames(temp_clus_df_filtered) <- c("group1", "group2", "value")

                        temp_clus_df_filtered <- temp_clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(group1_size = sum(value))
                        temp_clus_df_filtered <- temp_clus_df_filtered %>%
                            group_by(group2) %>%
                            mutate(group2_size = sum(value))

                        temp_clus_df_filtered <- temp_clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(weight = value) # /group2_size)

                        temp_clus_df_filtered$group1 <- sub("^", paste0(group1_name, "_"), temp_clus_df_filtered[["group1"]])
                        temp_clus_df_filtered$group2 <- sub("^", paste0(group2_name, "_"), temp_clus_df_filtered[["group2"]])

                        clus_df_filtered <- rbind(clus_df_filtered, temp_clus_df_filtered)
                    }
                }
                # compared <- c(compared, comp1, comp2)
            }
        }
    }

    clus_df_extra_filtered <- clus_df_filtered[, c("group1", "group2", "value")]
    g <- igraph::graph_from_data_frame(d = clus_df_extra_filtered, directed = FALSE)
    if (coloring_algorithm_advanced_option == "louvain") {
        partition <- igraph::cluster_louvain(g, weights = igraph::E(g)$value, resolution = resolution)
    } else if (coloring_algorithm_advanced_option == "leiden") {
        partition <- igraph::cluster_leiden(g, weights = igraph::E(g)$value, resolution_parameter = resolution)
    } else {
        stop(sprintf("coloring_algorithm_advanced_option '%s' is not recognized. Please choose from 'leiden' (default) or 'louvain'.", coloring_algorithm_advanced_option))
    }

    clus_df_leiden <- data.frame(group_name = partition$names, leiden = partition$membership)
    clus_df_leiden <- clus_df_leiden %>% tidyr::separate_wider_delim(group_name, names = c("col", "trash", "group"), delim = "_")
    clus_df_leiden <- clus_df_leiden %>% tidyr::separate_wider_delim(col, names = c("trash2", "col"), delim = "col")
    clus_df_leiden <- clus_df_leiden[, c("col", "group", "leiden")]
    clus_df_leiden[["colors"]] <- unlist(Map(function(x) ditto_colors[x], clus_df_leiden$leiden))

    final_colors <- c()
    for (col_int in unique(clus_df_leiden$col)) {
        temp_df <- clus_df_leiden[clus_df_leiden$col == col_int, ]
        final_colors <- c(final_colors, temp_df[rev(order(temp_df$group)), ]$colors)
    }

    return(final_colors)
}

find_group2_colors <- function(clus_df_gather, ditto_colors, unused_colors, current_g1_colors,
                               group1_name = "col1_int", group2_name = "col2_int",
                               cutoff = .5) {
    num_levels <- length(levels(clus_df_gather[[group2_name]]))

    clus_df_ungrouped <- clus_df_gather[, c(group1_name, group2_name, "value")]
    clus_df_filtered <- clus_df_ungrouped %>%
        add_count(!!rlang::sym(group1_name), !!rlang::sym(group2_name), wt = value) %>%
        select(!!rlang::sym(group1_name), !!rlang::sym(group2_name), n)
    clus_df_filtered <- distinct(clus_df_filtered)
    colnames(clus_df_filtered) <- c(group1_name, group2_name, "value")


    clus_df_filtered <- clus_df_filtered %>%
        group_by(!!rlang::sym(group1_name)) %>%
        mutate(group1_size = sum(value))
    clus_df_filtered <- clus_df_filtered %>%
        group_by(!!rlang::sym(group2_name)) %>%
        mutate(group2_size = sum(value))
    clus_df_filtered <- clus_df_filtered %>%
        group_by(!!rlang::sym(group1_name)) %>%
        mutate(weight = value / group2_size)

    parent_df <- clus_df_filtered %>%
        group_by(!!rlang::sym(group2_name)) %>%
        filter(weight == max(weight)) # %>% select(!!rlang::sym(group1_name))
    colnames(parent_df)[which(names(parent_df) == group1_name)] <- "parent"
    parent_df <- parent_df %>% mutate(parent = (ifelse(weight > cutoff, as.character(parent), "")))
    parent_df <- parent_df[, c(group2_name, "parent")]
    clus_df_filtered <- merge(clus_df_filtered, parent_df, by = group2_name)

    clus_df_filtered[["colors"]] <- unlist(Map(function(x) ifelse(x == "", "", current_g1_colors[x]), clus_df_filtered$parent))
    clus_df_parent <- clus_df_filtered[clus_df_filtered$parent != "", ]
    group2_colors <- vector("character", num_levels)
    for (col_int in unique(clus_df_gather[[group2_name]])) {
        if (col_int %in% unique(clus_df_parent[[group2_name]])) {
            temp_df <- clus_df_parent[clus_df_parent[[group2_name]] == col_int, ]
            temp_df <- temp_df[temp_df[[group1_name]] == temp_df[["parent"]], ]
            group2_colors[[as.integer(col_int)]] <- temp_df[["colors"]]
        } else {
            group2_colors[[as.integer(col_int)]] <- unused_colors[1]
            unused_colors <- unused_colors[2:length(unused_colors)]
        }
    }

    return(group2_colors)
}

#' Just the plotting
#'
#' Just the plotting.......
#'
#' @param df A data frame, tibble, or CSV file path. Must be in the format as the output of wompwomp::data_sort.
#' @param graphing_columns Character vector. Vector of column names from \code{df} to be used in graphing (i.e., alluvial plotting). Mutually exclusive with \code{column1} and \code{column2}.
#' @param column_weights Character. Column name from \code{df} that contains the weights of each combination of groupings if \code{df} is in format (2) (see above).
#' @param color_boxes Logical. Whether to color the strata/boxes (representing groups).
#' @param color_bands Logical. Whether to color the alluvia/edges (connecting the strata).
#' @param color_list Optional named list or vector of colors to override default group colors.
#' @param color_band_list Optional named list or vector of colors to override default band colors.
#' @param color_band_column Optional Character. Which column to use for coloring bands.
#' @param color_val Optional named list where the entries are colors and the names correspond to values of the dataframe that should use those colors
#' @param color_band_boundary Logical. Whether or not to color boundaries between bands
#' @param coloring_algorithm Character. Matching colors methods. Choices are 'advanced' (default), 'none', 'left', 'right', or any value in \code{graphing_columns}.
#' @param coloring_algorithm_advanced_option Character. If \code{coloring_algorithm == 'advanced'}, then choose graph clustering algorithm. Choices are 'leiden' (default) or 'louvain'.
#' @param resolution Numeric If \code{coloring_algorithm == 'advanced'}, then choose resolution for the graph clustering algorithm. Affects coloring of both bands and boxes.
#' @param cutoff Numeric If \code{coloring_algorithm != 'none' and coloring_algorithm != 'advanced'}, sets the cutoff for color matching, below which a new color will be assigned.
#' @param alluvial_alpha Numeric between 0 and 1. Transparency level for the alluvial bands.
#' @param include_labels_in_boxes Logical. Whether to include text labels inside the rectangular group boxes.
#' @param include_axis_titles Logical. Whether to display axis titles for column1 and column2.
#' @param include_group_sizes Logical. If \code{TRUE}, includes group sizes in the labels (e.g., "Group A (42)").
#' @param output_plot_path Character. File path to save the plot (e.g., "plot.png"). If \code{NULL}, then will not be saved.
#' @param box_width Numeric between 0 and 1. Box width
#' @param text_width Numeric between 0 and 1. Text width
#' @param min_text Integer greater than 0. Min text
#' @param text_size Integer greater than 0. Text size (works whether auto_adjust_text is TRUE or FALSE).
#' @param auto_adjust_text Logical. Whether to automatically adjust text size to fit in box.
#' @param axis_text_size Integer greater than 0. Axis text size
#' @param axis_text_vjust Integer. Axis text vjust
#' @param save_height Integer greater than 0. Save height, in inches
#' @param save_width Integer greater than 0. Save width, in inches
#' @param dpi Integer greater than 0. DPI for \code{output_plot_path}, if \code{output_plot_path} is a raster image or \code{rasterise_alluvia} is TRUE
#' @param rasterise_alluvia Logical. Whether to rasterize the alluvia if \code{output_plot_path} is a PDF. Can save space if DPI low enough
#' @param keep_y_labels Keep y labels
#' @param box_line_width Box line width
#' @param environment Character. Python environment (if applicable). Default: 'wompwomp_env'
#' @param use_conda Logical. Whether or not to use conda for Python (if applicable)
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param print_params Logical. If TRUE, will print function params.
#' @param add_legend Logical. If TRUE, will generate a legend of the colors of boxes and alluvial
#' @param legend_loc Character. Location of legend. Only applies if \code{add_legened == TRUE}. Choices are 'right' (default), 'left', 'bottom', 'top'
#' @param do_compute_alluvial_statistics Internal flag; not recommended to modify.
#'
#' @return A \code{ggplot2} object representing the alluvial plot.
#'
#' @examples
#' # Example:
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- df |>
#'     dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
#'     dplyr::group_by_all() |>
#'     dplyr::count(name = "value")
#' clus_df_gather <- data_sort(
#'     clus_df_gather,
#'     graphing_columns = c("method1", "method2"),
#'     column_weights = "value",
#'     sorting_algorithm = "none"
#' )
#' p <- plot_alluvial_internal(
#'     clus_df_gather,
#'     graphing_columns = c("method1", "method2"),
#'     column_weights = "value",
#'     coloring_algorithm = "right"
#' )
#'
#' @export
plot_alluvial_internal <- function(df, graphing_columns, column_weights,
                                   color_list = NULL, color_boxes = TRUE,
                                   color_bands = FALSE, color_band_list = NULL,
                                   color_band_column = NULL, color_band_boundary = FALSE,
                                   alluvial_alpha = 0.5, coloring_algorithm = "advanced", cutoff = .5,
                                   output_plot_path = NULL, color_val = NULL,
                                   include_labels_in_boxes = TRUE, include_axis_titles = FALSE,
                                   include_group_sizes = FALSE, verbose = FALSE, print_params = FALSE,
                                   box_width = 1 / 3, text_width = 1 / 4, min_text = 4, text_size = 14,
                                   auto_adjust_text = TRUE, axis_text_size = 2, axis_text_vjust = 0,
                                   save_height = 6, save_width = 6, dpi = 300, rasterise_alluvia = FALSE,
                                   keep_y_labels = FALSE, box_line_width = 1, do_compute_alluvial_statistics = TRUE,
                                   coloring_algorithm_advanced_option = "leiden", resolution = 1, environment = "wompwomp_env", use_conda = TRUE,
                                   add_legend = FALSE, legend_loc = "right") {
    if (print_params) print_function_params()
    lowercase_args(c("coloring_algorithm", "coloring_algorithm_advanced_option", "legend_loc"))

    if (verbose && do_compute_alluvial_statistics) compute_alluvial_statistics(clus_df_gather = df, graphing_columns = graphing_columns, column_weights = column_weights)

    geom_alluvium <- if (rasterise_alluvia) {
        function(...) ggrastr::rasterise(ggalluvial::geom_alluvium(...), dpi = dpi)
    } else {
        ggalluvial::geom_alluvium
    }

    if (column_weights != "value") {
        df <- df %>% dplyr::rename(value = !!sym(column_weights))
        column_weights <- "value"
    }

    df <- ggforce::gather_set_data(df, 1:2)
    if (!is.numeric(df$x)) {
        df$x <- match(as.character(df$x), graphing_columns)
    } # weird Docker issue
    df <- df[df$x == 1, ]

    if (!is.null(color_list)) {
        ditto_colors <- color_list
    } else {
        ditto_colors <- default_colors
    }

    # remove user-defined colors from available color list
    if (!(is.null(color_val))) {
        # convert named list into named vector
        if (is.list(color_val)) {
            color_val <- unlist(color_val)
        }
        ditto_colors <- ditto_colors[!(ditto_colors %in% color_val)]
    }

    match_colors <- (coloring_algorithm != "none") # we want to match color if coloring_algorithm is not none

    if (!(coloring_algorithm %in% c("none", "left", "right", "advanced", graphing_columns))) {
        stop("Invalid coloring_algorithm. Options are 'none', 'left', 'right', 'advanced', or any value in graphing_columns.")
    }

    # warning if coloring_algorithm is not none but color_boxes is False
    if (!color_boxes && !color_bands && match_colors) {
        if (verbose) message("Warning: color_boxes and color_bands are False but coloring_algorithm is specified. boxes will not be colored.")
    }

    # Extract colors for each factor, assuming ditto_colors is long enough
    if (match_colors) {
        unused_colors <- ditto_colors
        first <- TRUE
        final_colors <- c()

        if (coloring_algorithm == "left") {
            n <- 1
            for (col_group in graphing_columns) {
                num_levels <- length(levels(df[[col_group]]))
                if (first) {
                    old_colors <- unused_colors[1:num_levels]
                    names(old_colors) <- 1:num_levels
                    unused_colors <- unused_colors[!(unused_colors %in% old_colors)]
                    final_colors <- c(final_colors, rev(old_colors))
                    first <- FALSE
                } else {
                    temp_colors <- find_group2_colors(df, ditto_colors, unused_colors, old_colors,
                        group1_name = paste0("col", n - 1, "_int"), group2_name = paste0("col", n, "_int"),
                        cutoff = cutoff
                    )
                    unused_colors <- unused_colors[!(unused_colors %in% old_colors)]
                    old_colors <- temp_colors
                    names(old_colors) <- 1:num_levels
                    final_colors <- c(final_colors, rev(old_colors))
                }
                n <- n + 1
            }
        } else if (coloring_algorithm == "right") {
            n <- length(graphing_columns)
            for (col_group in rev(graphing_columns)) {
                num_levels <- length(levels(df[[col_group]]))
                if (first) {
                    old_colors <- unused_colors[1:num_levels]
                    names(old_colors) <- 1:num_levels
                    unused_colors <- unused_colors[!(unused_colors %in% old_colors)]
                    final_colors <- c(final_colors, rev(old_colors))
                    first <- FALSE
                } else {
                    temp_colors <- find_group2_colors(df, ditto_colors, unused_colors, old_colors,
                        group1_name = paste0("col", n, "_int"), group2_name = paste0("col", n - 1, "_int"),
                        cutoff = cutoff
                    )
                    unused_colors <- unused_colors[!(unused_colors %in% old_colors)]
                    old_colors <- temp_colors
                    names(old_colors) <- 1:num_levels
                    final_colors <- c(rev(old_colors), final_colors)
                    n <- n - 1
                }
            }
        } else if (coloring_algorithm == "advanced") {
            check_python_setup_with_necessary_packages(necessary_packages_for_this_step = c("igraph", "leidenalg"), additional_message = "do not set coloring_algorithm to 'advanced'", environment = environment, use_conda = use_conda)
            final_colors <- find_colors_advanced(df, ditto_colors, graphing_columns, coloring_algorithm_advanced_option = coloring_algorithm_advanced_option, resolution = resolution)
        } else {
            col_group <- coloring_algorithm
            num_levels <- length(levels(df[[col_group]]))
            old_colors <- unused_colors[1:num_levels]
            names(old_colors) <- 1:num_levels
            unused_colors <- unused_colors[!(unused_colors %in% old_colors)]
            final_colors <- c(final_colors, rev(old_colors))
            ref_colors <- rev(old_colors)
            really_final_colors <- c()
            ref_group_n <- which(coloring_algorithm == graphing_columns)[[1]]

            for (col_group in graphing_columns) {
                if (!(col_group == coloring_algorithm)) {
                    col_group_n <- which(col_group == graphing_columns)[[1]]
                    num_levels <- length(levels(df[[col_group]]))
                    temp_colors <- find_group2_colors(df, ditto_colors, unused_colors, old_colors,
                        group1_name = paste0("col", ref_group_n, "_int"), group2_name = paste0("col", col_group_n, "_int"),
                        cutoff = cutoff
                    )
                    unused_colors <- unused_colors[!(unused_colors %in% temp_colors)]
                    final_colors <- c(final_colors, rev(temp_colors))
                    really_final_colors <- c(really_final_colors, rev(temp_colors))
                } else {
                    really_final_colors <- c(really_final_colors, ref_colors)
                }
            }
            final_colors <- really_final_colors
        }
    } else {
        remaining_colors <- ditto_colors
        final_colors <- c()
        for (col_group in graphing_columns) {
            num_levels <- length(levels(df[[col_group]]))
            if (length(remaining_colors) < num_levels) {
                if (verbose) message("Warning: Some colors will be recycled.")
                remaining_colors <- ditto_colors
            }
            old_colors <- remaining_colors[1:num_levels]
            final_colors <- c(final_colors, rev(old_colors))
            remaining_colors <- remaining_colors[(1 + num_levels):length(remaining_colors)]
        }
    }

    if (!(is.null(color_val))) {
        final_value_order <- c()
        for (col_int in seq_along(graphing_columns)) {
            int_name <- paste0("col", col_int, "_int")
            group_name <- graphing_columns[[col_int]]

            curr_label <- as.character(unique(df[order(df[[int_name]]), ][[group_name]]))

            final_value_order <- c(final_value_order, rev(curr_label))
        }
        names(final_colors) <- final_value_order

        for (box_val in names(color_val)) {
            if (box_val %in% names(final_colors)) {
                box_val_color <- color_val[names(color_val) == box_val][1]
                # find where value is in final colors
                val_match <- which(box_val == names(final_colors))

                # JR added
                to_change <- c()
                for (vm in val_match) {
                    matched_color <- final_colors[vm]
                    maybe_to_change <- which(matched_color == final_colors)

                    for (matched in maybe_to_change) {
                        testing <- names(final_colors)[matched]
                        if (!(testing %in% names(color_val)) || (testing == box_val)) {
                            to_change <- c(to_change, matched)
                        }
                    }
                }
                # JR added

                # JR commented out
                # matched_color <- final_colors[val_match[1]]
                # maybe_to_change <- which(matched_color == final_colors)
                # to_change <- c()
                # for (matched in maybe_to_change) {
                #     testing <- names(final_colors)[matched]
                #     if (!(testing %in% names(color_val)) | (testing == box_val)) {
                #         to_change <- c(to_change, matched)
                #     }
                # }
                # JR commented out

                final_colors[as.integer(to_change)] <- box_val_color
            }
        }
    }

    remaining_colors <- ditto_colors[!(ditto_colors %in% final_colors)]

    final_colors_legend <- final_colors
    # generate label names
    final_label_names <- c()
    for (col_int in seq_along(graphing_columns)) {
        int_name <- paste0("col", col_int, "_int")
        group_name <- graphing_columns[[col_int]]

        curr_label <- as.character(unique(df[order(df[[int_name]]), ][[group_name]]))

        final_label_names <- c(final_label_names, rev(curr_label))
    }
    names(final_colors_legend) <- final_label_names

    # if color_bands, add to named list

    if (!is.null(color_band_list)) {
        final_colors_legend <- c(color_band_list, final_colors_legend)
    } else {
        if (!is.null(color_band_column)) {
            if (!(color_band_column %in% graphing_columns)) {
                num_levels <- length(unique(df[[color_band_column]]))
                color_band_list <- remaining_colors[1:num_levels]
                names(color_band_list) <- unique(df[[color_band_column]])
                final_colors_legend <- c(color_band_list, final_colors_legend)
            }
        }
    }
    # remove duplicate names
    final_colors_legend <- final_colors_legend[!duplicated(names(final_colors_legend))]
    # remove duplicate dims
    temp_df <- df # [1:as.integer(dim(df)[1]/2),1:dim(df)[2]]

    # uncomment to attempt mapping
    p <- ggplot(data = temp_df, aes(y = value), )
    for (x in seq_along(graphing_columns)) {
        p$mapping[[paste0("axis", x)]] <- sym(paste0("col", x, "_int"))
    }

    if (color_bands) {
        if (is.null(color_band_column)) {
            color_band_column <- graphing_columns[1]
        }
        if (color_band_boundary) {
            p <- p +
                geom_alluvium(aes(fill = !!sym(color_band_column), color = !!sym(color_band_column)),
                    alpha = alluvial_alpha
                )
        } else {
            p <- p +
                geom_alluvium(aes(fill = !!sym(color_band_column)), alpha = alluvial_alpha)
        }
        p <- p + scale_fill_manual(values = final_colors_legend)
    } else {
        if (color_band_boundary) {
            p <- p + geom_alluvium(color = "grey2", alpha = alluvial_alpha)
        } else {
            p <- p + geom_alluvium(alpha = alluvial_alpha)
        }
    }

    if (color_band_boundary) {
        if (add_legend) {
            p <- p + scale_color_manual(values = final_colors_legend, guide = "none")
        } else {
            p <- p + scale_color_manual(values = color_band_list, guide = "none")
        }
    }


    if (color_boxes) {
        if (add_legend) {
            p <- p + geom_stratum(width = box_width, aes(fill = after_stat(!!sym("final_label_names"))), linewidth = box_line_width) + scale_fill_manual(values = final_colors_legend)
        } else {
            p <- p + geom_stratum(width = box_width, fill = final_colors, linewidth = box_line_width)
        }
    } else {
        p <- p + geom_stratum(width = box_width, linewidth = box_line_width)
    }


    if (!(include_labels_in_boxes == FALSE)) {
        if (auto_adjust_text) {
            p <- p +
                geom_fit_text(
                    reflow = TRUE, stat = StatStratum, width = text_width, min.size = min_text, size = text_size,
                    aes(label = after_stat(final_label_names))
                )
        } else {
            p <- p +
                geom_text(stat = StatStratum, aes(label = after_stat(final_label_names)), size = text_size)
        }
    }

    top_y <- 0
    for (test_x in unique(df$x)) {
        curr_y <- df %>%
            filter(x == test_x) %>%
            group_by(y) %>%
            summarise(total = sum(value), .groups = "drop") %>%
            arrange(desc(total)) %>%
            mutate(cum_y = cumsum(total)) %>%
            pull(cum_y) %>%
            max()
        top_y <- max(curr_y, top_y)
    } # top_y1 and top_y2 are probably the same

    if (include_axis_titles) {
        p <- p +
            scale_x_continuous(
                breaks = 1:length(graphing_columns), labels = graphing_columns,
                position = "top"
            )
    }

    if (include_group_sizes) {
        offset_below <- top_y * 0.075
        x <- 1
        for (col_group in graphing_columns) {
            p <- p +
                annotate("text", x = x, y = -offset_below, label = length(levels(df[[col_group]])), hjust = 0.5, size = 5) # Adjust x, y for Scanpy
            x <- x + 1
        }
    }

    p <- p +
        theme_void() +
        theme(
            text = element_text(family = "sans"),
            axis.ticks.x = element_blank()
        )

    if (add_legend) {
        p <- p + theme(legend.position = legend_loc) + labs(fill = "")
    } else {
        p <- p + theme(legend.position = "none")
    }

    if (include_axis_titles) {
        p <- p + theme(axis.text.x = element_text(size = axis_text_size, vjust = axis_text_vjust)) # vjust adjusts the vertical position of column titles)
    } else {
        p <- p + theme(axis.text.x = element_blank())
    }

    if (keep_y_labels) {
        p <- p + theme(axis.text.y = element_text(size = axis_text_size))
    }

    if (!is.null(output_plot_path)) {
        if (verbose) message(sprintf("Saving plot to=%s", output_plot_path))
        ggsave(output_plot_path,
            plot = p,
            height = save_height, width = save_width, dpi = dpi, bg = "white"
        )
    }

    return(p)
}


add_int_columns <- function(df, graphing_columns, default_sorting = "alphabetical") {
    n <- 1
    for (col in graphing_columns) {
        col_int_name <- paste0("col", n, "_int")
        n <- n + 1

        # factorize input columns
        if (default_sorting == "alphabetical") {
            df[[col]] <- factor(as.character(df[[col]]), levels = sort(unique(as.character(df[[col]])), method = "radix"))
        } else if (default_sorting == "reverse_alphabetical") {
            df[[col]] <- factor(as.character(df[[col]]), levels = sort(unique(as.character(df[[col]])), method = "radix", decreasing = TRUE))
        } else if (default_sorting == "increasing") {
            df[[col]] <- factor(df[[col]], levels = names(sort(table(df[[col]]), decreasing = FALSE)))
        } else if (default_sorting == "decreasing") {
            df[[col]] <- factor(df[[col]], levels = names(sort(table(df[[col]]), decreasing = TRUE)))
        } else if (default_sorting == "random") {
            df[[col]] <- factor(df[[col]], levels = sample(unique(as.character(df[[col]]))))
        } else if (default_sorting == "fixed") {
            df[[col]] <- factor(df[[col]], levels = unique(as.character(df[[col]])))
        } else {
            stop(sprintf("default_sorting '%s' is not recognized. Please choose from 'alphabetical' (default), 'reverse_alphabetical', 'increasing', 'decreasing', or 'random'.", default_sorting))
        }

        if (!(col_int_name %in% colnames(df))) {
            # make columns integer for sorting
            df[[col_int_name]] <- as.integer(df[[col]])
        }
    }
    return(df)
}

get_alluvial_df <- function(df, column_weights = "value", do_gather_set_data = FALSE) {
    if (is.null(column_weights)) {
        column_weights <- "value"
    }
    # Convert numeric clustering columns to ordered factors
    df <- df |>
        dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
        dplyr::group_by_all() |>
        dplyr::count(name = column_weights)
    if (do_gather_set_data) {
        df <- ggforce::gather_set_data(df, 1:2)
    }
    return(df)
}

load_in_df <- function(df, graphing_columns = NULL, column_weights = NULL) {
    if (is.character(df) && grepl("\\.csv$", df)) {
        df <- read.csv(df) # load in CSV as dataframe
    } else if (tibble::is_tibble(df)) {
        df <- as.data.frame(df) # convert tibble to dataframe
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

randomly_map_int_columns <- function(clus_df_gather) {
    cols_to_shuffle <- grep("^col\\d+_int$", names(clus_df_gather), value = TRUE)

    for (col in cols_to_shuffle) {
        old_vals <- unique(clus_df_gather[[col]])
        new_vals <- sample(old_vals)
        clus_df_gather[[col]] <- match(clus_df_gather[[col]], old_vals)
        clus_df_gather[[col]] <- new_vals[clus_df_gather[[col]]]
    }

    return(clus_df_gather)
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
#' @param default_sorting Character. Default column sorting in data_preprocess if integer columns do not exist. Will not affect output if \code{sorting_algorithm == 'neighbornet'}. Options are 'alphabetical' (default), 'reverse_alphabetical', 'increasing', 'decreasing', 'random'.
#' @param output_df_path Optional character. Output path for the output data frame, in CSV format. If \code{NULL}, then will not be saved.
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param print_params Logical. If TRUE, will print function params.
#' @param load_df Internal flag; not recommended to modify.
#' @param do_gather_set_data Internal flag; not recommended to modify.
#' @param color_band_column Internal flag; not recommended to modify.
#'
#' @return A data frame where each row represents a combination of groupings, each column from \code{graphing_columns} represents a grouping, and the column \code{column_weights} ('value' if \code{column_weights} == NULL) represents the number of entities in that combination of groupings. For each column in \code{graphing_columns}, there will be an additional column \code{col1_int}, \code{col2_int}, etc. where each column corresponds to a position mapping of groupings in the respective entry of \code{graphing_columns} - for example, \code{col1_int} corresponds to \code{graphing_columns[1]}, \code{col2_int} corresponds to \code{graphing_columns[2]}, etc.
#'
#' @examples
#' # Example 1: df format 1
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- data_preprocess(df, graphing_columns = c("method1", "method2"))
#'
#' # Example 2: df format 2
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- df |>
#'     dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
#'     dplyr::group_by_all() |>
#'     dplyr::count(name = "value")
#' clus_df_gather <- data_preprocess(
#'     clus_df_gather,
#'     graphing_columns = c("method1", "method2"),
#'     column_weights = "value"
#' )
#'
#' @export
data_preprocess <- function(df, graphing_columns, column_weights = NULL, default_sorting = "alphabetical", output_df_path = NULL, verbose = FALSE, print_params = FALSE, load_df = TRUE, do_gather_set_data = FALSE, color_band_column = NULL) {
    if (print_params) print_function_params()
    lowercase_args(c("default_sorting"))

    if (load_df) {
        df <- load_in_df(df = df, graphing_columns = graphing_columns, column_weights = column_weights)
    }

    # remove all columns outside of graphing_columns, column_weights, and int columns
    cols_to_keep <- c(
        graphing_columns,
        column_weights,
        color_band_column,
        grep("^col[0-9]+_int$", names(df), value = TRUE)
    )
    cols_to_keep <- intersect(cols_to_keep, names(df))
    df <- df[, cols_to_keep, drop = FALSE]

    # # # Fill in NA values with "Missing"
    # # df[is.na(df)] <- "Missing"
    # df <- df %>%
    #     mutate(across(all_of(graphing_columns), ~ {
    #         if (is.factor(.x) || is.character(.x)) {
    #             replace(as.character(.x), is.na(.x), "Missing")
    #         } else {
    #             .x  # leave numeric and other types untouched
    #         }
    #     }))
    for (col in graphing_columns) {
        if (col %in% colnames(df)) {
            if (is.factor(df[[col]]) || is.character(df[[col]])) {
                df[[col]] <- replace(as.character(df[[col]]), is.na(df[[col]]), "Missing")
            }
        }
    }

    df <- add_int_columns(df, graphing_columns = graphing_columns, default_sorting = default_sorting)

    # sort columns according to graphing_columns
    # df <- df %>% dplyr::relocate(all_of(graphing_columns))  # put graphing_columns in front
    if (!all(intersect(colnames(df), graphing_columns) == graphing_columns)) {
        df <- reorder_and_rename_columns(df, graphing_columns)
    }

    if (is.null(column_weights) || !(column_weights %in% colnames(df))) {
        clus_df_gather <- get_alluvial_df(df, column_weights = column_weights, do_gather_set_data = do_gather_set_data)
    } else {
        clus_df_gather <- df
    }

    if ((is.character(output_df_path) && grepl("\\.csv$", output_df_path, ignore.case = TRUE))) {
        if (verbose) message(sprintf("Saving dataframe to=%s", output_df_path))
        write.csv(clus_df_gather, file = output_df_path, row.names = FALSE)
    }

    return(clus_df_gather)
}



sort_neighbornet <- function(clus_df_gather, graphing_columns = NULL, column_weights = "value", optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE, matrix_initialization_value = 1e6, same_side_matrix_initialization_value = 1e6, weight_scalar = 5e5, matrix_initialization_value_column_order = 1e6, weight_scalar_column_order = 1, column_sorting_metric = "edge_crossing", sorting_algorithm = "neighbornet", column_sorting_algorithm = "tsp", cycle_start_positions = NULL, verbose = FALSE, make_intermediate_neighbornet_plots = FALSE, environment = "wompwomp_env", use_conda = TRUE) {
    if (sorting_algorithm == "neighbornet" || column_sorting_algorithm == "neighbornet") {
        check_python_setup_with_necessary_packages(necessary_packages_for_this_step = c("splitspy", "numpy"), additional_message = "do not set sorting_algorithm or column_sorting_algorithm to 'neighbornet'", environment = environment, use_conda = use_conda)
    }
    if (verbose) message("Running neighbornet")
    cycle <- run_neighbornet(clus_df_gather, graphing_columns = graphing_columns, column_weights = column_weights, matrix_initialization_value = matrix_initialization_value, same_side_matrix_initialization_value = same_side_matrix_initialization_value, weight_scalar = weight_scalar, sorting_algorithm = sorting_algorithm, verbose = verbose)
    if (verbose) message("Cycle: ", paste(cycle, collapse = ", "))
    if (verbose) message("Determining optimal cycle start")
    res <- determine_optimal_cycle_start(clus_df_gather, cycle, graphing_columns = graphing_columns, column_weights = column_weights, optimize_column_order = optimize_column_order, optimize_column_order_per_cycle = optimize_column_order_per_cycle, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_sorting_metric = column_sorting_metric, column_sorting_algorithm = column_sorting_algorithm, cycle_start_positions = cycle_start_positions, verbose = verbose, make_intermediate_neighbornet_plots = make_intermediate_neighbornet_plots, environment = environment, use_conda = use_conda)
    clus_df_gather_neighbornet <- res$clus_df_gather
    # graphing_columns_neighbornet <- res$graphing_columns
    if (verbose) message(sprintf("crossing edges objective = %s", res$neighbornet_objective))
    return(clus_df_gather_neighbornet)
}

sort_greedy_wolf <- function(clus_df_gather, graphing_columns = NULL, column1 = NULL, column2 = NULL, fixed_column = NULL, random_initializations = 1, column_weights = "value", sorting_algorithm = "greedy_wblf", verbose = FALSE, environment = "wompwomp_env", use_conda = TRUE) {
    if (length(graphing_columns) != 2) {
        stop(sprintf("graphing_columns must be of length 2 for greedy_wblf/greedy_wolf"))
    }

    if (is.null(fixed_column)) {
        fixed_column <- column1
    } else if ((is.integer(fixed_column) || (is.double(fixed_column)))) {
        if (fixed_column > length(colnames(clus_df_gather))) {
            stop(sprintf("fixed_column index '%s' is not a column in the dataframe.", fixed_column))
        } else {
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


    length_clus_df_gather_original <- nrow(clus_df_gather)
    if (length(setdiff(c("id", "x", "y"), colnames(clus_df_gather))) > 0) {
        clus_df_gather <- ggforce::gather_set_data(clus_df_gather, 1:2)
    }

    for (i in seq_len(random_initializations)) {
        clus_df_gather_tmp <- clus_df_gather
        if (sorting_algorithm == "greedy_wblf") {
            # randomize clus_df_gather order
            for (column_num in c("col1_int", "col2_int")) {
                clus_df_gather_tmp[[column_num]] <- as.factor(clus_df_gather_tmp[[column_num]])
                clus_df_gather_tmp[[column_num]] <- factor(clus_df_gather_tmp[[column_num]], levels = sample(levels(clus_df_gather_tmp[[column_num]])))
                # clus_df_gather_tmp[[column_num]] = as.integer(clus_df_gather_tmp[[column_num]])
            }
            # WBLF
            clus_df_gather_tmp <- sort_clusters_by_agreement(clus_df_gather_tmp,
                stable_column = "col1_int",
                reordered_column = "col2_int"
            )
            clus_df_gather_tmp <- sort_clusters_by_agreement(clus_df_gather_tmp,
                stable_column = "col2_int",
                reordered_column = "col1_int"
            )
        } else if (sorting_algorithm == "greedy_wolf") {
            # randomize clus_df_gather order
            column_num <- reordered_column
            clus_df_gather_tmp[[column_num]] <- as.factor(clus_df_gather_tmp[[column_num]])
            clus_df_gather_tmp[[column_num]] <- factor(clus_df_gather_tmp[[column_num]], levels = sample(levels(clus_df_gather_tmp[[column_num]])))
            # clus_df_gather_tmp[[column_num]] = as.integer(clus_df_gather_tmp[[column_num]])
            # WOLF

            clus_df_gather_tmp <- sort_clusters_by_agreement(clus_df_gather_tmp,
                stable_column = fixed_column,
                reordered_column = reordered_column
            )
        }
        if (random_initializations > 1) {
            crossing_edges_objective <- determine_crossing_edges(clus_df_gather_tmp,
                column1 = column1, column2 = column2,
                column_weights = column_weights, load_df = FALSE, preprocess_data = FALSE, # verbose = verbose,
                output_df_path = NULL, return_weighted_layer_free_objective = TRUE, environment = environment, use_conda = use_conda
            )
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
            slice(1:(n() %/% 2)) %>% # keep first half of rows
            select(-id, -x, -y) # drop columns
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
#' @param sorting_algorithm Character. Algorithm with which to sort the values in the dataframe. Can choose from {'neighbornet', 'tsp', 'greedy_wolf', 'greedy_wblf', 'none'}. 'neighbornet' performs sorting with NeighborNet (Bryant and Moulton, 2004). 'tsp' performs Traveling Salesman Problem solver from the TSP package. greedy_wolf' implements a custom greedy algorithm where one layer is fixed, and the other layer is sorted such that each node is positioned as close to its largest parent from the fixed side as possible in a greedy fashion. 'greedy_wblf' implements the 'greedy_wolf' algorithm described previously twice, treating each column as fixed in one iteration and free in the other iteration. 'greedy_wolf' and 'greedy_wblf' are only valid when \code{graphing_columns} has exactly two entries. 'random' randomly maps blocks. 'none' keeps the mappings as-is when passed into the function.
#' @param optimize_column_order Logical. If TRUE, will optimize the order of \code{graphing_columns} to minimize edge overlap. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{length(graphing_columns) > 2}.
#' @param optimize_column_order_per_cycle Logical. If TRUE, will optimize the order of \code{graphing_columns} to minimize edge overlap upon each cycle. If FALSE, will optimize the order of \code{graphing_columns} to minimize edge overlap on the beginning cycle only. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{length(graphing_columns) > 2}.
#' @param matrix_initialization_value Positive integer. Initialized value in distance matrix for nodes in different layers without a shared edge/path. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param same_side_matrix_initialization_value Positive integer. Initialized value in distance matrix for nodes in the same layer. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param weight_scalar Positive integer. Scalar with which to multiply edge weights after taking their -log in the distance matrix for nodes with a nonzero edge. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param matrix_initialization_value_column_order Positive integer. Initialized value in distance matrix for optimizing column order. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param weight_scalar_column_order Positive integer. Scalar with which to loss function after taking their log1p in the distance matrix for optimizing column order. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param column_sorting_metric Character. Metric to use for determining column order. Options are "edge_crossing" (default) or "ARI". Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param column_sorting_algorithm Character. Algorithm to use for determining column order. Options are "tsp" (default) or "neighbornet". Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param cycle_start_positions Set. Cycle start positions to consider. Anything outside this set will be skipped. Only applies when \code{sorting_algorithm == 'neighbornet'}.
#' @param fixed_column Character or Integer. Name or position of the column in \code{graphing_columns} to keep fixed during sorting. Only applies when \code{sorting_algorithm == 'greedy_wolf'}.
#' @param random_initializations Integer. Number of random initializations for the positions of each grouping in \code{graphing_columns}. Only applies when \code{sorting_algorithm == 'greedy_wolf' or sorting_algorithm == 'greedy_wblf'}.
#' @param output_df_path Optional character. Output path for the output data frame, in CSV format. If \code{NULL}, then will not be saved.
#' @param preprocess_data Logical. If TRUE, will preprocess the data with the \code{data_preprocess} function.
#' @param default_sorting Character. Default column sorting in data_preprocess if integer columns do not exist. Will not affect output if \code{sorting_algorithm == 'neighbornet'}. Options are 'alphabetical' (default), 'reverse_alphabetical', 'increasing', 'decreasing', 'random'.
#' @param return_updated_graphing_columns Logical. If FALSE, will only return the updated data frame. If TRUE, will return both the updated data frame and the updated graphing_columns parameter in the order in which the columns should be graphed.
#' @param environment Character. Python environment (if applicable). Default: 'wompwomp_env'
#' @param use_conda Logical. Whether or not to use conda for Python (if applicable)
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param print_params Logical. If TRUE, will print function params.
#' @param load_df Internal flag; not recommended to modify.
#' @param make_intermediate_neighbornet_plots Internal flag; not recommended to modify.
#' @param do_compute_alluvial_statistics Internal flag; not recommended to modify.
#'
#' @return
#' If return_updated_graphing_columns == FALSE (default): A data frame where each row represents a combination of groupings, each column from \code{graphing_columns} represents a grouping, and the column \code{column_weights} ('value' if \code{column_weights} == NULL) represents the number of entities in that combination of groupings. For each column in \code{graphing_columns}, there will be an additional column \code{col1_int}, \code{col2_int}, etc. where each column corresponds to a position mapping of groupings in the respective entry of \code{graphing_columns} - for example, \code{col1_int} corresponds to \code{graphing_columns[1]}, \code{col2_int} corresponds to \code{graphing_columns[2]}, etc. The position mappings in these columns, as well as the order of the columns (if \code{optimize_column_order} is TRUE), will be sorted according to \code{sorting_algorithm}.
#' If return_updated_graphing_columns == TRUE: A list of the data frame described above and the sorted \code{graphing_columns}, in the keys 'clus_df_gather' and 'graphing_columns', respectively.
#'
#' @examples
#' # Example 1: df format 1
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- data_sort(
#'     df,
#'     graphing_columns = c("method1", "method2"),
#'     sorting_algorithm = "tsp",
#'     column_sorting_algorithm = "tsp"
#' )
#'
#' # Example 2: df format 2
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- df |>
#'     dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
#'     dplyr::group_by_all() |>
#'     dplyr::count(name = "value")
#' clus_df_gather <- data_sort(
#'     clus_df_gather,
#'     graphing_columns = c("method1", "method2"),
#'     column_weights = "value",
#'     sorting_algorithm = "tsp",
#'     column_sorting_algorithm = "tsp"
#' )
#'
#' @export
data_sort <- function(df, graphing_columns = NULL, column1 = NULL, column2 = NULL, column_weights = NULL, sorting_algorithm = "neighbornet", optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE, matrix_initialization_value = 1e6, same_side_matrix_initialization_value = 1e6, weight_scalar = 5e5, matrix_initialization_value_column_order = 1e6, weight_scalar_column_order = 1, column_sorting_metric = "edge_crossing", column_sorting_algorithm = "tsp", cycle_start_positions = NULL, fixed_column = NULL, random_initializations = 1, output_df_path = NULL, preprocess_data = TRUE, default_sorting = "alphabetical", return_updated_graphing_columns = FALSE, verbose = FALSE, print_params = FALSE, load_df = TRUE, make_intermediate_neighbornet_plots = FALSE, do_compute_alluvial_statistics = TRUE, environment = "wompwomp_env", use_conda = TRUE) {
    if (print_params) print_function_params()
    lowercase_args(c("sorting_algorithm", "column_sorting_metric", "column_sorting_algorithm", "default_sorting"))

    #* Type Checking Start
    valid_algorithms <- c("neighbornet", "tsp", "greedy_wolf", "greedy_wblf", "none", "random")
    if (!(sorting_algorithm %in% valid_algorithms)) {
        stop(sprintf(
            "Invalid sorting_algorithm: '%s'. Must be one of: %s",
            sorting_algorithm, paste(valid_algorithms, collapse = ", ")
        ))
    }

    if (sorting_algorithm == "neighbornet" || sorting_algorithm == "tsp") {
        for (col in graphing_columns) {
            if (grepl("~~", col)) {
                stop(sprintf("No entry of graphing_columns can contain '~~' when sorting_algorithm == neighbornet or tsp. Issue with column '%s'.", col))
            }
        }
    }

    if (((sorting_algorithm == "none") || (sorting_algorithm == "random") || (sorting_algorithm == "neighbornet") || (sorting_algorithm == "tsp")) && (random_initializations > 1)) {
        sprintf("random_initializations > 1 but sorting algorithm is %s Setting random_initializations to 1.", sorting_algorithm)
        random_initializations <- 1
    }

    if (!is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        stop("Specify either graphing_columns or column1/column2, not both.")
    }

    if (is.vector(graphing_columns) && length(graphing_columns) < 2) {
        stop("graphing_columns must have at least 2 entries.")
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
    } else { # length 2
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

    if ((is.null(fixed_column)) && (sorting_algorithm == "greedy_wolf")) {
        if (verbose) message(sprintf("Using column %s as fixed_column for greedy_wolf by default", column1))
        fixed_column <- column1
    }

    if (sorting_algorithm == "greedy_wolf") {
        default_sorting <- "fixed"
    }
    #* Type Checking End


    # Preprocess (i.e., add int columns and do the grouping)
    if (preprocess_data) {
        if (verbose) message("Preprocessing data before sorting")
        clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns, column_weights = column_weights, default_sorting = default_sorting, load_df = load_df, do_gather_set_data = FALSE)
        if (is.null(column_weights)) {
            column_weights <- "value" # is set during data_preprocess
        }
    } else {
        clus_df_gather <- df
    }

    if (verbose && do_compute_alluvial_statistics) compute_alluvial_statistics(clus_df_gather = clus_df_gather, graphing_columns = graphing_columns, column_weights = column_weights)

    if (sorting_algorithm == "neighbornet" || sorting_algorithm == "tsp") {
        # O(n^3) complexity, where n is the sum of blocks across all layers
        clus_df_gather_sorted <- sort_neighbornet(clus_df_gather = clus_df_gather, graphing_columns = graphing_columns, column_weights = column_weights, optimize_column_order = optimize_column_order, optimize_column_order_per_cycle = optimize_column_order_per_cycle, matrix_initialization_value = matrix_initialization_value, same_side_matrix_initialization_value = same_side_matrix_initialization_value, weight_scalar = weight_scalar, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_sorting_metric = column_sorting_metric, sorting_algorithm = sorting_algorithm, column_sorting_algorithm = column_sorting_algorithm, cycle_start_positions = cycle_start_positions, verbose = verbose, make_intermediate_neighbornet_plots = make_intermediate_neighbornet_plots, environment = environment, use_conda = use_conda)
    } else if (sorting_algorithm == "greedy_wblf" || sorting_algorithm == "greedy_wolf") {
        # O(n_1 * n_2) complexity, where n1 is the number of blocks in layer 1, and n2 is the number of blocks in layer 2
        clus_df_gather_sorted <- sort_greedy_wolf(clus_df_gather = clus_df_gather, graphing_columns = graphing_columns, column1 = column1, column2 = column2, column_weights = column_weights, fixed_column = fixed_column, random_initializations = random_initializations, sorting_algorithm = sorting_algorithm, verbose = verbose, environment = environment, use_conda = use_conda)
    } else if (sorting_algorithm == "random") {
        clus_df_gather_sorted <- randomly_map_int_columns(clus_df_gather)
    } else if (sorting_algorithm == "none") {
        clus_df_gather_sorted <- clus_df_gather
    } else {
        stop(sprintf("Invalid sorting_algorithm: '%s'. Must be one of: %s", sorting_algorithm, paste(valid_algorithms, collapse = ", ")))
    }

    # print objective - don't do for neighbornet because I did it right before
    if ((verbose) && (sorting_algorithm != "neighbornet")) {
        message("Determining crossing edges objective (to disable, use verbose==FALSE)")
        objective <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns = graphing_columns, column_weights = column_weights, load_df = FALSE, preprocess_data = FALSE, return_weighted_layer_free_objective = TRUE, default_sorting = default_sorting, environment = environment, use_conda = use_conda) # verbose = verbose
        message(sprintf("crossing edges objective = %s", objective))
    }

    if (verbose) message("Complete with sorting")

    # Save if desired
    if ((is.character(output_df_path) && grepl("\\.csv$", output_df_path, ignore.case = TRUE))) {
        if (verbose) message(sprintf("Saving sorted dataframe to=%s", output_df_path))
        write.csv(clus_df_gather_sorted, file = output_df_path, row.names = FALSE)
    }

    if (return_updated_graphing_columns) {
        graphing_columns_sorted <- graphing_columns[order(match(graphing_columns, names(clus_df_gather_sorted)))] # reorder graphing_columns to match any changed order in clus_df_gather_sorted
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
#' @param sorting_algorithm Character. Algorithm with which to sort the values in the dataframe. Can choose from {'neighbornet', 'tsp', 'greedy_wolf', 'greedy_wblf', 'random', 'none'}. 'neighbornet' performs sorting with NeighborNet (Bryant and Moulton, 2004). 'tsp' performs Traveling Salesman Problem solver from the TSP package. 'greedy_wolf' implements a custom greedy algorithm where one layer is fixed, and the other layer is sorted such that each node is positioned as close to its largest parent from the fixed side as possible in a greedy fashion. 'greedy_wblf' implements the 'greedy_wolf' algorithm described previously twice, treating each column as fixed in one iteration and free in the other iteration. 'greedy_wolf' and 'greedy_wblf' are only valid when \code{graphing_columns} has exactly two entries. 'random' randomly maps blocks. 'none' keeps the mappings as-is when passed into the function.
#' @param optimize_column_order Logical. If TRUE, will optimize the order of \code{graphing_columns} to minimize edge overlap. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{length(graphing_columns) > 2}.
#' @param optimize_column_order_per_cycle Logical. If TRUE, will optimize the order of \code{graphing_columns} to minimize edge overlap upon each cycle. If FALSE, will optimize the order of \code{graphing_columns} to minimize edge overlap on the beginning cycle only. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{length(graphing_columns) > 2}.
#' @param matrix_initialization_value Positive integer. Initialized value in distance matrix for nodes in different layers without a shared edge/path. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param same_side_matrix_initialization_value Positive integer. Initialized value in distance matrix for nodes in the same layer. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param weight_scalar Positive integer. Scalar with which to multiply edge weights after taking their -log in the distance matrix for nodes with a nonzero edge. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param matrix_initialization_value_column_order Positive integer. Initialized value in distance matrix for optimizing column order. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param weight_scalar_column_order Positive integer. Scalar with which to loss function after taking their log1p in the distance matrix for optimizing column order. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param column_sorting_metric Character. Metric to use for determining column order. Options are "edge_crossing" (default) or "ARI". Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param column_sorting_algorithm Character. Algorithm to use for determining column order. Options are "tsp" (default) or "neighbornet". Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param cycle_start_positions Set. Cycle start positions to consider. Anything outside this set will be skipped. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param fixed_column Character or Integer. Name or position of the column in \code{graphing_columns} to keep fixed during sorting. Only applies when \code{sorting_algorithm == 'greedy_wolf'}.
#' @param random_initializations Integer. Number of random initializations for the positions of each grouping in \code{graphing_columns}. Only applies when \code{sorting_algorithm == 'greedy_wolf' or sorting_algorithm == 'greedy_wblf'}.
#' @param color_boxes Logical. Whether to color the strata/boxes (representing groups).
#' @param color_bands Logical. Whether to color the alluvia/edges (connecting the strata).
#' @param color_list Optional named list or vector of colors to override default group colors.
#' @param color_band_list Optional named list or vector of colors to override default band colors.
#' @param color_band_column Optional Character. Which column to use for coloring bands.
#' @param color_val Optional named list where the entries are colors and the names correspond to values of the dataframe that should use those colors
#' @param color_band_boundary Logical. Whether or not to color boundaries between bands
#' @param coloring_algorithm Character. Matching colors methods. Choices are 'advanced' (default), 'none', 'left', 'right', or any value in \code{graphing_columns}.
#' @param coloring_algorithm_advanced_option Character. If \code{coloring_algorithm == 'advanced'}, then choose graph clustering algorithm. Choices are 'leiden' (default) or 'louvain'.
#' @param resolution Numeric If \code{coloring_algorithm == 'advanced'}, then choose resolution for the graph clustering algorithm. Affects coloring of both bands and boxes.
#' @param cutoff Numeric If \code{coloring_algorithm != 'none' and coloring_algorithm != 'advanced'}, sets the cutoff for color matching, below which a new color will be assigned.
#' @param alluvial_alpha Numeric between 0 and 1. Transparency level for the alluvial bands.
#' @param include_labels_in_boxes Logical. Whether to include text labels inside the rectangular group boxes.
#' @param include_axis_titles Logical. Whether to display axis titles for column1 and column2.
#' @param include_group_sizes Logical. If \code{TRUE}, includes group sizes in the labels (e.g., "Group A (42)").
#' @param output_plot_path Character. File path to save the plot (e.g., "plot.png"). If \code{NULL}, then will not be saved.
#' @param output_df_path Optional character. Output path for the output data frame, in CSV format. If \code{NULL}, then will not be saved.
#' @param preprocess_data Logical. If TRUE, will preprocess the data with the \code{data_preprocess} function.
#' @param default_sorting Character. Default column sorting in data_preprocess if integer columns do not exist. Will not affect output if \code{sorting_algorithm == 'neighbornet'}. Options are 'alphabetical' (default), 'reverse_alphabetical', 'increasing', 'decreasing', 'random'.
#' @param box_width Numeric between 0 and 1. Box width
#' @param text_width Numeric between 0 and 1. Text width
#' @param min_text Integer greater than 0. Min text
#' @param text_size Integer greater than 0. Text size (works whether auto_adjust_text is TRUE or FALSE).
#' @param auto_adjust_text Logical. Whether to automatically adjust text size to fit in box.
#' @param axis_text_size Integer greater than 0. Axis text size
#' @param axis_text_vjust Integer. Axis text vjust
#' @param save_height Integer greater than 0. Save height, in inches
#' @param save_width Integer greater than 0. Save width, in inches
#' @param dpi Integer greater than 0. DPI for \code{output_plot_path}, if \code{output_plot_path} is a raster image or \code{rasterise_alluvia} is TRUE
#' @param rasterise_alluvia Logical. Whether to rasterize the alluvia if \code{output_plot_path} is a PDF. Can save space if DPI low enough
#' @param keep_y_labels Keep y labels
#' @param box_line_width Box line width
#' @param environment Character. Python environment (if applicable). Default: 'wompwomp_env'
#' @param use_conda Logical. Whether or not to use conda for Python (if applicable)
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param print_params Logical. If TRUE, will print function params.
#' @param make_intermediate_neighbornet_plots Internal flag; not recommended to modify.
#' @param add_legend Logical. If TRUE, will generate a legend of the colors of boxes and alluvial
#' @param legend_loc Character. Location of legend. Only applies if \code{add_legened == TRUE}. Choices are 'right' (default), 'left', 'bottom', 'top'
#'
#' @return A \code{ggplot2} object representing the alluvial plot.
#'
#' @examples
#' # Example 1: df format 1
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' p <- plot_alluvial(df,
#'     graphing_columns = c("method1", "method2"),
#'     sorting_algorithm = "tsp",
#'     column_sorting_algorithm = "tsp",
#'     coloring_algorithm = "right"
#' )
#'
#' # Example 2: df format 2
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- df |>
#'     dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
#'     dplyr::group_by_all() |>
#'     dplyr::count(name = "value")
#' p <- plot_alluvial(
#'     clus_df_gather,
#'     graphing_columns = c("method1", "method2"),
#'     column_weights = "value",
#'     sorting_algorithm = "tsp",
#'     column_sorting_algorithm = "tsp",
#'     coloring_algorithm = "right"
#' )
#'
#' @export
plot_alluvial <- function(df, graphing_columns = NULL, column1 = NULL, column2 = NULL,
                          column_weights = NULL, sorting_algorithm = "neighbornet",
                          optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE,
                          matrix_initialization_value = 1e6, same_side_matrix_initialization_value = 1e6,
                          weight_scalar = 5e5, matrix_initialization_value_column_order = 1e6,
                          weight_scalar_column_order = 1, column_sorting_metric = "edge_crossing",
                          column_sorting_algorithm = "tsp", cycle_start_positions = NULL, fixed_column = NULL,
                          random_initializations = 1, color_boxes = TRUE, color_bands = FALSE,
                          color_list = NULL, color_band_list = NULL, color_band_column = NULL, color_val = NULL,
                          color_band_boundary = FALSE, coloring_algorithm = "advanced", coloring_algorithm_advanced_option = "leiden", resolution = 1, cutoff = .5,
                          alluvial_alpha = 0.5,
                          include_labels_in_boxes = TRUE, include_axis_titles = TRUE, include_group_sizes = FALSE,
                          output_plot_path = NULL, output_df_path = NULL, preprocess_data = TRUE,
                          default_sorting = "alphabetical", box_width = 1 / 3, text_width = 1 / 4, min_text = 4, text_size = 14,
                          auto_adjust_text = TRUE, axis_text_size = 2, axis_text_vjust = 0, save_height = 6, save_width = 6, dpi = 300, rasterise_alluvia = FALSE,
                          keep_y_labels = FALSE, box_line_width = 1, verbose = FALSE, print_params = FALSE,
                          make_intermediate_neighbornet_plots = FALSE, environment = "wompwomp_env", use_conda = TRUE,
                          add_legend = FALSE, legend_loc = "right") {
    if (print_params) print_function_params()
    lowercase_args(c("sorting_algorithm", "column_sorting_metric", "column_sorting_algorithm", "coloring_algorithm", "coloring_algorithm_advanced_option", "default_sorting", "legend_loc"))

    #* Type Checking Start
    # ensure someone doesn't specify both graphing_columns and column1/2
    if (!is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        stop("Specify either graphing_columns or column1/column2, not both.")
    }

    if (is.vector(graphing_columns) && length(graphing_columns) < 2) {
        stop("graphing_columns must have at least 2 entries.")
    }

    if (preprocess_data) {
        if (verbose) message("Loading in data")
        df <- load_in_df(df = df, graphing_columns = graphing_columns, column_weights = column_weights)
    }

    if (nrow(df) == 0) {
        stop("df has no rows.")
    }

    if (!is.null(graphing_columns) && any(!graphing_columns %in% colnames(df))) {
        stop("Some graphing_columns are not present in the dataframe.")
    }

    if (ncol(df) < 2) {
        stop("Dataframe must have at least 2 columns when column_weights is NULL.")
    } else if (ncol(df) > 2) {
        if (is.null(graphing_columns) && is.null(column1) && is.null(column2)) {
            stop("graphing_columns must be specified when dataframe has more than 2 columns and column_weights is NULL.")
        }
    } else { # length 2
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
    } else if ((is.integer(fixed_column) || (is.double(fixed_column)))) {
        if (fixed_column > length(colnames(df))) {
            stop(sprintf("fixed_column index '%s' is not a column in the dataframe.", fixed_column))
        } else {
            fixed_column <- colnames(df)[fixed_column]
        }
    } else if (!(fixed_column %in% colnames(df))) {
        stop(sprintf("fixed_column '%s' is not a column in the dataframe.", fixed_column))
    }

    if (!is.null(color_band_column)) {
        color_bands <- TRUE
    }
    #* Type Checking End

    if (sorting_algorithm == "greedy_wolf") {
        default_sorting <- "fixed"
    }
    # Preprocess
    if (preprocess_data) {
        if (verbose) message("Preprocessing data before sorting")
        clus_df_gather_unsorted <- data_preprocess(df = df, graphing_columns = graphing_columns, column_weights = column_weights, color_band_column = color_band_column, default_sorting = default_sorting, load_df = FALSE, do_gather_set_data = FALSE)
        if (is.null(column_weights)) {
            column_weights <- "value" # is set during data_preprocess
        }
    } else {
        clus_df_gather_unsorted <- df
    }

    do_compute_alluvial_statistics <- TRUE
    if (verbose) {
        compute_alluvial_statistics(clus_df_gather = clus_df_gather_unsorted, graphing_columns = graphing_columns, column_weights = column_weights)
        do_compute_alluvial_statistics <- FALSE
    }

    # Sort
    if (verbose) message(sprintf("Sorting data with sorting_algorithm=%s", sorting_algorithm))
    data_sort_output <- data_sort(df = clus_df_gather_unsorted, graphing_columns = graphing_columns, column_weights = column_weights, sorting_algorithm = sorting_algorithm, optimize_column_order = optimize_column_order, optimize_column_order_per_cycle = optimize_column_order_per_cycle, matrix_initialization_value = matrix_initialization_value, same_side_matrix_initialization_value = same_side_matrix_initialization_value, weight_scalar = weight_scalar, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_sorting_metric = column_sorting_metric, column_sorting_algorithm = column_sorting_algorithm, cycle_start_positions = cycle_start_positions, fixed_column = fixed_column, random_initializations = random_initializations, output_df_path = output_df_path, return_updated_graphing_columns = TRUE, preprocess_data = FALSE, load_df = FALSE, verbose = verbose, make_intermediate_neighbornet_plots = make_intermediate_neighbornet_plots, do_compute_alluvial_statistics = do_compute_alluvial_statistics, environment = environment, use_conda = use_conda)
    clus_df_gather <- data_sort_output$clus_df_gather
    graphing_columns <- data_sort_output$graphing_columns

    # Plot
    if (verbose) message("Plotting data")
    alluvial_plot <- plot_alluvial_internal(clus_df_gather,
        graphing_columns = graphing_columns, column_weights = column_weights,
        color_list = color_list, color_boxes = color_boxes,
        color_bands = color_bands, color_band_list = color_band_list,
        color_band_column = color_band_column, color_band_boundary = color_band_boundary,
        coloring_algorithm = coloring_algorithm, cutoff = cutoff, color_val = color_val,
        alluvial_alpha = alluvial_alpha,
        include_labels_in_boxes = include_labels_in_boxes, include_axis_titles = include_axis_titles,
        include_group_sizes = include_group_sizes,
        output_plot_path = output_plot_path, verbose = verbose,
        box_width = box_width, text_width = text_width, min_text = min_text, text_size = text_size,
        auto_adjust_text = auto_adjust_text, axis_text_size = axis_text_size, axis_text_vjust = axis_text_vjust,
        save_height = save_height, save_width = save_width, dpi = dpi, rasterise_alluvia = rasterise_alluvia,
        keep_y_labels = keep_y_labels, box_line_width = box_line_width, do_compute_alluvial_statistics = do_compute_alluvial_statistics,
        coloring_algorithm_advanced_option = coloring_algorithm_advanced_option, resolution = resolution, environment = environment, use_conda = use_conda,
        add_legend = add_legend, legend_loc = legend_loc
    )

    return(alluvial_plot)
}
