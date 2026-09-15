#' wompwomp: Cluster-matching alluvial plots
#'
#' @name wompwomp-imports
#' @rdname wompwomp
#' @importFrom dplyr mutate select group_by summarise desc ungroup slice n pull across all_of arrange
#' @importFrom igraph V cluster_louvain cluster_leiden E
#' @importFrom tibble is_tibble
#' @importFrom utils read.csv write.csv combn
#' @importFrom stats setNames
#' @importFrom rlang sym .data
#' @importFrom tidyselect eval_select

utils::globalVariables(c(
    ".data", ":=", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "stratum", "total", "cum_y", "best_cluster_agreement", "neighbor_net", "alluvium", "pos", "count", "group1", "group2", "value", "group1_size", "group2_size", "weight", "parent", "group_name",
    "default_sorting", "print_params", "preprocess_data", "do_compute_alluvial_statistics",
    "optimize_column_order_per_cycle", "weight_scalar", "matrix_initialization_value", "same_side_matrix_initialization_value",
    "matrix_initialization_value_column_order", "weight_scalar_column_order", "column_metric",
    "cycle_start_positions", "weighted_metric", "valid_algorithms"
))

compute_alluvial_statistics <- function(clus_df_gather, cols, wt = "value") {
    message(sprintf("Alluvial statistics: n = number of elements; m = number of graphing columns; a = number of alluvia/edges; k_i = number of blocks in layer i (where i goes from 1:m); K_sum = number of blocks across all layers; K_prod = product of blocks across all layers"))
    message(sprintf("n = %s", sum(clus_df_gather[[wt]], na.rm = TRUE)))
    message(sprintf("m = %s", length(cols)))
    message(sprintf("a = %s", nrow(clus_df_gather)))
    
    K_sum <- 0
    K_prod <- 1
    for (i in 1:length(cols)) {
        colname <- paste0("col", i, "_int")
        k_i <- length(unique(clus_df_gather[[colname]]))
        message(sprintf("k_%s = %s", i, k_i))
        
        K_sum <- K_sum + k_i
        K_prod <- K_prod * k_i
    }
    
    message(sprintf("K_sum = %s", K_sum))
    message(sprintf("K_prod = %s", K_prod))
}

determine_column_order <- function(clus_df_gather_neighbornet, cols, wt = "value", matrix_initialization_value_column_order = 1e6, weight_scalar_column_order = 1, column_metric = "edge_crossing", column_method = "tsp", verbose = FALSE, weighted_metric = TRUE) {
    # sort_to_uncross_options() match.arg()s column_metric to lowercase, so
    # accept either spelling of "ari" rather than only the uppercase one.
    column_metric <- tolower(column_metric)
    if (column_method == "none") {
        return(cols)
    } else if (column_method == "random") {
        return(sample(cols))
    }
    
    # this doesn't strictly need its own condition (2 choose 2 is 1 anyways), but does avoid a little overhead
    if (length(cols) == 2) {
        return(cols)
    }
    
    if (column_metric == "ari") {
        if (!requireNamespace("mclust", quietly = TRUE)) {
            stop("The 'mclust' package is required to compute Adjusted Rand Index (ARI) with column_metric == 'ARI'. Please install it with install.packages('mclust').")
        }
    }
    
    column_dist_matrix <- matrix(matrix_initialization_value_column_order,
                                 nrow = length(cols), ncol = length(cols),
                                 dimnames = list(cols, cols)
    )
    
    pairs <- combn(cols, 2)
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
        cols_to_keep <- c(column1, column2, col1_int, col2_int, wt)
        clus_df_gather_neighbornet_tmp <- clus_df_gather_neighbornet[, cols_to_keep]
        
        # Step 4: Rename columns col{integer1}_int → col1_int, col{integer2}_int → col2_int
        names(clus_df_gather_neighbornet_tmp)[match(c(col1_int, col2_int), names(clus_df_gather_neighbornet_tmp))] <- c("col1_int", "col2_int")
        graphing_columns_tmp <- c(column1, column2)
        
        if (column_metric == "ari") {
            # NOTE: rep() coerces the weights to integer, so non-integer weights
            # are silently truncated when expanding rows for the ARI.
            expanded_df <- clus_df_gather_neighbornet_tmp[rep(seq_len(nrow(clus_df_gather_neighbornet_tmp)), clus_df_gather_neighbornet_tmp[[wt]]), ]
            neighbornet_objective <- mclust::adjustedRandIndex(expanded_df$col1_int, expanded_df$col2_int)
            # ARI in [-0.5, 1]; map to a distance in [0, 1.5] where perfect
            # agreement (ARI = 1) is distance 0.
            neighbornet_objective <- -neighbornet_objective + 1
            # Rescale to a range comparable with the log1p edge-crossing metric so
            # one TSP tolerance works for both.
            ARI_DISTANCE_SCALE <- 50
            neighbornet_objective <- weight_scalar_column_order * ARI_DISTANCE_SCALE * neighbornet_objective
        } else if (column_metric == "edge_crossing") {
            neighbornet_objective <- compute_crossing_objective(
                clus_df_gather_neighbornet_tmp,
                cols = graphing_columns_tmp,
                wt = wt,
                weighted_metric = weighted_metric
            )$output_objective
            neighbornet_objective <- weight_scalar_column_order * log1p(neighbornet_objective) # log1p to avoid issue of log(0)
        } else {
            stop(sprintf("column_metric '%s' is not a valid option.", column_metric))
        }
        
        column_dist_matrix[column1, column2] <- neighbornet_objective
        column_dist_matrix[column2, column1] <- neighbornet_objective
    }
    labels <- cols # assuming this is a character vector
    if (verbose) message(sprintf("Running '%s' for column order", column_method))
    
    if (column_method == "tsp") {
        tsp_instance <- TSP::TSP(column_dist_matrix)
        tour <- TSP::solve_TSP(tsp_instance)
        cycle <- as.integer(tour)
    } else if (column_method == "neighbornet") {
        cycle <- neighbor_net_cycle(labels, column_dist_matrix)
    } else {
        stop(sprintf("column_method '%s' is not a valid option.", column_method))
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

run_neighbornet <- function(data, cols, wt = "value", matrix_initialization_value = 1e6, same_side_matrix_initialization_value = 1e6, weight_scalar = 5e5, method = "neighbornet", verbose = FALSE, fixed_orders = list(), fixed_order_bias = 10) {
    # map from string to int if needed
    if (is.null(wt) || length(wt) == 0 || !(wt %in% colnames(data))) {
        clus_df_gather <- get_alluvial_df(data, wt = wt)
    } else {
        clus_df_gather <- data
    }
    
    
    # Add prefixes to distinguish node types
    
    # prefix is "tissue~~"
    for (col in cols) {
        clus_df_gather[[col]] <- paste0(col, "~~", clus_df_gather[[col]])
    }
    
    # prefix is "column1_"
    # for (i in seq_along(cols)) {
    #     col <- cols[i]
    #     clus_df_gather[[col]] <- paste0("column", i, "_", clus_df_gather[[col]])
    # }
    
    # Get all node names
    all_nodes <- sort(unique(unlist(clus_df_gather[cols])))
    
    # Compute full distance matrix based on -log(edge weight)
    # Initialize distance matrix
    full_dist_matrix <- matrix(matrix_initialization_value,
                               nrow = length(all_nodes), ncol = length(all_nodes),
                               dimnames = list(all_nodes, all_nodes)
    )
    
    # Same-side initialization
    if (same_side_matrix_initialization_value != matrix_initialization_value) {
        for (col in cols) {
            prefix <- paste0(col, "~~")
            node_indices <- which(startsWith(all_nodes, prefix))
            full_dist_matrix[node_indices, node_indices] <- same_side_matrix_initialization_value
        }
    }
    
    # Within a fixed axis, grow the same-side distance with the gap between two
    # blocks' fixed positions, so the cycle tends to visit them in that order.
    # A ramp of 10 * weight_scalar gave the fewest crossings on synthetic
    # clustered data; weaker ramps are ignored and stronger ones distort the cycle.
    for (col in names(fixed_orders)) {
        fixed_nodes <- paste0(col, "~~", fixed_orders[[col]])
        n_fixed <- length(fixed_nodes)
        if (n_fixed > 1 && fixed_order_bias > 0) {
            gaps <- abs(outer(seq_len(n_fixed), seq_len(n_fixed), "-")) / (n_fixed - 1)
            full_dist_matrix[fixed_nodes, fixed_nodes] <- same_side_matrix_initialization_value + fixed_order_bias * weight_scalar * gaps
        }
    }

    # Get all 2-column combinations
    pairwise_groupings <- combn(cols, 2, simplify = FALSE)

    # For each combination, group and summarize, then fill the (symmetric)
    # distance matrix for that pair in one vectorized assignment instead of a
    # per-row strsplit()+scalar lookup loop.
    for (columns in pairwise_groupings) {
        summarized <- clus_df_gather |>
            group_by(across(all_of(columns))) |>
            summarise(total_value = sum(!!sym(wt)), .groups = "drop")

        n1 <- as.character(summarized[[columns[1]]])
        n2 <- as.character(summarized[[columns[2]]])
        w <- summarized$total_value

        valid <- w > 0
        n1 <- n1[valid]
        n2 <- n2[valid]
        vals <- weight_scalar * -log(w[valid])

        full_dist_matrix[cbind(n1, n2)] <- vals
        full_dist_matrix[cbind(n2, n1)] <- vals # symmetric since graph is undirected
    }
    
    # Translate so every entry is positive, but only for `tsp`.
    # TSP::solve_TSP does not terminate on a matrix containing negative
    # entries -- measured 0/8 completions within 8s on matrices with
    # negatives vs 8/8 on the same matrices after translating -- and a TSP
    # tour visits the same number of edges whichever order it takes, so
    # adding a constant leaves the optimal tour unchanged.
    # NeighborNet accepts negative distances, so it is left untranslated.
    if (method == "tsp") {
        min_val_abs <- abs(min(full_dist_matrix))
        full_dist_matrix <- full_dist_matrix + (min_val_abs + 1)
    }
    
    labels <- all_nodes # assuming this is a character vector
    mat <- full_dist_matrix
    mat[is.infinite(mat)] <- 1e6
    mat[is.na(mat)] <- 1e6

    if (verbose) message(sprintf("Running '%s' for stratum order", method))
    
    if (method == "tsp") {
        tsp_instance <- TSP::TSP(mat)
        tour <- TSP::solve_TSP(tsp_instance)
        cycle <- as.integer(tour)
    } else if (method == "neighbornet") {
        cycle <- neighbor_net_cycle(labels, mat)
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

# swap_graphing_column_order_based_on_graphing_column_int_order <- function(cols, graphing_columns_int) {
#     # Get the index of each graphing_columns_int entry (e.g., "col2_int" → 2)
#     int_positions <- as.integer(gsub("col([0-9]+)_int", "\\1", graphing_columns_int))
#
#     # Create an empty character vector of the correct length
#     reordered_graphing_columns <- character(length(cols))
#
#     # Place each graphing column at its new position
#     reordered_graphing_columns[int_positions] <- cols
#
#     return(reordered_graphing_columns)
# }

# # example:
# cols          <- c("tissue", "sex", "cluster")
# graphing_columns_int      <- c("col2_int", "col3_int", "col1_int")
# output: c("sex", "cluster", "tissue")
swap_graphing_column_order_based_on_graphing_column_int_order <- function(cols, graphing_columns_int) {
    # Extract suffixes to determine new order
    suffixes <- as.integer(gsub("col([0-9]+)_int", "\\1", graphing_columns_int))
    
    # Match suffix to index in original cols (which are col1 = cols[1], etc.)
    return(cols[suffixes])
}



determine_optimal_cycle_start <- function(data, cycle, cols = NULL, wt = "value", optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE, matrix_initialization_value_column_order = 1e6, weight_scalar_column_order = 1, column_metric = "edge_crossing", column_method = "tsp", cycle_start_positions = NULL, verbose = FALSE, weighted_metric = TRUE, fixed_orders = list()) {
    # # Commented out because I'm not sold on ARI
    # if (optimize_column_order_per_cycle && (column_metric == "edge_crossing")) {
    #     if (verbose) message("column_metric == 'edge_crossing' and optimize_column_order_per_cycle is TRUE. This might be a bit slow. Consider setting column_metric == 'ARI' and/or optimize_column_order_per_cycle to FALSE.")
    # }
    
    # factorize input columns
    for (col in cols) {
        data[[col]] <- as.factor(as.character(data[[col]]))
    }
    
    neighbornet_objective_minimum <- Inf
    # p_best_neighbornet <- NULL
    cycle_best <- NULL
    clus_df_gather_best <- NULL
    graphing_columns_best <- NULL
    objective_matrix_vector <- c()

    n <- length(cycle)

    # `data`/`wt` do not change across cycle-start iterations, so the base
    # alluvial dataframe is loop-invariant; compute it once instead of
    # recomputing (potentially an expensive group_by_all()+count()) on every
    # one of the n iterations below.
    if (is.null(wt) || length(wt) == 0 || !(wt %in% colnames(data))) {
        clus_df_gather_base <- get_alluvial_df(data, wt = wt)
    } else {
        clus_df_gather_base <- data
    }

    graphing_columns_tmp <- cols
    # A fixed axis keeps its given order whatever the cycle says, so reversing
    # the cycle is no longer equivalent to flipping the plot; try both directions.
    orientations <- if (length(fixed_orders) > 0) list(cycle, rev(cycle)) else list(cycle)
    for (orientation in seq_along(orientations)) {
    for (i in 0:(n - 1)) {
        if ((!is.null(cycle_start_positions)) && !((i + 1) %in% cycle_start_positions)) {
            next
        }
        if (verbose) message(sprintf("Starting iteration %s / %s%s", i + 1, n, if (orientation == 2) " (reversed cycle)" else ""))
        # if (i == 0) {
        #     if (verbose) message(sprintf("Starting iteration 1"))
        # } else if (i == 1) {
        #     if (verbose) message(sprintf("Starting subsequent iterations (should go much faster than iteration 1 if optimize_column_order is FALSE and/or optimize_column_order_per_cycle is FALSE)"))
        # }

        cycle_shifted <- rotate_left(orientations[[orientation]], i)
        graphs_list <- get_graph_groups(cycle_shifted)

        # remove prefix (column1_, etc)
        graphs_list_stripped <- lapply(graphs_list, function(x) {
            sub("^.*?~~", "", x)
        })
        graphs_list_stripped[names(fixed_orders)] <- fixed_orders

        clus_df_gather_neighbornet <- clus_df_gather_base

        graphing_columns_int <- c()
        for (j in seq_along(cols)) {
            col_name <- cols[j]
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
            if ((optimize_column_order_per_cycle) || (i == 0 && orientation == 1)) {
                # verbose_tmp <- verbose
                verbose_tmp <- if (i == 0 && orientation == 1) verbose else FALSE # only have the option for verbose on first iteration
                graphing_columns_tmp <- determine_column_order(clus_df_gather_neighbornet, cols = cols, wt = wt, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_metric = column_metric, column_method = column_method, verbose = verbose_tmp, weighted_metric = weighted_metric)
            }
        }
        
        clus_df_gather_neighbornet_tmp <- reorder_and_rename_columns(clus_df_gather_neighbornet, graphing_columns_tmp)
        
        neighbornet_objective <- compute_crossing_objective(
            clus_df_gather_neighbornet_tmp,
            cols = graphing_columns_tmp,
            wt = wt,
            weighted_metric = weighted_metric
        )$output_objective
        
        # print(neighbornet_objective)
        if (verbose) message(sprintf("objective for iteration %s = %s", i + 1, neighbornet_objective))
        if (neighbornet_objective < neighbornet_objective_minimum) {
            neighbornet_objective_minimum <- neighbornet_objective
            cycle_best <- cycle_shifted
            individual_graphs <- graphs_list_stripped
            # p_best_neighbornet <- p_neighbornet
            graphing_columns_best <- graphing_columns_tmp
            clus_df_gather_best <- clus_df_gather_neighbornet_tmp
        }
    }
    }
    
    # clus_df_gather_best <- reorder_and_rename_columns(clus_df_gather_best, graphing_columns_best)  # done earlier now
    
    # make factors
    for (j in seq_along(graphing_columns_best)) {
        int_col_name <- paste0("col", j, "_int")
        clus_df_gather_best[[int_col_name]] <- factor(clus_df_gather_best[[int_col_name]])
    }
    
    return(list(cycle = cycle_best, individual_graphs = individual_graphs, neighbornet_objective = neighbornet_objective_minimum, clus_df_gather = clus_df_gather_best, cols = graphing_columns_best))
}

increment_if_zeros <- function(clus_df_gather, column) {
    group_numeric <- as.numeric(as.character(clus_df_gather[[column]]))

    if (any(group_numeric == 0, na.rm = TRUE)) {
        group_numeric <- group_numeric + 1
        clus_df_gather[[column]] <- factor(group_numeric)
    }

    clus_df_gather
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

        # Base-R subsetting/assignment throughout this function instead of
        # dplyr pipes (mutate/slice/ungroup/select): this only ever runs over
        # a table capped at n_categories^2 rows, so the dplyr per-call
        # dispatch/NSE overhead (profiled dominating wall time here) is pure
        # fixed cost with no payoff at this scale. ungroup() is kept since a
        # caller could in principle pass a grouped tibble; nothing below
        # relies on grouping metadata.
        clus_df_gather <- dplyr::ungroup(clus_df_gather)
        subset_data <- clus_df_gather[(half_rows + 1):nrow(clus_df_gather), ]

        # Order matters: reordered_column_original_clusters_name captures the
        # pre-negation value, and best_cluster_agreement must read
        # reordered_column *after* it's negated on the next line (mirrors
        # dplyr::mutate()'s left-to-right sequential evaluation).
        subset_data[[reordered_column_original_clusters_name]] <- as.numeric(as.character(subset_data[[reordered_column]]))
        subset_data[[reordered_column]] <- -as.numeric(as.character(subset_data[[reordered_column]]))
        subset_data$y <- -as.numeric(as.character(subset_data$y))
        subset_data$best_cluster_agreement <- subset_data[[reordered_column]]

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
        
        subset_data[[reordered_column]] <- mapping[as.character(subset_data[[reordered_column]])]
        
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


add_int_columns <- function(data, cols, default_sorting = "alphabetical") {
    n <- 1
    for (col in cols) {
        col_int_name <- paste0("col", n, "_int")
        n <- n + 1
        
        # factorize input columns
        if (default_sorting == "alphabetical") {
            data[[col]] <- factor(as.character(data[[col]]), levels = sort(unique(as.character(data[[col]])), method = "radix"))
        } else if (default_sorting == "reverse_alphabetical") {
            data[[col]] <- factor(as.character(data[[col]]), levels = sort(unique(as.character(data[[col]])), method = "radix", decreasing = TRUE))
        } else if (default_sorting == "increasing") {
            data[[col]] <- factor(data[[col]], levels = names(sort(table(data[[col]]), decreasing = FALSE)))
        } else if (default_sorting == "decreasing") {
            data[[col]] <- factor(data[[col]], levels = names(sort(table(data[[col]]), decreasing = TRUE)))
        } else if (default_sorting == "random") {
            data[[col]] <- factor(data[[col]], levels = sample(unique(as.character(data[[col]]))))
        } else if (default_sorting == "fixed") {
            data[[col]] <- factor(data[[col]], levels = unique(as.character(data[[col]])))
        } else {
            stop(sprintf("default_sorting '%s' is not recognized. Please choose from 'alphabetical' (default), 'reverse_alphabetical', 'increasing', 'decreasing', or 'random'.", default_sorting))
        }
        
        if (!(col_int_name %in% colnames(data))) {
            # make columns integer for sorting
            data[[col_int_name]] <- as.integer(data[[col]])
        }
    }
    return(data)
}

# from ggforce
gather_set_data <- function (data, x, id_name = "id") 
{
    columns <- tidyselect::eval_select(rlang::enquo(x), data)
    data[[id_name]] <- seq_len(nrow(data))
    vctrs::vec_rbind(!!!lapply(names(columns), function(n) {
        data$x <- n
        data$y <- data[[n]]
        data
    }))
}

get_alluvial_df <- function(data, wt = "value", do_gather_set_data = FALSE) {
    if (is.null(wt) || length(wt) == 0) {
        wt <- "value"
    }
    # Convert numeric clustering columns to ordered factors
    data <- data |>
        dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
        dplyr::group_by_all() |>
        dplyr::count(name = wt)
    if (do_gather_set_data) {
        data <- gather_set_data(data, 1:2)
    }
    return(data)
}

# reorders int columns to match cols - eg if data has tissue, cluster, col1_int (for tissue), col2_int (for cluster) and cols = (cluster, tissue), then the output data will have cluster, tissue, col1_int (for cluster), col2_int (for tissue)
reorder_and_rename_columns <- function(data, cols) {
    # Find the order in data of the columns listed in cols
    original_graphing_columns <- intersect(colnames(data), cols)
    
    # Get original colX_int names based on original order
    original_int_cols <- paste0("col", seq_along(original_graphing_columns), "_int")
    
    # Target int column names based on desired new graphing order
    new_int_cols <- paste0("col", seq_along(cols), "_int")
    
    # Fix: old col name → new col name
    old_int_cols <- original_int_cols[match(cols, original_graphing_columns)]
    old_to_new_int_names <- setNames(new_int_cols, old_int_cols)
    
    # Rename data
    names(data)[names(data) %in% names(old_to_new_int_names)] <-
        old_to_new_int_names[names(data)[names(data) %in% names(old_to_new_int_names)]]
    
    
    graphing_and_int_columns <- union(cols, new_int_cols)
    
    # Final column order
    everything_else <- setdiff(names(data), graphing_and_int_columns)
    data <- data[, c(graphing_and_int_columns, everything_else)]
    
    return(data)
}

randomly_map_int_columns <- function(clus_df_gather, exclude = character(0)) {
    cols_to_shuffle <- setdiff(grep("^col\\d+_int$", names(clus_df_gather), value = TRUE), exclude)
    
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
#' Preprocess data (load in, add integer columns, reorder columns to match cols, and group as needed)
#'
#' @param data A data frame, tibble, or CSV file path. Must be in one of two formats:
#' (1) wt == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two cols).
#' (2) wt != NULL: Each row represents a combination of groupings, each column from \code{cols} represents a grouping, and the column \code{wt} represents the number of entities in that combination of groupings. Must contain at least three columns (two \code{cols}, one \code{wt}).
#' @param cols Character vector. Vector of column names from \code{data} to be used in graphing (i.e., alluvial plotting).
#' @param wt Optional character. Column name from \code{data} that contains the weights of each combination of groupings if \code{data} is in format (2) (see above).
#' @param default_sorting Character. Default column sorting in [prep_for_lodes()] if integer columns do not exist. Options are 'alphabetical' (default), 'reverse_alphabetical', 'increasing', 'decreasing', 'random'.
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param print_params Logical. If TRUE, will print function params.
#' @param do_gather_set_data Internal flag; not recommended to modify.
#' @param color_band_column Internal flag; not recommended to modify.
#' @param do_add_int_columns Internal flag; not recommended to modify.
#'
#' @return A data frame where each row represents a combination of groupings, each column from \code{cols} represents a grouping, and the column \code{wt} ('value' if \code{wt} == NULL) represents the number of entities in that combination of groupings. For each column in \code{cols}, there will be an additional column \code{col1_int}, \code{col2_int}, etc. where each column corresponds to a position mapping of groupings in the respective entry of \code{cols} - for example, \code{col1_int} corresponds to \code{cols[1]}, \code{col2_int} corresponds to \code{cols[2]}, etc.
#'
#' @examples
#' set.seed(235488)
#' data <- data.frame(
#'   method1 = sample(1:3, 100, TRUE),
#'   method2 = sample(4:6, 100, TRUE)
#' )
#' head(data)
#' lapply(data, unique)
#' 
#' # Example 1: data format 1
#' clus_df_gather <- prep_for_lodes(
#'   data,
#'   cols = c("method1", "method2")
#' )
#' print(clus_df_gather)
#' lapply(clus_df_gather[, 1:2], levels)
#'
#' # Example 2: data format 2
#' clus_df_gather <- data |>
#'     dplyr::mutate_if(
#'       is.numeric,
#'       function(x) factor(x, levels = as.character(sort(unique(x))))
#'     ) |>
#'     dplyr::group_by_all() |>
#'     dplyr::count(name = "value")
#' print(clus_df_gather)
#' lapply(clus_df_gather[, 1:2], unique)
#' clus_df_gather <- prep_for_lodes(
#'   clus_df_gather,
#'   cols = c("method1", "method2"),
#'   wt = "value"
#' )
#' print(clus_df_gather)
#' lapply(clus_df_gather[, 1:2], unique)
#'
#' @export
prep_for_lodes <- function(data, cols, wt = NULL, default_sorting = "alphabetical", verbose = FALSE, print_params = FALSE, do_gather_set_data = FALSE, color_band_column = NULL, do_add_int_columns = FALSE) {
    if (print_params) print_function_params()
    lowercase_args(c("default_sorting"))
    
    # remove all columns outside of cols, wt, and int columns
    cols_to_keep <- c(
        cols,
        wt,
        color_band_column,
        grep("^col[0-9]+_int$", names(data), value = TRUE)
    )
    cols_to_keep <- intersect(cols_to_keep, names(data))
    data <- data[, cols_to_keep, drop = FALSE]
    
    for (col in cols) {
        if (!(col %in% colnames(data))) {
            stop(sprintf("column '%s' is not a column in the dataframe.", col))
        }
        
        # convert to factor
        if (!is.factor(data[[col]])) {
            data[[col]] <- as.factor(data[[col]])
        }
        
        # fill in na with "Missing"
        data[[col]] <- replace(as.character(data[[col]]), is.na(data[[col]]), "Missing")
    }
    
    if (do_add_int_columns) {
        data <- add_int_columns(data, cols = cols, default_sorting = default_sorting)
    }
    
    # sort columns according to cols
    # data <- data |> dplyr::relocate(all_of(cols))  # put cols in front
    if (!all(intersect(colnames(data), cols) == cols)) {
        data <- reorder_and_rename_columns(data, cols)
    }
    
    wt_col <- get_col_name(data, wt)
    wt_col_in_data <- check_col(data, wt)
    if (!wt_col_in_data) {
        if (is.null(wt_col)) {
            wt_col <- "value"
        }
        clus_df_gather <- get_alluvial_df(data, wt = wt_col, do_gather_set_data = do_gather_set_data)
    } else {
        clus_df_gather <- data
    }
    
    clus_df_gather <- clus_df_gather |> dplyr::ungroup()
    
    # if ((is.character(output_df_path) && grepl("\\.rds$", output_df_path, ignore.case = TRUE))) {
    #     if (verbose) message(sprintf("Saving dataframe to=%s", output_df_path))
    #     saveRDS(clus_df_gather, output_df_path)
    # }
    
    return(clus_df_gather)
}



sort_neighbornet <- function(clus_df_gather, cols = NULL, wt = "value", optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE, matrix_initialization_value = 1e6, same_side_matrix_initialization_value = 1e6, weight_scalar = 5e5, matrix_initialization_value_column_order = 1e6, weight_scalar_column_order = 1, column_metric = "edge_crossing", method = "neighbornet", column_method = "tsp", cycle_start_positions = NULL, verbose = FALSE, weighted_metric = TRUE, fixed_column = character(0)) {
    if (verbose) message(sprintf("Running %s", method))
    fixed_orders <- lapply(setNames(fixed_column, fixed_column), function(col) {
        positions <- as.numeric(as.character(clus_df_gather[[int_col(match(col, cols))]]))
        unique(as.character(clus_df_gather[[col]])[order(positions)])
    })
    cycle <- run_neighbornet(clus_df_gather, cols = cols, wt = wt, matrix_initialization_value = matrix_initialization_value, same_side_matrix_initialization_value = same_side_matrix_initialization_value, weight_scalar = weight_scalar, method = method, verbose = verbose, fixed_orders = fixed_orders)
    if (verbose) message("Cycle: ", paste(cycle, collapse = ", "))
    if (verbose) message("Determining optimal cycle start")
    res <- determine_optimal_cycle_start(clus_df_gather, cycle, cols = cols, wt = wt, optimize_column_order = optimize_column_order, optimize_column_order_per_cycle = optimize_column_order_per_cycle, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_metric = column_metric, column_method = column_method, cycle_start_positions = cycle_start_positions, verbose = verbose, weighted_metric = weighted_metric, fixed_orders = fixed_orders)
    clus_df_gather_neighbornet <- res$clus_df_gather
    # graphing_columns_neighbornet <- res$cols
    if (verbose) message(sprintf("crossing edges objective = %s", res$neighbornet_objective))
    return(clus_df_gather_neighbornet)
}

int_col <- function(j) sprintf("col%d_int", as.integer(j))

deprecated_methods <- list(
    greedy_wolf = list(method = "greedy", one_sided = TRUE),
    greedy_wblf = list(method = "greedy", one_sided = FALSE),
    barycenter_one_sided = list(method = "barycenter", one_sided = TRUE),
    median_one_sided = list(method = "median", one_sided = TRUE)
)

# Resolves `fixed_column` (names or positions in `cols`) to unique column names.
resolve_fixed_columns <- function(fixed_column, cols) {
    if (is.null(fixed_column) || length(fixed_column) == 0) {
        return(character(0))
    }
    if (is.numeric(fixed_column)) {
        bad <- fixed_column[is.na(fixed_column) | fixed_column < 1 | fixed_column > length(cols) | fixed_column != round(fixed_column)]
        if (length(bad) > 0) {
            stop(sprintf("fixed_column position(s) %s are not positions in cols (1 to %d).", paste(bad, collapse = ", "), length(cols)))
        }
        fixed_column <- cols[fixed_column]
    } else if (is.character(fixed_column)) {
        bad <- setdiff(fixed_column, cols)
        if (length(bad) > 0) {
            stop(sprintf("fixed_column entries %s are not in cols.", paste0("'", bad, "'", collapse = ", ")))
        }
    } else {
        stop("fixed_column must be a character vector of names in cols or an integer vector of positions in cols.")
    }
    unique(fixed_column)
}

# Reassigns integer positions in `reordered_column` so that each of its blocks
# sits as close as possible to its heaviest neighbor in `stable_column`
# (greedy, O(n1 * n2)). The pair is first collapsed to one row per
# (stable, reordered) block combination, since with more than two axes a pair
# of adjacent blocks is split across many rows.
reorder_by_greedy_agreement <- function(clus_df_gather, stable_column, reordered_column, wt = "value") {
    stable <- as.numeric(as.character(clus_df_gather[[stable_column]]))
    reordered <- as.numeric(as.character(clus_df_gather[[reordered_column]]))
    # Dense ranks keep sort_clusters_by_agreement()'s processing order while
    # guaranteeing its factor assignments stay within the existing levels.
    reordered_dense <- match(reordered, sort(unique(reordered)))

    key <- paste(stable, reordered_dense)
    key_factor <- factor(key, levels = unique(key))
    first <- !duplicated(key)
    pair <- data.frame(
        col1_int = factor(stable[first]),
        col2_int = factor(reordered_dense[first]),
        value = as.numeric(rowsum(clus_df_gather[[wt]], key_factor, reorder = FALSE))
    )

    n_pairs <- nrow(pair)
    doubled <- sort_clusters_by_agreement(gather_set_data(pair, 1:2), stable_column = "col1_int", reordered_column = "col2_int")
    new_positions <- as.integer(as.character(doubled$col2_int[(n_pairs + 1):(2 * n_pairs)]))
    lookup <- setNames(new_positions, as.character(reordered_dense[first]))

    clus_df_gather[[reordered_column]] <- factor(unname(lookup[as.character(reordered_dense)]))
    clus_df_gather
}

weighted_median <- function(values, weights) {
    ord <- order(values)
    values <- values[ord]
    weights <- weights[ord]
    cum_w <- cumsum(weights)
    half <- sum(weights) / 2
    values[which(cum_w >= half)[1]]
}

# Reassigns integer positions in `reordered_column` by ranking each of its
# levels by the weighted mean (barycenter) or weighted median of the
# positions of its neighbors in `stable_column`, weighted by `wt`. This is
# the classic Sugiyama-style two-layer crossing-reduction heuristic: a cheap
# O(n log n) proxy for the crossing-count objective, rather than optimizing
# it directly like TSP or the greedy method (O(n1*n2) pairwise search).
reorder_by_neighbor_stat <- function(clus_df_gather, stable_column, reordered_column, wt = "value", stat = c("barycenter", "median")) {
    stat <- match.arg(stat)

    positions <- as.numeric(as.character(clus_df_gather[[stable_column]]))
    free_ids <- as.character(clus_df_gather[[reordered_column]])
    weights <- clus_df_gather[[wt]]

    split_idx <- split(seq_along(free_ids), free_ids)
    node_stats <- vapply(split_idx, function(idx) {
        if (stat == "barycenter") {
            sum(positions[idx] * weights[idx]) / sum(weights[idx])
        } else {
            weighted_median(positions[idx], weights[idx])
        }
    }, numeric(1))

    # Break ties deterministically by original position rather than relying
    # on split()'s (locale-dependent) name ordering.
    orig_ids <- as.numeric(names(node_stats))
    ord <- order(node_stats, orig_ids)
    new_rank <- integer(length(ord))
    new_rank[ord] <- seq_along(ord)
    new_positions <- setNames(new_rank, names(node_stats))

    clus_df_gather[[reordered_column]] <- factor(unname(new_positions[free_ids]))
    clus_df_gather
}


# The (reordered, stable) axis pairs of one forward-then-backward sweep. With no
# fixed axes every axis is reordered against its left neighbor going forward and
# its right neighbor going back. With fixed axes, order only propagates away
# from them: a free axis is reordered against its left (right) neighbor only if
# some fixed axis lies to its left (right). For two axes this reproduces the
# former greedy_wblf/barycenter (no fixed axis) and greedy_wolf/*_one_sided
# (one fixed axis) passes exactly.
sweep_passes <- function(n_cols, fixed_idx) {
    unconstrained <- length(fixed_idx) == 0
    free <- !(seq_len(n_cols) %in% fixed_idx)
    forward <- Filter(function(j) free[j] && (unconstrained || any(fixed_idx < j)), seq_len(n_cols)[-1])
    backward <- Filter(function(j) free[j] && (unconstrained || any(fixed_idx > j)), rev(seq_len(n_cols - 1)))
    c(
        lapply(forward, function(j) c(reordered = j, stable = j - 1)),
        lapply(backward, function(j) c(reordered = j, stable = j + 1))
    )
}

relabel_randomly <- function(x) {
    vals <- as.numeric(as.character(x))
    uniq <- sort(unique(vals))
    factor(sample(length(uniq))[match(vals, uniq)])
}

sort_by_sweep <- function(clus_df_gather, cols, wt = "value", method = c("greedy", "barycenter", "median"), fixed_column = NULL, column_method = "none", random_initializations = 1, weighted_metric = TRUE, matrix_initialization_value_column_order = 1e6, weight_scalar_column_order = 1, column_metric = "edge_crossing", verbose = FALSE) {
    method <- match.arg(method)
    reorder_pair <- if (method == "greedy") {
        function(df, stable_column, reordered_column) reorder_by_greedy_agreement(df, stable_column, reordered_column, wt = wt)
    } else {
        function(df, stable_column, reordered_column) reorder_by_neighbor_stat(df, stable_column, reordered_column, wt = wt, stat = method)
    }

    run_sweeps <- function(df, cols) {
        fixed_idx <- match(fixed_column, cols)
        free_idx <- setdiff(seq_along(cols), fixed_idx)
        passes <- sweep_passes(length(cols), fixed_idx)
        best <- NULL
        best_objective <- Inf
        for (i in seq_len(random_initializations)) {
            candidate <- df
            # The first initialization starts from the incoming order.
            if (i > 1) {
                for (j in free_idx) {
                    candidate[[int_col(j)]] <- relabel_randomly(candidate[[int_col(j)]])
                }
            }
            for (p in passes) {
                if (verbose) message(sprintf("Reordering %s (%s) by %s against %s (%s)", int_col(p[["reordered"]]), cols[p[["reordered"]]], method, int_col(p[["stable"]]), cols[p[["stable"]]]))
                candidate <- reorder_pair(candidate, int_col(p[["stable"]]), int_col(p[["reordered"]]))
            }
            if (random_initializations == 1) {
                return(candidate)
            }
            objective <- compute_crossing_objective(candidate, cols = cols, wt = wt, weighted_metric = weighted_metric)$output_objective
            if (verbose) message(sprintf("Initialization %d / %d: crossing edges objective = %s", i, random_initializations, objective))
            if (objective < best_objective) {
                best_objective <- objective
                best <- candidate
            }
        }
        best
    }

    if (length(fixed_column) == length(cols)) {
        if (verbose) message("Every column is fixed; nothing to sort")
        return(clus_df_gather)
    }

    clus_df_gather <- run_sweeps(clus_df_gather, cols)

    # Which strata cross depends on which axes are adjacent, so after choosing
    # an axis order the strata are swept again against their new neighbors.
    if (length(cols) > 2 && column_method != "none") {
        cols_ordered <- determine_column_order(clus_df_gather, cols = cols, wt = wt, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_metric = column_metric, column_method = column_method, verbose = verbose, weighted_metric = weighted_metric)
        if (!identical(cols_ordered, cols)) {
            if (verbose) message(sprintf("Column order: %s", paste(cols_ordered, collapse = ", ")))
            clus_df_gather <- run_sweeps(reorder_and_rename_columns(clus_df_gather, cols_ordered), cols_ordered)
        }
    }

    clus_df_gather
}

#' Control Options for `sort_to_uncross()`
#'
#' Creates a list of control parameters that modify the behavior of
#' [sort_to_uncross()]. These options allow tuning algorithmic behavior without
#' cluttering the main function arguments.
#'
#' @param optimize_column_order_per_cycle Logical. If TRUE, will optimize the order of \code{cols} to minimize edge overlap upon each cycle. If FALSE, will optimize the order of \code{cols} to minimize edge overlap on the beginning cycle only. Only applies when \code{method \%in\% c('neighbornet', 'tsp')} and \code{length(cols) > 2}.
#' @param weight_scalar Advanced. Positive number \eqn{c} by which \eqn{-\log(\text{edge weight})} is multiplied in the block distance matrix. Because the NeighborNet cycle is invariant under rescaling of the whole matrix (up to floating-point ties), \code{weight_scalar} only sets the absolute scale of the matrix; the ratios that actually shape the cycle are \code{alpha} and \code{beta} in [sort_to_uncross()]. Only applies when \code{method \%in\% c('neighbornet', 'tsp')}.
#' @param matrix_initialization_value Advanced. Positive number \eqn{d_{\max}}: the distance placed between two blocks in different axes that share no observations. \code{NULL} (default) derives it as \code{alpha * weight_scalar}; supplying a value overrides \code{alpha}. Only applies when \code{method \%in\% c('neighbornet', 'tsp')}.
#' @param same_side_matrix_initialization_value Advanced. Positive number \eqn{d_{\text{same}}}: the distance placed between two distinct blocks of the same axis. \code{NULL} (default) derives it as \code{beta * weight_scalar}; supplying a value overrides \code{beta}. Only applies when \code{method \%in\% c('neighbornet', 'tsp')}.
#' @param matrix_initialization_value_column_order Positive integer. Initialized value in distance matrix for optimizing column order. Only applies when \code{column_method != 'none'}.
#' @param weight_scalar_column_order Positive integer. Scalar with which to loss function after taking their log1p in the distance matrix for optimizing column order. Only applies when \code{column_method != 'none'}.
#' @param column_metric Character. Metric to use for determining column order. Options are "edge_crossing" (default) or "ARI". Only applies when \code{column_method != 'none'}.
#' @param weighted_metric Logical. Determines if the objective is total number of edge crossings (weighted_metric=FALSE) or sum of product of overlapping edge weights (weighted_metric=TRUE).
#' @param cycle_start_positions Set. Cycle start positions to consider. Anything outside this set will be skipped. Only applies when \code{method \%in\% c('neighbornet', 'tsp')}.
#' @param random_initializations Integer. Number of initializations of the stratum positions of the non-fixed axes: the first starts from the incoming order and each further one from a random order, and the one with the lowest crossing objective is kept. Only applies when \code{method \%in\% c('greedy', 'barycenter', 'median')}.
#' @param preprocess_data Logical. If TRUE, will preprocess the data with the [prep_for_lodes()] function.
#' @param default_sorting Character. Default column sorting in [prep_for_lodes()] if integer columns do not exist. Options are 'alphabetical' (default), 'reverse_alphabetical', 'increasing', 'decreasing', 'random'.
#' @param print_params Logical. If TRUE, will print function params.
#' @param do_compute_alluvial_statistics Internal flag; not recommended to modify.
#' 
#' @return A named list of control parameters, to be passed into [sort_to_uncross()]
#'   via the `options` argument.
#'
#' @examples
#' data <- data.frame(
#'   method1 = LETTERS[sample(1:3, 100, TRUE)],
#'   method2 = LETTERS[27 - sample(1:3, 100, TRUE)]
#' )
#' opts <- sort_to_uncross_options(
#'   default_sorting = "reverse_alphabetical",
#'   weighted_metric = FALSE
#' )
#' sort_to_uncross(data = data, cols = c('method1', 'method2'), options = opts)
#'
#' @export
sort_to_uncross_options <- function(
        optimize_column_order_per_cycle = FALSE,
        weight_scalar = 5e5,
        matrix_initialization_value = NULL,
        same_side_matrix_initialization_value = NULL,
        matrix_initialization_value_column_order = 1e6,
        weight_scalar_column_order = 1,
        column_metric = c("edge_crossing", "ari"),
        weighted_metric = TRUE,
        cycle_start_positions = NULL,
        random_initializations = 1,
        preprocess_data = TRUE,
        default_sorting = c("alphabetical", "reverse_alphabetical", "increasing", "decreasing", "random", "fixed"),
        print_params = FALSE,
        do_compute_alluvial_statistics = FALSE
) {
    # enforce match.arg for any fixed-choice settings
    column_metric <- match.arg(column_metric)
    default_sorting <- match.arg(default_sorting)
    
    list(
        optimize_column_order_per_cycle = optimize_column_order_per_cycle,
        weight_scalar = weight_scalar,
        matrix_initialization_value = matrix_initialization_value,
        same_side_matrix_initialization_value = same_side_matrix_initialization_value,
        matrix_initialization_value_column_order = matrix_initialization_value_column_order,
        weight_scalar_column_order = weight_scalar_column_order,
        column_metric = column_metric,
        weighted_metric = weighted_metric,
        cycle_start_positions = cycle_start_positions,
        random_initializations = random_initializations,
        preprocess_data = preprocess_data,
        default_sorting = default_sorting,
        print_params = print_params,
        do_compute_alluvial_statistics = do_compute_alluvial_statistics
    )
}

sort_to_uncross_internal <- function(data, cols, wt = NULL, method = c("neighbornet", "tsp", "greedy", "barycenter", "median", "none", "random"), column_method = c("tsp", "neighbornet", 'none', 'random'), alpha = 2, beta = alpha, fixed_column = NULL, output_df_path = NULL, verbose = FALSE, options = NULL) {
    default_opt <- sort_to_uncross_options()
    if (!is.null(options)) {
        if (!is.list(options)) stop("`options` must be a list.")
        for (nm in names(default_opt)) {
            if (!nm %in% names(options)) {
                val <- ifelse(!is.null(default_opt[[nm]]), default_opt[[nm]], "NULLTMP")
                options[[nm]] <- val
            }
        }
    } else {
        options <- default_opt
        for (nm in names(options)) {
            if (is.null(options[[nm]])) {
                options[[nm]] <- "NULLTMP"
            }
        }
    }
    for (nm in names(options)) {
        if (is.null(options[[nm]]) || options[[nm]] == "NULLTMP") {
            val <- NULL
        } else {
            val <- options[[nm]]
        }
        assign(nm, val, envir = environment())
    }
    
    if (print_params) print_function_params()
    # lowercase_args(c("method", "column_metric", "column_method", "default_sorting"))

    #* Type Checking Start
    if (length(method) == 1 && method %in% names(deprecated_methods)) {
        replacement <- deprecated_methods[[method]]
        fixed_hint <- if (replacement$one_sided) sprintf(", fixed_column = %s", if (is.null(fixed_column)) "cols[1]" else "<fixed_column>") else ""
        warning(sprintf("method = '%s' is deprecated; use method = '%s'%s instead.", method, replacement$method, fixed_hint), call. = FALSE)
        if (replacement$one_sided && is.null(fixed_column)) {
            fixed_column <- cols[1]
        } else if (!replacement$one_sided) {
            fixed_column <- NULL
        }
        method <- replacement$method
    }
    method <- match.arg(method)
    column_method <- match.arg(column_method)
    fixed_column <- resolve_fixed_columns(fixed_column, cols)

    if (method == "tsp" || method == "neighbornet") {
        for (col in cols) {
            if (grepl("~~", col)) {
                stop(sprintf("No entry of cols can contain '~~' when method == tsp. Issue with column '%s'.", col))
            }
        }
    }
    
    if (!(method %in% c("greedy", "barycenter", "median")) && (random_initializations > 1)) {
        if (verbose) message(sprintf("random_initializations > 1 but sorting algorithm is %s. Setting random_initializations to 1.", method))
        random_initializations <- 1
    }

    if (ncol(data) < 2) {
        stop("Dataframe must have at least 2 columns.")
    }

    if (length(cols) < 2) {
        stop("cols must have at least 2 entries.")
    }
    
    if (any(!cols %in% colnames(data))) {
        stop("Some cols are not present in the dataframe.")
    }
    
    optimize_column_order <- (column_method != "none")

    # The block distance matrix (Eq. distmat of the paper) has three constants:
    # c = weight_scalar, d_max = matrix_initialization_value and
    # d_same = same_side_matrix_initialization_value. The cycle depends on them
    # only through the ratios d_max / c and d_same / c, so the user-facing
    # parameters are those ratios (alpha, beta); the constants themselves are
    # advanced overrides in sort_to_uncross_options().
    if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) stop("`alpha` must be a single positive number.")
    if (!is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) stop("`beta` must be a single positive number.")
    if (!is.numeric(weight_scalar) || length(weight_scalar) != 1 || !is.finite(weight_scalar) || weight_scalar <= 0) stop("`weight_scalar` must be a single positive number.")
    if (is.null(matrix_initialization_value)) {
        matrix_initialization_value <- alpha * weight_scalar
    } else if (verbose) {
        message("matrix_initialization_value supplied; ignoring `alpha`")
    }
    if (is.null(same_side_matrix_initialization_value)) {
        same_side_matrix_initialization_value <- beta * weight_scalar
    } else if (verbose) {
        message("same_side_matrix_initialization_value supplied; ignoring `beta`")
    }
    #* Type Checking End
    
    # Preprocess (i.e., add int columns and do the grouping)
    if (preprocess_data) {
        if (verbose) message("Preprocessing data before sorting")
        clus_df_gather <- prep_for_lodes(data = data, cols = cols, wt = wt, default_sorting = default_sorting, do_gather_set_data = FALSE, do_add_int_columns = TRUE)
        if (is.null(wt) || length(wt) == 0) {
            wt <- "value" # is set during prep_for_lodes
        }
    } else {
        clus_df_gather <- data
    }
    
    if (verbose && do_compute_alluvial_statistics) compute_alluvial_statistics(clus_df_gather = clus_df_gather, cols = cols, wt = wt)
    if (method == "neighbornet" || method == "tsp") {
        # O(n^3) complexity, where n is the sum of blocks across all layers
        clus_df_gather_sorted <- sort_neighbornet(clus_df_gather = clus_df_gather, cols = cols, wt = wt, optimize_column_order = optimize_column_order, optimize_column_order_per_cycle = optimize_column_order_per_cycle, matrix_initialization_value = matrix_initialization_value, same_side_matrix_initialization_value = same_side_matrix_initialization_value, weight_scalar = weight_scalar, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_metric = column_metric, method = method, column_method = column_method, cycle_start_positions = cycle_start_positions, verbose = verbose, weighted_metric = weighted_metric, fixed_column = fixed_column)
    } else if (method %in% c("greedy", "barycenter", "median")) {
        # Per adjacent-axis pass: O(n_i * n_j) for greedy, O(a log a) for barycenter/median (a = number of alluvia)
        clus_df_gather_sorted <- sort_by_sweep(clus_df_gather = clus_df_gather, cols = cols, wt = wt, method = method, fixed_column = fixed_column, column_method = column_method, random_initializations = random_initializations, weighted_metric = weighted_metric, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_metric = column_metric, verbose = verbose)
    } else if (method == "random") {
        clus_df_gather_sorted <- randomly_map_int_columns(clus_df_gather, exclude = int_col(match(fixed_column, cols)))
        #!!! check this
        if (optimize_column_order) {
            graphing_columns_tmp <- determine_column_order(clus_df_gather_sorted, cols = cols, wt = wt, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_metric = column_metric, column_method = column_method, verbose = verbose, weighted_metric = weighted_metric)
            clus_df_gather_tmp <- reorder_and_rename_columns(clus_df_gather_sorted, graphing_columns_tmp)
            # make factors
            for (j in seq_along(graphing_columns_tmp)) {
                int_col_name <- paste0("col", j, "_int")
                clus_df_gather_tmp[[int_col_name]] <- factor(clus_df_gather_tmp[[int_col_name]])
            }
            clus_df_gather_sorted <- clus_df_gather_tmp
        }
        #!!! check this
    } else if (method == "none") {
        clus_df_gather_sorted <- clus_df_gather
        #!!! check this
        if (optimize_column_order) {
            graphing_columns_tmp <- determine_column_order(clus_df_gather_sorted, cols = cols, wt = wt, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_metric = column_metric, column_method = column_method, verbose = verbose, weighted_metric = weighted_metric)
            clus_df_gather_tmp <- reorder_and_rename_columns(clus_df_gather_sorted, graphing_columns_tmp)
            # make factors
            for (j in seq_along(graphing_columns_tmp)) {
                int_col_name <- paste0("col", j, "_int")
                clus_df_gather_tmp[[int_col_name]] <- factor(clus_df_gather_tmp[[int_col_name]])
            }
            clus_df_gather_sorted <- clus_df_gather_tmp
        }
        #!!! check this
    } else {
        stop(sprintf("Invalid method: '%s'. Must be one of: %s", method, paste(valid_algorithms, collapse = ", ")))
    }
    
    # print objective - don't do for neighbornet because I did it right before
    if ((verbose) && (method != "neighbornet")) {
        message("Determining crossing edges objective (to disable, use verbose==FALSE)")
        objective <- compute_crossing_objective(clus_df_gather_sorted, cols = cols, wt = wt, weighted_metric = weighted_metric)$output_objective
        message(sprintf("crossing edges objective = %s", objective))
    }
    
    if (verbose) message("Complete with sorting")
        
    # reorder cols to match any changed order in clus_df_gather_sorted
    graphing_columns_sorted <- cols[order(match(cols, names(clus_df_gather_sorted)))]
    clus_df_gather <- generalized_reorder(clus_df_gather=clus_df_gather, clus_df_gather_sorted=clus_df_gather_sorted, cols=graphing_columns_sorted)
    clus_df_gather <- clus_df_gather |> dplyr::ungroup()
    
    # Save if desired
    if ((is.character(output_df_path) && grepl("\\.rds$", output_df_path, ignore.case = TRUE))) {
        if (verbose) message(sprintf("Saving sorted dataframe to=%s", output_df_path))
        saveRDS(clus_df_gather, output_df_path)
    }
    
    return(clus_df_gather)
}


#' Sorts a dataframe to minimize crossings in a parallel sets / alluvial plot.
#'
#' Sorts a dataframe with the algorithm specified by \code{method}.
#'
#' @param data A data frame/tibble. Must be in one of two formats:
#' (1) wt == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two cols).
#' (2) wt != NULL: Each row represents a combination of groupings, each column from \code{cols} represents a grouping, and the column \code{wt} represents the number of entities in that combination of groupings. Must contain at least three columns (two \code{cols}, one \code{wt}).
#' @param cols Character vector. Vector of column names from \code{data} to be used in graphing (i.e., alluvial plotting).
#' @param wt Optional character. Column name from \code{data} that contains the weights of each combination of groupings if \code{data} is in format (2) (see above).
#' @param method Character. Algorithm with which to sort the values in the dataframe. Can choose from: 'neighbornet' (default), 'tsp', 'greedy', 'barycenter', 'median', 'random', 'none'. All of them accept any number of \code{cols} and honor \code{fixed_column}. 'neighbornet' builds the block distance matrix and orders the blocks with the NeighborNet algorithm (the W_POMP method described in the paper). 'tsp' is identical except the block ordering is produced by the Traveling Salesman Problem solver from the TSP package rather than NeighborNet. 'greedy', 'barycenter' and 'median' are layer-by-layer sweeps: one pass forward across the axes and one back, each pass reordering an axis against a neighboring axis. 'greedy' places each block as close as possible to its heaviest neighbor in the neighboring axis (O(n_i * n_j) per pass); 'barycenter' and 'median' are the classic Sugiyama-style heuristics that place each block at the weighted mean or weighted median position of its neighbors (O(a log a) per pass, for a alluvia), typically at the cost of noisier (less minimized) crossing counts. 'random' randomly maps blocks. 'none' keeps the mappings as-is when passed into the function. The former names 'greedy_wolf', 'greedy_wblf', 'barycenter_one_sided' and 'median_one_sided' are deprecated: 'greedy_wblf' is 'greedy', and the others are 'greedy', 'barycenter' or 'median' with \code{fixed_column} (defaulting to \code{cols[1]}).
#' @param column_method Character. Algorithm to use for determining column order. Options are 'tsp' (default), 'neighbornet', 'random', and 'none'. Only applies when \code{cols} has more than two entries. With 'greedy', 'barycenter' and 'median', the strata are swept once in the given axis order, the axis order is chosen from that result, and the strata are swept again in the new order.
#' @param alpha Positive number (default 2). Ratio \eqn{d_{\max} / c} between the distance assigned to two blocks in different axes that share no observations and the scale \eqn{c} of the \eqn{-\log(\text{edge weight})} distances between blocks that do. Larger values hold unconnected blocks further apart relative to the spread induced by differences in edge weight. Together with \code{beta} this is the only tuning knob of the block distance matrix: the NeighborNet (and TSP) cycle is unchanged when the whole matrix is rescaled, so the absolute values \eqn{c}, \eqn{d_{\max}} and \eqn{d_{\text{same}}} matter only through the ratios \code{alpha} and \code{beta} (see \code{weight_scalar}, \code{matrix_initialization_value} and \code{same_side_matrix_initialization_value} in [sort_to_uncross_options()] for the underlying constants). Only applies when \code{method \%in\% c('neighbornet', 'tsp')}.
#' @param beta Positive number (default \code{alpha}). Ratio \eqn{d_{\text{same}} / c} between the distance assigned to two distinct blocks of the same axis and the scale \eqn{c} of the edge-weight distances. Only applies when \code{method \%in\% c('neighbornet', 'tsp')}.
#' @param fixed_column Optional character or integer vector. Names, or positions in \code{cols}, of the axes whose stratum order is kept as it comes in (the order \code{method = 'none'} would return, set by \code{default_sorting} in [sort_to_uncross_options()] unless the data already carry \code{col*_int} columns); only the other axes are sorted. With 'greedy', 'barycenter' and 'median', order propagates outward from the fixed axes. With 'neighbornet' and 'tsp', which order a distance matrix rather than axes, the fixed order is encouraged in the distance matrix and then imposed exactly on the resulting cycle, whose start and direction are chosen to minimize crossings with the fixed axes in place. Because the cycle is only nudged toward the fixed order, 'greedy', 'barycenter' and 'median' usually reach fewer crossings than 'neighbornet' and 'tsp' when axes are fixed. \code{NULL} (default) sorts every axis. Does not affect the order of the axes themselves (see \code{column_method}).
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param options Additional arguments. See [sort_to_uncross_options()].
#'
#' @return
#' A data frame where each row represents an alluvium and each column represents an axis. Each column of \code{cols} represents an axis, stored as a factor ordered in ascending order of strata. There is an additional column \code{wt} ('value' if NULL) that represents the size of the alluvium in that row. The order of columns represents the recommended order of axes.
#'
#' @examples
#' # Example 1: data format 1 (uncounted)
#' set.seed(429144)
#' data <- data.frame(
#'   method1 = factor(LETTERS[sample(1:3, 100, TRUE)]),
#'   method2 = factor(LETTERS[27 - sample(1:3, 100, TRUE)])
#' )
#' head(data)
#' lapply(data, levels)
#' clus_df_gather <- sort_to_uncross(
#'   data,
#'   cols = c("method1", "method2"),
#'   method = "tsp",
#'   column_method = "tsp"
#' )
#' print(clus_df_gather)
#' lapply(clus_df_gather[, 1:2], levels)
#'
#' # Example 2: data format 2 (counted)
#' set.seed(806949)
#' data <- data.frame(
#'   method1 = factor(LETTERS[sample(1:3, 100, TRUE)]),
#'   method2 = factor(LETTERS[27 - sample(1:3, 100, TRUE)])
#' )
#' clus_df_gather <- data |>
#'   dplyr::mutate_if(
#'     is.numeric,
#'     function(x) factor(x, levels = as.character(sort(unique(x))))
#'   ) |>
#'   dplyr::group_by_all() |>
#'   dplyr::count(name = "value")
#' print(clus_df_gather)
#' lapply(clus_df_gather[, 1:2], levels)
#' clus_df_gather <- sort_to_uncross(
#'   clus_df_gather,
#'   cols = c("method1", "method2"),
#'   wt = "value",
#'   method = "tsp",
#'   column_method = "tsp"
#' )
#' print(clus_df_gather)
#' lapply(clus_df_gather[, 1:2], levels)
#'
#' @export
sort_to_uncross <- function(data, cols, wt = NULL, method = c("neighbornet", "tsp", "greedy", "barycenter", "median", "none", "random"), column_method = c("tsp", "neighbornet", 'none', 'random'), alpha = 2, beta = alpha, fixed_column = NULL, verbose = FALSE, options = NULL) {
    default_opt <- sort_to_uncross_options()
    if (!is.null(options)) {
        if (!is.list(options)) stop("`options` must be a list.")
        for (nm in names(default_opt)) {
            if (!nm %in% names(options)) {
                val <- ifelse(!is.null(default_opt[[nm]]), default_opt[[nm]], "NULLTMP")
                options[[nm]] <- val
            }
        }
    } else {
        options <- default_opt
        for (nm in names(options)) {
            if (is.null(options[[nm]])) {
                options[[nm]] <- "NULLTMP"
            }
        }
    }
    for (nm in names(options)) {
        if (is.null(options[[nm]]) || options[[nm]] == "NULLTMP") {
            val <- NULL
        } else {
            val <- options[[nm]]
        }
        assign(nm, val, envir = environment())
    }

    if (missing(wt)) {
        data <- prep_for_lodes(data = data, cols = cols, default_sorting = default_sorting, do_gather_set_data = FALSE, do_add_int_columns = TRUE)
        wt <- "value" # is set during prep_for_lodes
    }
    
    # `cols` and `wt` are documented as character vectors; select by value with
    # all_of() so a variable holding the name works (ensym() resolved the
    # variable name itself, not its value) and no external-vector warning fires.
    cols_pos <- tidyselect::eval_select(tidyselect::all_of(cols), data = data)
    wt_pos <- tidyselect::eval_select(tidyselect::all_of(wt), data = data)
    # Pass the data through unchanged (rather than subsetting to cols + wt) so
    # that pre-computed `col*_int` columns, a `color_band_column`, and any other
    # metadata survive when `preprocess_data = FALSE`.
    sort_to_uncross_internal(data = data, cols = names(cols_pos), wt = names(wt_pos), method = method, column_method = column_method, alpha = alpha, beta = beta, fixed_column = fixed_column, verbose = verbose, options = options)
}
