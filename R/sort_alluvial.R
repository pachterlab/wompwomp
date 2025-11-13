#' wompwomp: Sorting alluvia
#'
#' Sorting alluvia
#' @docType package
#' @name wompwomp
#'
#' @importFrom dplyr mutate select group_by summarise desc ungroup slice n pull bind_rows across all_of
#' @importFrom purrr map
#' @importFrom igraph V cluster_louvain cluster_leiden E
#' @importFrom tibble is_tibble
#' @importFrom utils read.csv write.csv combn
#' @importFrom stats setNames
#' @importFrom rlang sym .data
#' @importFrom magrittr %>%
#' @importFrom data.table :=

utils::globalVariables(c(
    ".data", ":=", "%>%", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "stratum", "total", "cum_y", "best_cluster_agreement", "neighbor_net", "alluvium", "pos", "count", "group1", "group2", "value", "group1_size", "group2_size", "weight", "parent", "group_name"
))

# devtools::document() needs package = "wompwomp"; but R CMD build wompwomp needs it not there stored globally - so I make this function
get_neighbornet_script_path <- function() {
    neighbornet_script_path <- system.file("scripts", "run_neighbornet.py", package = "wompwomp")
    if (neighbornet_script_path == "") {
        # Fallback to development location
        neighbornet_script_path <- file.path(here::here("inst", "scripts", "run_neighbornet.py"))
        if (!grepl("wompwomp", neighbornet_script_path)) {
            # if here::here isn't putting it inside wompwomp
            neighbornet_script_path <- file.path("inst", "scripts", "run_neighbornet.py")
        }
    }
    neighbornet_script_path <- normalizePath(neighbornet_script_path, mustWork = TRUE)
    stopifnot(file.exists(neighbornet_script_path))
    return (neighbornet_script_path)
}

#reticulate::source_python(neighbornet_script_path)  # Error: Unable to access object (object is from previous session and is now invalid)

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

determine_column_order <- function(clus_df_gather_neighbornet, graphing_columns, column_weights = "value", matrix_initialization_value_column_order = 1e6, weight_scalar_column_order = 1, column_sorting_metric = "edge_crossing", column_sorting_algorithm = "tsp", verbose = FALSE, weighted = TRUE) {
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
                weighted = weighted
            )$output_objective
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
        check_python_setup_with_necessary_packages(necessary_packages_for_this_step = c("splitspy", "numpy"), additional_message = "do not set column_sorting_algorithm to 'neighbornet'")
        reticulate::source_python(get_neighbornet_script_path())
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

run_neighbornet <- function(df, graphing_columns = NULL, column1 = NULL, column2 = NULL, column_weights = "value", matrix_initialization_value = 1e6, same_side_matrix_initialization_value = 1e6, weight_scalar = 5e5, sorting_algorithm = "tsp", verbose = FALSE) {
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
        check_python_setup_with_necessary_packages(necessary_packages_for_this_step = c("splitspy", "numpy"), additional_message = "do not set sorting_algorithm to 'neighbornet'")
        reticulate::source_python(get_neighbornet_script_path())
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



determine_optimal_cycle_start <- function(df, cycle, graphing_columns = NULL, column1 = NULL, column2 = NULL, column_weights = "value", optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE, matrix_initialization_value_column_order = 1e6, weight_scalar_column_order = 1, column_sorting_metric = "edge_crossing", column_sorting_algorithm = "tsp", cycle_start_positions = NULL, verbose = FALSE, weighted = TRUE) {
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
                graphing_columns_tmp <- determine_column_order(clus_df_gather_neighbornet, graphing_columns = graphing_columns, column_weights = column_weights, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_sorting_metric = column_sorting_metric, column_sorting_algorithm = column_sorting_algorithm, verbose = verbose_tmp, weighted = weighted)
            }
        }
        
        clus_df_gather_neighbornet_tmp <- reorder_and_rename_columns(clus_df_gather_neighbornet, graphing_columns_tmp)
        
        neighbornet_objective <- determine_crossing_edges(
            clus_df_gather_neighbornet_tmp,
            graphing_columns = graphing_columns_tmp,
            column_weights = column_weights,
            weighted = weighted
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
    
    
    graphing_and_int_columns <- union(graphing_columns, new_int_cols)
    
    # Final column order
    everything_else <- setdiff(names(df), graphing_and_int_columns)
    df <- df[, c(graphing_and_int_columns, everything_else)]
    
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
#' @param default_sorting Character. Default column sorting in data_preprocess if integer columns do not exist. Options are 'alphabetical' (default), 'reverse_alphabetical', 'increasing', 'decreasing', 'random'.
#' @param output_df_path Optional character. Output path for the output data frame, in CSV format. If \code{NULL}, then will not be saved.
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param print_params Logical. If TRUE, will print function params.
#' @param load_df Internal flag; not recommended to modify.
#' @param do_gather_set_data Internal flag; not recommended to modify.
#' @param color_band_column Internal flag; not recommended to modify.
#' @param do_add_int_columns Internal flag; not recommended to modify.
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
data_preprocess <- function(df, graphing_columns, column_weights = NULL, default_sorting = "alphabetical", output_df_path = NULL, verbose = FALSE, print_params = FALSE, load_df = TRUE, do_gather_set_data = FALSE, color_band_column = NULL, do_add_int_columns = TRUE) {
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
    for (col in graphing_columns) {
        if (col %in% colnames(df)) {
            if (is.factor(df[[col]]) || is.character(df[[col]])) {
                df[[col]] <- replace(as.character(df[[col]]), is.na(df[[col]]), "Missing")
            }
        }
    }
    
    if (do_add_int_columns) {
        df <- add_int_columns(df, graphing_columns = graphing_columns, default_sorting = default_sorting)
    }
    
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



sort_neighbornet <- function(clus_df_gather, graphing_columns = NULL, column_weights = "value", optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE, matrix_initialization_value = 1e6, same_side_matrix_initialization_value = 1e6, weight_scalar = 5e5, matrix_initialization_value_column_order = 1e6, weight_scalar_column_order = 1, column_sorting_metric = "edge_crossing", sorting_algorithm = "tsp", column_sorting_algorithm = "tsp", cycle_start_positions = NULL, verbose = FALSE, weighted = TRUE) {
    if (verbose) message(sprintf("Running %s", sorting_algorithm))
    cycle <- run_neighbornet(clus_df_gather, graphing_columns = graphing_columns, column_weights = column_weights, matrix_initialization_value = matrix_initialization_value, same_side_matrix_initialization_value = same_side_matrix_initialization_value, weight_scalar = weight_scalar, sorting_algorithm = sorting_algorithm, verbose = verbose)
    if (verbose) message("Cycle: ", paste(cycle, collapse = ", "))
    if (verbose) message("Determining optimal cycle start")
    res <- determine_optimal_cycle_start(clus_df_gather, cycle, graphing_columns = graphing_columns, column_weights = column_weights, optimize_column_order = optimize_column_order, optimize_column_order_per_cycle = optimize_column_order_per_cycle, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_sorting_metric = column_sorting_metric, column_sorting_algorithm = column_sorting_algorithm, cycle_start_positions = cycle_start_positions, verbose = verbose, weighted = weighted)
    clus_df_gather_neighbornet <- res$clus_df_gather
    # graphing_columns_neighbornet <- res$graphing_columns
    if (verbose) message(sprintf("crossing edges objective = %s", res$neighbornet_objective))
    return(clus_df_gather_neighbornet)
}

sort_greedy_wolf <- function(clus_df_gather, graphing_columns = NULL, column1 = NULL, column2 = NULL, fixed_column = NULL, weighted = TRUE, random_initializations = 1, column_weights = "value", sorting_algorithm = "greedy_wblf", verbose = FALSE) {
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
            crossing_edges_objective <- determine_crossing_edges(clus_df_gather_tmp, graphing_columns = graphing_columns, column_weights = column_weights, weighted = weighted)$output_objective
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
#' @param sorting_algorithm Character. Algorithm with which to sort the values in the dataframe. Can choose from: 'tsp', 'greedy_wolf', 'greedy_wblf', 'none'. 'tsp' performs Traveling Salesman Problem solver from the TSP package. greedy_wolf' implements a custom greedy algorithm where one layer is fixed, and the other layer is sorted such that each node is positioned as close to its largest parent from the fixed side as possible in a greedy fashion. 'greedy_wblf' implements the 'greedy_wolf' algorithm described previously twice, treating each column as fixed in one iteration and free in the other iteration. 'greedy_wolf' and 'greedy_wblf' are only valid when \code{graphing_columns} has exactly two entries. 'random' randomly maps blocks. 'none' keeps the mappings as-is when passed into the function.
#' @param optimize_column_order Logical. If TRUE, will optimize the order of \code{graphing_columns} to minimize edge overlap. Only applies when \code{sorting_algorithm == 'tsp'} and \code{length(graphing_columns) > 2}.
#' @param optimize_column_order_per_cycle Logical. If TRUE, will optimize the order of \code{graphing_columns} to minimize edge overlap upon each cycle. If FALSE, will optimize the order of \code{graphing_columns} to minimize edge overlap on the beginning cycle only. Only applies when \code{sorting_algorithm == 'tsp'} and \code{length(graphing_columns) > 2}.
#' @param matrix_initialization_value Positive integer. Initialized value in distance matrix for nodes in different layers without a shared edge/path. Only applies when \code{sorting_algorithm == 'tsp'}.
#' @param same_side_matrix_initialization_value Positive integer. Initialized value in distance matrix for nodes in the same layer. Only applies when \code{sorting_algorithm == 'tsp'}.
#' @param weight_scalar Positive integer. Scalar with which to multiply edge weights after taking their -log in the distance matrix for nodes with a nonzero edge. Only applies when \code{sorting_algorithm == 'tsp'}.
#' @param matrix_initialization_value_column_order Positive integer. Initialized value in distance matrix for optimizing column order. Only applies when \code{sorting_algorithm == 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param weight_scalar_column_order Positive integer. Scalar with which to loss function after taking their log1p in the distance matrix for optimizing column order. Only applies when \code{sorting_algorithm == 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param column_sorting_metric Character. Metric to use for determining column order. Options are "edge_crossing" (default) or "ARI". Only applies when \code{sorting_algorithm == 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param column_sorting_algorithm Character. Algorithm to use for determining column order. Options are "tsp" (default). Only applies when \code{sorting_algorithm == 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param weighted Logical. Weighted objective
#' @param cycle_start_positions Set. Cycle start positions to consider. Anything outside this set will be skipped. Only applies when \code{sorting_algorithm == 'tsp'}.
#' @param fixed_column Character or Integer. Name or position of the column in \code{graphing_columns} to keep fixed during sorting. Only applies when \code{sorting_algorithm == 'greedy_wolf'}.
#' @param random_initializations Integer. Number of random initializations for the positions of each grouping in \code{graphing_columns}. Only applies when \code{sorting_algorithm == 'greedy_wolf' or sorting_algorithm == 'greedy_wblf'}.
#' @param output_df_path Optional character. Output path for the output data frame, in CSV format. If \code{NULL}, then will not be saved.
#' @param preprocess_data Logical. If TRUE, will preprocess the data with the \code{data_preprocess} function.
#' @param default_sorting Character. Default column sorting in data_preprocess if integer columns do not exist. Options are 'alphabetical' (default), 'reverse_alphabetical', 'increasing', 'decreasing', 'random'.
#' @param return_updated_graphing_columns Logical. If FALSE, will only return the updated data frame. If TRUE, will return both the updated data frame and the updated graphing_columns parameter in the order in which the columns should be graphed.
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param print_params Logical. If TRUE, will print function params.
#' @param load_df Internal flag; not recommended to modify.
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
data_sort <- function(df, graphing_columns = NULL, column1 = NULL, column2 = NULL, column_weights = NULL, sorting_algorithm = "tsp", optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE, matrix_initialization_value = 1e6, same_side_matrix_initialization_value = 1e6, weight_scalar = 5e5, matrix_initialization_value_column_order = 1e6, weight_scalar_column_order = 1, column_sorting_metric = "edge_crossing", column_sorting_algorithm = "tsp", weighted = TRUE, cycle_start_positions = NULL, fixed_column = NULL, random_initializations = 1, output_df_path = NULL, preprocess_data = TRUE, default_sorting = "alphabetical", return_updated_graphing_columns = FALSE, verbose = FALSE, print_params = FALSE, load_df = TRUE, do_compute_alluvial_statistics = TRUE) {
    if (print_params) print_function_params()
    lowercase_args(c("sorting_algorithm", "column_sorting_metric", "column_sorting_algorithm", "default_sorting"))
    
    #* Type Checking Start
    valid_algorithms <- c("tsp", "neighbornet", "greedy_wolf", "greedy_wblf", "none", "random")
    if (!(sorting_algorithm %in% valid_algorithms)) {
        stop(sprintf(
            "Invalid sorting_algorithm: '%s'. Must be one of: %s",
            sorting_algorithm, paste(valid_algorithms, collapse = ", ")
        ))
    }
    
    if (sorting_algorithm == "tsp" || sorting_algorithm == "neighbornet") {
        for (col in graphing_columns) {
            if (grepl("~~", col)) {
                stop(sprintf("No entry of graphing_columns can contain '~~' when sorting_algorithm == tsp. Issue with column '%s'.", col))
            }
        }
    }
    
    if (((sorting_algorithm == "none") || (sorting_algorithm == "random") || sorting_algorithm == "neighbornet" || (sorting_algorithm == "tsp")) && (random_initializations > 1)) {
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
        clus_df_gather_sorted <- sort_neighbornet(clus_df_gather = clus_df_gather, graphing_columns = graphing_columns, column_weights = column_weights, optimize_column_order = optimize_column_order, optimize_column_order_per_cycle = optimize_column_order_per_cycle, matrix_initialization_value = matrix_initialization_value, same_side_matrix_initialization_value = same_side_matrix_initialization_value, weight_scalar = weight_scalar, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_sorting_metric = column_sorting_metric, sorting_algorithm = sorting_algorithm, column_sorting_algorithm = column_sorting_algorithm, cycle_start_positions = cycle_start_positions, verbose = verbose, weighted = weighted)
    } else if (sorting_algorithm == "greedy_wblf" || sorting_algorithm == "greedy_wolf") {
        # O(n_1 * n_2) complexity, where n1 is the number of blocks in layer 1, and n2 is the number of blocks in layer 2
        clus_df_gather_sorted <- sort_greedy_wolf(clus_df_gather = clus_df_gather, graphing_columns = graphing_columns, column1 = column1, column2 = column2, column_weights = column_weights, fixed_column = fixed_column, random_initializations = random_initializations, sorting_algorithm = sorting_algorithm, verbose = verbose, weighted = weighted)
    } else if (sorting_algorithm == "random") {
        clus_df_gather_sorted <- randomly_map_int_columns(clus_df_gather)
        #!!! check this
        if (optimize_column_order) {
            graphing_columns_tmp <- determine_column_order(clus_df_gather_sorted, graphing_columns = graphing_columns, column_weights = column_weights, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_sorting_metric = column_sorting_metric, column_sorting_algorithm = column_sorting_algorithm, verbose = verbose, weighted = weighted)
            clus_df_gather_tmp <- reorder_and_rename_columns(clus_df_gather_sorted, graphing_columns_tmp)
            # make factors
            for (j in seq_along(graphing_columns_tmp)) {
                int_col_name <- paste0("col", j, "_int")
                clus_df_gather_tmp[[int_col_name]] <- factor(clus_df_gather_tmp[[int_col_name]])
            }
            clus_df_gather_sorted <- clus_df_gather_tmp
        }
        #!!! check this
    } else if (sorting_algorithm == "none") {
        clus_df_gather_sorted <- clus_df_gather
        #!!! check this
        if (optimize_column_order) {
            graphing_columns_tmp <- determine_column_order(clus_df_gather_sorted, graphing_columns = graphing_columns, column_weights = column_weights, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_sorting_metric = column_sorting_metric, column_sorting_algorithm = column_sorting_algorithm, verbose = verbose, weighted = weighted)
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
        stop(sprintf("Invalid sorting_algorithm: '%s'. Must be one of: %s", sorting_algorithm, paste(valid_algorithms, collapse = ", ")))
    }
    
    # print objective - don't do for neighbornet because I did it right before
    if ((verbose) && (sorting_algorithm != "neighbornet")) {
        message("Determining crossing edges objective (to disable, use verbose==FALSE)")
        objective <- determine_crossing_edges(clus_df_gather_sorted, graphing_columns = graphing_columns, column_weights = column_weights, weighted = weighted)$output_objective
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
