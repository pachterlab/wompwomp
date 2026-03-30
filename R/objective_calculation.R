#' wompwomp: Cluster-matching alluvial plots
#'
#' @name wompwomp-imports
#' @rdname wompwomp
#' @importFrom dplyr mutate group_by ungroup n row_number left_join
#' @importFrom rlang sym .data
#' @importFrom R6 R6Class
#' 
utils::globalVariables(c(
    ".data", ":=", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "total", "cum_y", "best_cluster_agreement"
))

# ---- Binary Indexed Tree (Fenwick Tree) ----
BIT <- R6::R6Class("BIT",
                   public = list(
                       tree_count = NULL,
                       tree_weight = NULL,
                       
                       initialize = function(size) {
                           self$tree_count <- integer(size + 1)
                           self$tree_weight <- numeric(size + 1)
                       },
                       
                       update = function(index, weight) {
                           index <- index + 1L
                           n <- length(self$tree_count)
                           while (index < n) {
                               self$tree_count[index] <- self$tree_count[index] + 1L
                               self$tree_weight[index] <- self$tree_weight[index] + weight
                               index <- index + bitwAnd(index, -index)
                           }
                       },
                       
                       query = function(index) {
                           count <- 0L
                           weight_sum <- 0.0
                           index <- index + 1L
                           while (index > 0L) {
                               count <- count + self$tree_count[index]
                               weight_sum <- weight_sum + self$tree_weight[index]
                               index <- index - bitwAnd(index, -index)
                           }
                           list(count = count, weight_sum = weight_sum)
                       },
                       
                       query_range = function(low, high) {
                           q_high <- self$query(high)
                           q_low <- self$query(low - 1L)
                           list(
                               count = q_high$count - q_low$count,
                               weight_sum = q_high$weight_sum - q_low$weight_sum
                           )
                       }
                   )
)

calculate_objective_fenwick <- function(data, y1 = "y1", y2 = "y2", wt = 'value', weighted_metric = TRUE) {
    # Step 1: Sort by y1
    df_sorted <- data[order(data[[y1]]), ]
    rownames(df_sorted) <- NULL
    
    # Step 2: Rank-compress y2 (higher y2 → higher rank)
    df_sorted$y2_rank <- match(df_sorted[[y2]], sort(unique(df_sorted[[y2]])))
    max_rank <- max(df_sorted$y2_rank)
    
    # Step 3: Initialize BIT
    bit <- BIT$new(size = max_rank + 2L)
    total_cross_weight <- 0.0
    
    for (i in seq_len(nrow(df_sorted))) {
        y2_rank <- df_sorted$y2_rank[i]
        weight <- if (weighted_metric) df_sorted[[wt]][i] else 1.0
        
        # Count previous y2s > current (strictly greater)
        q <- bit$query_range(y2_rank + 1L, max_rank)
        
        if (weighted_metric) {
            total_cross_weight <- total_cross_weight + (weight * q$weight_sum)
        } else {
            total_cross_weight <- total_cross_weight + q$count
        }
        
        # Add current y2_rank to BIT
        bit$update(y2_rank, if (weighted_metric) weight else 1.0)
    }
    
    total_cross_weight
}


make_lode_df <- function(data, cols = NULL, wt = "value") {
    lode_df <- data
    lode_df$alluvium <- seq_len(nrow(lode_df))

    # create a temp column
    wt_was_null <- is.null(wt)
    if (wt_was_null) {
        lode_df$.wt_internal <- 1
        wt <- ".wt_internal"
    }

    n_cols <- length(cols)
    for (x in seq_len(n_cols)) {
        i <- cols[x]
        # Within each stratum, order edges by the adjacent axis so that edges
        # sharing a stratum are not counted as crossings by the Fenwick tree.
        # Use the next column as tiebreaker (or the previous for the last axis).
        if (x < n_cols) {
            tiebreaker <- cols[x + 1]
        } else {
            tiebreaker <- cols[x - 1]
        }
        ordered_df <- lode_df[order(lode_df[[i]], lode_df[[tiebreaker]]), ]
        ordered_df[[paste0('y', x)]] <- cumsum(ordered_df[[wt]])
        lode_df <- dplyr::left_join(
            lode_df,
            ordered_df[, c('alluvium', paste0('y', x))],
            by = 'alluvium'
        )
    }

    # remove temp column if created
    if (wt_was_null) {
        lode_df$.wt_internal <- NULL
    }

    return(lode_df)
}

# # if uncommenting, then move ggplot2 and ggalluvial from Suggests to Imports in DESCRIPTION
# make_lode_df_old <- function(data, cols = NULL, wt = "value") {
#     if (wt != "value") {
#         data <- data |> dplyr::rename(value = !!sym(wt))
#         wt <- "value"
#     }
#     
#     p <- ggplot2::ggplot(data = data, ggplot2::aes(y = value), )
#     for (x in seq_along(cols)) {
#         int_col <- paste0("col", x, "_int")
#         if (!(int_col %in% colnames(data))) {
#             stop(sprintf("%s not in columns. Please run data_preprocess first.", int_col))
#         }
#         p$mapping[[paste0("axis", x)]] <- sym(int_col)
#     }
#     p <- p + ggalluvial::stat_alluvium(geom = "blank")
#     columns_to_keep <- c("alluvium", "x", "y", "stratum", "count")
#     lode_df_long_full <- ggplot2::ggplot_build(p)$data[[1]][columns_to_keep]
#     
#     # Initialize result list and seen pair tracker
#     crossing_edges <- list()
#     row_index <- 1
#     
#     output_objective <- 0
#     
#     # Get unique x values, sorted
#     x_vals <- sort(unique(lode_df_long_full$x))
#     n_x <- length(x_vals)
#     
#     # make the full lode_df
#     lode_df_long_indexed_full <- lode_df_long_full |>
#         group_by(alluvium) |>
#         mutate(pos = row_number()) |>
#         ungroup()
#     
#     # Pivot each of x, y, stratum into wide format
#     lode_df_full <- lode_df_long_indexed_full |>
#         select(alluvium, pos, x, y, stratum, count) |>
#         tidyr::pivot_wider(
#             id_cols = c(alluvium, count),
#             names_from = pos,
#             values_from = c(x, y, stratum),
#             names_glue = "{.value}{pos}"
#         )
#     
#     # add the actual character values
#     for (i in seq_along(cols)) {
#         int_col <- paste0("col", i, "_int") # e.g. col1_int
#         label_col <- cols[i]
#         stratum_col <- paste0("stratum", i) # e.g. stratum1
#         stratum_char_col <- paste0(stratum_col, "_char") # e.g. stratum1_char
#         
#         mapping <- setNames(data[[label_col]], data[[int_col]])
#         mapping <- mapping[!duplicated(names(mapping))]
#         lode_df_full[[stratum_char_col]] <- mapping[as.character(lode_df_full[[stratum_col]])]
#     }
#     
#     lode_df_full <- lode_df_full |> dplyr::rename(value = "count")
#     return(lode_df_full)
# }


#' Determine overlapping edges
#'
#' Determine overlapping edges of k-partite graph.
#'
#' @param data A data frame, tibble, or CSV file path. Must be in one of two formats:
#' (1) wt == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two cols).
#' (2) wt != NULL: Each row represents a combination of groupings, each column from \code{cols} represents a grouping, and the column \code{wt} represents the number of entities in that combination of groupings. Must contain at least three columns (two \code{cols}, one \code{wt}).
#' @param cols Optional character vector. Vector of column names from \code{data} to be used in graphing (i.e., alluvial plotting). Mutually exclusive with \code{column1} and \code{column2}.
#' @param wt Optional character. Column name from \code{data} that contains the weights of each combination of groupings if \code{data} is in format (2) (see above). If null, then sets \code{weighted_metric} to FALSE.
#' @param weighted_metric Logical. Determines if the objective is total number of edge crossings (weighted_metric=FALSE) or sum of product of overlapping edge weights (weighted_metric=TRUE).
#' @param verbose Logical. If TRUE, will display messages during the function.
#'
#' @return
#' If return_weighted_layer_free_objective is FALSE (default): A list of values, as follows:
#' 'lode_df': A data frame containing the following columns:
#'   - alluvium: A specific alluvium/edge.
#'   - count: The weight of the alluvium/edge.
#'   - x1, x2, ...: Each xi represents the x position of axis/layer i.
#'   - y1, y2, ...: Each yi represents the height of a lode in axis/layer i.
#'   - stratum1, stratum2, ...: Each stratumi represents the stratum through which the alluvial crosses in axis/layer i.
#'   - weight1: The weight of the first alluvium, corresponding to the 'count' column in \code{lode_df}.
#'   - weight2: The weight of the second alluvium, corresponding to the 'count' column in \code{lode_df}.
#' 'output_objective': An integer representing the sum of products of overlapping edge weights.
#'
#' @examples
#' data <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' data <- data_sort(data, cols = c("method1", "method2"), method = "tsp")
#' result <- determine_crossing_edges(data, cols = c("method1", "method2"))
#'
#' @export
determine_crossing_edges <- function(data, cols = NULL, wt = "value", weighted_metric = TRUE, verbose = FALSE) {
    if (!is.null(wt)) {
        if (!(wt %in% names(data))) {
            stop(sprintf("Column '%s' (wt) not found in data.", wt))
        }
    } else {
        # if wt is null, then weighted_metric is FALSE
        weighted_metric <- FALSE
    }
    
    col_ints <- c()
    for (h in seq_len(length(cols))) {
        col_ints <- c(col_ints, paste0('col', h, '_int'))
    }
    
    # add int columns if needed
    n_present <- sum(col_ints %in% names(data))
    if (n_present == 0) {
        data <- generalized_make_int_columns(data, cols)
    } else if (n_present == length(col_ints)) {
    } else {
        stop("Some int columns are present, but not all.")
    }
    
    lode_df <- make_lode_df(data, col_ints, wt)
    objective_val <- 0
    for (h in seq_len(length(cols) - 1)) {
        y1 <- paste0('y', h)
        y2 <- paste0('y', h+1)
        
        objective_val <- objective_val + calculate_objective_fenwick(lode_df, y1 = y1, y2 = y2, wt = wt, weighted_metric = weighted_metric)
    }
    return(list(lode_df = lode_df, output_objective = objective_val))
}