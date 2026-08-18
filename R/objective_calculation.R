#' wompwomp: Cluster-matching alluvial plots
#'
#' @name wompwomp-imports
#' @rdname wompwomp
#' @importFrom dplyr mutate group_by ungroup n row_number left_join
#' @importFrom rlang sym .data
#'
utils::globalVariables(c(
    ".data", ":=", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "total", "cum_y", "best_cluster_agreement"
))

# ---- Binary Indexed Tree (Fenwick Tree) ----
# (An earlier R6-class version of this function was ~5x slower end-to-end due
# to R6 method-dispatch overhead on the per-row update()/query() calls, which
# dominate runtime since this is called on every cycle-start iteration. The
# BIT build/query loop itself was later profiled at >95% of wall time on
# large sweeps -- bitwAnd() alone, called as a full R function per Fenwick
# step, was ~37% on its own -- so that loop now runs as compiled code; see
# calculate_objective_fenwick_cpp() in src/fenwick.cpp.)
calculate_objective_fenwick <- function(data, y1 = "y1", y2 = "y2", wt = 'value', weighted_metric = TRUE) {
    # Step 1: Sort by y1 (only the two columns actually needed below; avoids
    # reordering every column of `data`, which is a tibble and so pays
    # tibble-subsetting dispatch overhead on top of the copy).
    n <- nrow(data)
    ord <- order(data[[y1]])
    y2v <- data[[y2]][ord]

    # Step 2: Rank-compress y2 (higher y2 → higher rank)
    y2_rank <- match(y2v, sort(unique(y2v)))
    weight_vec <- if (weighted_metric) data[[wt]][ord] else rep(1.0, n)

    calculate_objective_fenwick_cpp(as.integer(y2_rank), as.numeric(weight_vec), weighted_metric)
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

    n <- nrow(lode_df)
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
        ord <- order(lode_df[[i]], lode_df[[tiebreaker]])
        y_vals <- cumsum(lode_df[[wt]][ord])

        # `alluvium` is exactly 1:n, so `ord` is a permutation of the row
        # indices; invert it to place each cumulative value back at its
        # original row directly, instead of a dplyr::left_join keyed on that
        # same permutation.
        inv_ord <- integer(n)
        inv_ord[ord] <- seq_len(n)
        lode_df[[paste0('y', x)]] <- y_vals[inv_ord]
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
#             stop(sprintf("%s not in columns. Please run prep_for_lodes first.", int_col))
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
#' data <- sort_to_uncross(data, cols = c("method1", "method2"), method = "tsp")
#' result <- compute_crossing_objective(data, cols = c("method1", "method2"))
#'
#' @export
compute_crossing_objective <- function(data, cols = NULL, wt = "value", weighted_metric = TRUE, verbose = FALSE) {
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