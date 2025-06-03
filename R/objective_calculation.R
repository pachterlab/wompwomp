#' alluvialmatch: Cluster-matching alluvial plots
#'
#' Main plotting function and helpers for bipartite-matching-based alluvial diagrams
#' @docType package
#' @name alluvialmatch
#'
#' @importFrom dplyr mutate group_by summarise arrange desc ungroup slice n pull filter
#' @importFrom tidyr pivot_wider
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
#' crossing_edges_output <- determine_crossing_edges(df, column1 = "method1", column2 = "method2")
#' objective <- determine_weighted_layer_free_objective(crossing_edges_output$crossing_edges_df)
#' }
#'
#' @export
determine_weighted_layer_free_objective <- function(crossing_edges_df, verbose = FALSE) {
    # Case 1: CSV
    if (is.character(crossing_edges_df) && length(crossing_edges_df) == 1 && file.exists(crossing_edges_df)) {
        # Read the file
        df_in <- read.csv(crossing_edges_df, stringsAsFactors = FALSE)
    # Case 2: Already a list of edge pairs
    } else if (is.data.frame(crossing_edges_df)) {
        # do nothing
    } else
        stop("Input must be either a file path or a list.")

    total_weighted_crossings <- sum(crossing_edges_df$weight1 * crossing_edges_df$weight2) / 2  # Correct for double-counting
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
#' result <- determine_crossing_edges(df, column1 = "col1_int", column2 = "col2_int")
#' }
#'
#' @export
determine_crossing_edges <- function(df, graphing_columns = NULL, column1 = NULL, column2 = NULL, column_weights = "value", minimum_edge_weight = 0, output_df_path = NULL, output_lode_df_path = NULL,  return_weighted_layer_free_objective = FALSE, verbose = FALSE) {
    #* Type Checking Start
    # ensure someone doesn't specify both graphing_columns and column1/2
    if (!is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        stop("Specify either graphing_columns or column1/column2, not both.")
    }

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

    # # set to factors if not already
    # if (!is.factor(df[[column1]])) df[[column1]] <- factor(df[[column1]])
    # if (!is.factor(df[[column2]])) df[[column2]] <- factor(df[[column2]])

    if (verbose) message("Preprocessing data")
    clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns, column_weights = column_weights, load_df = FALSE, do_gather_set_data = FALSE)

    p <- ggplot(data = clus_df_gather, aes(y = value),
    )
    for (x in seq_along(graphing_columns)) {
        p$mapping[[paste0('axis',x)]] = sym(paste0('col', x,'_int'))
    }
    p <- p + stat_alluvium(geom = "blank")

    columns_to_keep <- c("alluvium", "x", "y", "stratum", "count")
    lode_df_long_full <- ggplot_build(p)$data[[1]][columns_to_keep]

    # Initialize result list and seen pair tracker
    crossing_edges <- list()
    row_index <- 1

    # Get unique x values, sorted
    x_vals <- sort(unique(lode_df_long_full$x))
    n_x <- length(x_vals)

    if (verbose) message("Beginning loop through layers")
    for (h in 1:(n_x - 1)) {
        x1 <- h
        x2 <- h + 1

        lode_df_long <- lode_df_long_full %>% filter(x == x1 | x == x2)

        # Sort by alluvium and x to ensure consistent ordering
        if (verbose) message("Filtering, grouping, and widening lode_df")
        lode_df_long_sorted <- lode_df_long %>%
            arrange(alluvium, x)

        # Create index within each alluvium group to track position (1, 2, ..., n)
        lode_df_long_indexed <- lode_df_long_sorted %>%
            group_by(alluvium) %>%
            mutate(pos = row_number()) %>%
            ungroup()

        # Pivot each of x, y, stratum into wide format
        lode_df <- lode_df_long_indexed %>%
            select(alluvium, pos, x, y, stratum, count) %>%
            pivot_wider(
                id_cols = c(alluvium, count),
                names_from = pos,
                values_from = c(x, y, stratum),
                names_glue = "{.value}{pos}"
            )

        # browser()

        lode_df_length <- nrow(lode_df)

        # Compare each pair of edges
        if (verbose) message("Looping through alluvia")
        for (i in 1:lode_df_length) {
            for (j in 1:lode_df_length) {
                # only look at rows where i's stratum is immediately adjacent to left of j's stratum
                if (i != j) {
                    if ((lode_df$y1[i] < lode_df$y1[j] && lode_df$y2[i] > lode_df$y2[j]) | (lode_df$y1[i] > lode_df$y1[j] && lode_df$y2[i] < lode_df$y2[j])) {
                        id1 <- lode_df$alluvium[i]
                        id2 <- lode_df$alluvium[j]

                        strat_layer <- lode_df$x1[i]  # will be same for i and j
                        # stratum1_left <- load_df$stratum1[i]
                        # stratum1_right <- load_df$stratum2[i]  # apologies for the i/j and 1/2 confusion - this is correct though
                        # stratum2_left <- load_df$stratum1[j]
                        # stratum2_right <- load_df$stratum2[j]
                        w1 <- lode_df$count[i]
                        w2 <- lode_df$count[j]

                        if (minimum_edge_weight > 0 && (w1 < minimum_edge_weight || w2 < minimum_edge_weight)) {
                            next
                        }

                        new_row <- data.frame(id1 = id1, id2 = id2, strat_layer = strat_layer, weight1 = w1, weight2 = w2)
                        crossing_edges[[row_index]] <- new_row
                        row_index <- row_index + 1
                    }
                }
            }
        }
        if (verbose) message(sprintf("Complete with iteration=%s", h))
    }

    if (length(crossing_edges) > 0) {
        crossing_edges_df <- do.call(rbind, crossing_edges)
    } else {
        crossing_edges_df <- data.frame(id1 = numeric(), id2 = numeric(), strat_layer = numeric(), weight1 = numeric(), weight2 = numeric())
    }

    if (is.character(output_df_path) && grepl("\\.csv$", output_df_path, ignore.case = TRUE)) {
        if (verbose) message(sprintf("Saving crossing_edges_df dataframe to=%s", output_df_path))
        write.csv(crossing_edges_df, file = output_df_path, row.names = FALSE, quote = FALSE)

        if (is.null(output_lode_df_path)) {
            output_lode_df_path <- sub("\\.csv$", "_lodes.csv", output_df_path)
        }
        if (verbose) message(sprintf("Saving lode_df dataframe to=%s", output_lode_df_path))
        write.csv(lode_df, file = output_lode_df_path, row.names = FALSE, quote = FALSE)
    }

    WLF_objective <- NULL
    if (return_weighted_layer_free_objective) {
        if (verbose) message("Calculating WLF objective")
        WLF_objective <- determine_weighted_layer_free_objective(crossing_edges_df)
        return(WLF_objective)
    }

    return(list(crossing_edges_df = crossing_edges_df, lode_df = lode_df))
}
