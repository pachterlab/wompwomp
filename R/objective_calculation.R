#' alluvialmatch: Cluster-matching alluvial plots
#'
#' Main plotting function and helpers for bipartite-matching-based alluvial diagrams
#' @docType package
#' @name alluvialmatch
#'
#' @importFrom dplyr mutate group_by summarise arrange desc ungroup slice n pull filter row_number left_join rename_with
#' @importFrom tidyr pivot_wider
#' @importFrom tibble is_tibble
#' @importFrom utils read.csv
#' @importFrom rlang sym .data
#' @importFrom magrittr %>%
#' @importFrom data.table :=
#' @importFrom ggalluvial stat_alluvium
#' @importFrom ggplot2 ggplot_build

utils::globalVariables(c(
    ".data", ":=", "%>%", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "total", "cum_y", "best_cluster_agreement"
))


#' Compute crossing objective
#'
#' Determine the sum of products of overlapping edge weights.
#'
#' @param df A CSV path or data frame as outputted with \code{crossing_edges_df} (in R) or \code{output_df_path} (as a file) from \code{determine_crossing_edges}.
#' @param verbose Logical. If TRUE, will display messages during the function.
#'
#' @return A non-negative integer.
#'
#' @examples
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- data_sort(df, graphing_columns = c('method1', 'method2'))
#' crossing_edges_output <- determine_crossing_edges(clus_df_gather, column1 = "method1", column2 = "method2")
#' objective <- determine_weighted_layer_free_objective(crossing_edges_output$crossing_edges_df)
#'
#' @export
determine_weighted_layer_free_objective <- function(df, verbose = FALSE) {
    # Case 1: CSV
    if (is.character(df) && length(df) == 1 && file.exists(df)) {
        # Read the file
        df <- read.csv(df, stringsAsFactors = FALSE)
    # Case 2: Already a list of edge pairs
    } else if (is.data.frame(df)) {
        # do nothing
    } else
        stop("Input must be either a file path or a list.")

    total_weighted_crossings <- sum(df$weight1 * df$weight2) / 2  # Correct for double-counting
    return(total_weighted_crossings)
}


#' Determine overlapping edges
#'
#' Determine overlapping edges of k-partite graph.
#'
#' @param df A data frame, tibble, or CSV file path. Must be in one of two formats:
#' (1) column_weights == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two graphing_columns).
#' (2) column_weights != NULL: Each row represents a combination of groupings, each column from \code{graphing_columns} represents a grouping, and the column \code{column_weights} represents the number of entities in that combination of groupings. Must contain at least three columns (two \code{graphing_columns}, one \code{column_weights}).
#' @param graphing_columns Optional character vector. Vector of column names from \code{df} to be used in graphing (i.e., alluvial plotting). Mutually exclusive with \code{column1} and \code{column2}.
#' @param column1 Optional character. Can be used along with \code{column2} in place of \code{graphing_columns} if working with two columns only. Mutually exclusive with \code{graphing_columns}.
#' @param column2 Optional character. Can be used along with \code{column1} in place of \code{graphing_columns} if working with two columns only. Mutually exclusive with \code{graphing_columns}.
#' @param column_weights Optional character. Column name from \code{df} that contains the weights of each combination of groupings if \code{df} is in format (2) (see above).
#' @param output_df_path Optional character. Output path for the data frame containing crossing edges, in CSV format (see details below). If not provided, then nothing will be saved.
#' @param output_lode_df_path Optional character. Output path for the data frame containing lode information on each alluvium, in CSV format (see details below). If not provided, then nothing will be saved.
#' @param include_output_objective_matrix_vector Logical. Whether to return a vector of matrices, where each matrix is square with dimension equal to the number of alluvia, and where entry (i,j) of a matrix represents the product of weights of alluvium i and alluvium j if they cross, and 0 otherwise. There are (n-1) matrices in the vector, where n is the length of graphing_columns.
#' @param return_weighted_layer_free_objective Logical. Whether to return a list of overlapping edges (FALSE) or the sum of products of overlapping edges (TRUE)
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param stratum_column_and_value_to_keep Internal flag; not recommended to modify.
#' @param input_objective_matrix_vector Internal flag; not recommended to modify.
#' @param input_objective Internal flag; not recommended to modify.
#'
#' @return
#' If return_weighted_layer_free_objective is FALSE (default): A list of values, as follows:
#' 'crossing_edges_df': A data frame containing the following columns:
#'   - alluvium1: The ID of the first alluvium, corresponding to the 'alluvium' column in \code{lode_df}.
#'   - alluvium2: The ID of the second alluvium, corresponding to the 'alluvium' column in \code{lode_df}.
#'   - strat_layer: The region in which the overlap occurred, corresponding to the 'xi' column in \code{lode_df}, where i is an integer 1, 2, ..., length(graphing_columns).
#'   - weight1: The weight of the first alluvium, corresponding to the 'count' column in \code{lode_df}.
#'   - weight2: The weight of the second alluvium, corresponding to the 'count' column in \code{lode_df}.
#' 'lode_df': A data frame containing the following columns:
#'   - alluvium: A specific alluvium/edge.
#'   - count: The weight of the alluvium/edge.
#'   - x1, x2, ...: Each xi represents the x position of axis/layer i.
#'   - y1, y2, ...: Each yi represents the height of a lode in axis/layer i.
#'   - stratum1, stratum2, ...: Each stratumi represents the stratum through which the alluvial crosses in axis/layer i.
#'   - weight1: The weight of the first alluvium, corresponding to the 'count' column in \code{lode_df}.
#'   - weight2: The weight of the second alluvium, corresponding to the 'count' column in \code{lode_df}.
#' 'output_objective': An integer representing the sum of products of overlapping edge weights.
#' 'objective_matrix_vector' (if and only if \code{include_output_objective_matrix_vector} is TRUE): A vector of square symmetric matrices. Each matrix in index h of the vector has rank equal to the number of alluvia present, where entry (i,j) represents the product of edge weights between alluvium i and alluvium j between layers h and h+1 (where the first layer has h=1).
#' If return_weighted_layer_free_objective is TRUE: An integer representing the sum of products of overlapping edge weights.
#'
#'
#' @examples
#' \dontrun{
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' df <- data_sort(df)
#' result <- determine_crossing_edges(df, column1 = "col1_int", column2 = "col2_int")
#' }
#'
#' @export
determine_crossing_edges <- function(df, graphing_columns = NULL, column1 = NULL, column2 = NULL, column_weights = "value", output_df_path = NULL, output_lode_df_path = NULL, include_output_objective_matrix_vector = FALSE, return_weighted_layer_free_objective = FALSE, verbose = FALSE, stratum_column_and_value_to_keep = NULL, input_objective_matrix_vector = NULL, input_objective = NULL, preprocess_data = TRUE, load_df = TRUE) {
    #* Type Checking Start
    # ensure someone doesn't specify both graphing_columns and column1/2
    if (!is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        stop("Specify either graphing_columns or column1/column2, not both.")
    }

    if (load_df) {
        df <- load_in_df(df = df, graphing_columns = graphing_columns, column_weights = column_weights)
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
    if (preprocess_data) {
        clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns, column_weights = column_weights, load_df = FALSE, do_gather_set_data = FALSE)
    } else {
        clus_df_gather <- df
    }

    p <- ggplot(data = clus_df_gather, aes(y = value),
    )
    for (x in seq_along(graphing_columns)) {
        int_col <- paste0('col', x,'_int')
        if (!(int_col %in% colnames(clus_df_gather))) {
            stop(sprintf("%s not in columns. Please run data_preprocess first.", int_col))
        }
        p$mapping[[paste0('axis',x)]] = sym(int_col)
    }
    p <- p + stat_alluvium(geom = "blank")

    columns_to_keep <- c("alluvium", "x", "y", "stratum", "count")
    lode_df_long_full <- ggplot_build(p)$data[[1]][columns_to_keep]

    # Initialize result list and seen pair tracker
    crossing_edges <- list()
    row_index <- 1

    if (is.null(input_objective_matrix_vector)) {
        output_objective_matrix_vector <- c()
    } else {
        output_objective_matrix_vector <- input_objective_matrix_vector
    }

    if (is.null(input_objective)) {
        output_objective <- 0
    } else {
        output_objective <- input_objective
    }

    # Get unique x values, sorted
    x_vals <- sort(unique(lode_df_long_full$x))
    n_x <- length(x_vals)

    # make the full lode_df
    lode_df_long_indexed_full <- lode_df_long_full %>%
        group_by(alluvium) %>%
        mutate(pos = row_number()) %>%
        ungroup()

    # Pivot each of x, y, stratum into wide format
    lode_df_full <- lode_df_long_indexed_full %>%
        select(alluvium, pos, x, y, stratum, count) %>%
        pivot_wider(
            id_cols = c(alluvium, count),
            names_from = pos,
            values_from = c(x, y, stratum),
            names_glue = "{.value}{pos}"
        )

    # add the actual character values
    for (i in seq_along(graphing_columns)) {
        int_col <- paste0("col", i, "_int")         # e.g. col1_int
        label_col <- graphing_columns[i]
        stratum_col <- paste0("stratum", i)         # e.g. stratum1
        stratum_char_col <- paste0(stratum_col, "_char")  # e.g. stratum1_char

        mapping <- setNames(clus_df_gather[[label_col]], clus_df_gather[[int_col]])
        mapping <- mapping[!duplicated(names(mapping))]
        lode_df_full[[stratum_char_col]] <- mapping[as.character(lode_df_full[[stratum_col]])]
    }

    if (!is.null(stratum_column_and_value_to_keep)) {
        layer_number <- as.integer(names(stratum_column_and_value_to_keep)[1])  # the layer (eg 3 from stratum3)
        stratum_number <- stratum_column_and_value_to_keep[[1]]  # the stratum (eg value 28 in column stratum3)

        alluvium_values_in_stratum_to_keep <- lode_df_long_full %>%
            filter(x == layer_number, stratum == stratum_number) %>%
            pull(alluvium) %>%
            unique()
    }

    if (verbose) message("Beginning loop through layers")
    for (h in 1:(n_x - 1)) {
        x1 <- h
        x2 <- h + 1

        # if (!is.null(stratum_column_and_value_to_keep)) {
        #     # even if I have 6 layers, this function will soon filter out all but x1 and x2 - also, if I am doing the whole matrix thing, then I can rest easy that other layers won't matter
        #     if (layer_number == x1) {
        #         layer_number_in_lode_df <- 1
        #     } else if (layer_number == x2) {
        #         layer_number_in_lode_df <- 2
        #     } else {
        #         next
        #     }
        # }

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

        x_vec <- c(x1, x2)
        for (i in seq_along(x_vec)) {
            x <- x_vec[i]
            stratum_col <- paste0("stratum", x)
            stratum_char_col <- paste0("stratum", x, "_char")
            mapping <- setNames(
                lode_df_full[[stratum_char_col]],
                lode_df_full[[stratum_col]]
            )
            stratum_col_for_lode_df <- paste0("stratum", i)
            stratum_char_col_for_lode_df <- paste0("stratum", i, "_char")
            lode_df[[stratum_char_col_for_lode_df]] <- mapping[as.character(lode_df[[stratum_col_for_lode_df]])]
        }

        # stratum_column_and_value_to_keep could be like list("3" = 28), where 3 is the x/layer number and 28 is the stratum number
        if (!is.null(stratum_column_and_value_to_keep)) {
            # stratum_column_name <- paste0("stratum", layer_number_in_lode_df)
            lode_df_filtered_with_stratum_of_interest <- lode_df[lode_df$alluvium %in% alluvium_values_in_stratum_to_keep, ]
            lode_df_filtered_with_stratum_of_interest_length <- nrow(lode_df_filtered_with_stratum_of_interest)

            lode_df_filtered_without_stratum_of_interest <- lode_df[!lode_df$alluvium %in% alluvium_values_in_stratum_to_keep, ]
            lode_df_filtered_without_stratum_of_interest_length <- nrow(lode_df_filtered_without_stratum_of_interest)
        }

        lode_df_length <- nrow(lode_df)

        objective_matrix <- NULL
        if (include_output_objective_matrix_vector) {
            if (is.null(input_objective_matrix_vector)) {
                ids <- sort(unique(lode_df$alluvium))
                # number_alluvia <- length(ids)
                objective_matrix <- matrix(0, nrow = lode_df_length, ncol = lode_df_length, dimnames = list(as.character(ids), as.character(ids)))
            } else {
                objective_matrix <- input_objective_matrix_vector[[h]]
            }
        }

        # browser()

        # Compare each pair of edges
        if (verbose) message("Looping through alluvia")
        if (is.null(stratum_column_and_value_to_keep)) {
            for (i in 1:(lode_df_length - 1)) {
                for (j in (i + 1):lode_df_length) {
                    # Check crossing condition once per unordered pair
                    if ((lode_df$y1[i] - lode_df$y1[j]) * (lode_df$y2[i] - lode_df$y2[j]) < 0) {
                    # equivalent to below but faster
                    # if ((lode_df$y1[i] < lode_df$y1[j] && lode_df$y2[i] > lode_df$y2[j]) | (lode_df$y1[i] > lode_df$y1[j] && lode_df$y2[i] < lode_df$y2[j])) {
                        alluvium1 <- lode_df$alluvium[i]
                        alluvium2 <- lode_df$alluvium[j]
                        w1 <- lode_df$count[i]
                        w2 <- lode_df$count[j]
                        strat_layer <- lode_df$x1[i]  # will be same for i and j

                        # stratum1_left <- lode_df$stratum1_char[i]
                        # stratum1_right <- lode_df$stratum2_char[i]  # apologies for the i/j and 1/2 confusion - this is correct though
                        # stratum2_left <- lode_df$stratum1_char[j]
                        # stratum2_right <- lode_df$stratum2_char[j]

                        # Append (i, j)
                        new_row <- data.frame(alluvium1 = alluvium1, alluvium2 = alluvium2, strat_layer = strat_layer, weight1 = w1, weight2 = w2)
                        # new_row <- data.frame(alluvium1 = alluvium1, alluvium2 = alluvium2, strat_layer = strat_layer, stratum1_left = stratum1_left, stratum1_right = stratum1_right, stratum2_left = stratum2_left, stratum2_right = stratum2_right, weight1 = w1, weight2 = w2)
                        crossing_edges[[row_index]] <- new_row
                        row_index <- row_index + 1

                        # Append (j, i)
                        new_row <- data.frame(alluvium1 = alluvium2, alluvium2 = alluvium1, strat_layer = strat_layer, weight1 = w2, weight2 = w1)
                        # new_row <- data.frame(alluvium1 = alluvium2, alluvium2 = alluvium1, strat_layer = strat_layer, stratum1_left = stratum2_left, stratum1_right = stratum2_right, stratum2_left = stratum1_left, stratum2_right = stratum1_right, weight1 = w2, weight2 = w1)
                        crossing_edges[[row_index]] <- new_row
                        row_index <- row_index + 1

                        weight_product <- w1 * w2
                        output_objective <- output_objective + weight_product

                        if (include_output_objective_matrix_vector) {
                            objective_matrix[alluvium1, alluvium2] <- weight_product
                            objective_matrix[alluvium2, alluvium1] <- weight_product
                        }
                    }
                }
            }
        } else {
            if (!include_output_objective_matrix_vector) {
                stop("If stratum_column_and_value_to_keep is not NULL, then include_output_objective_matrix_vector must be provided")
            }
            for (i in 1:(lode_df_filtered_with_stratum_of_interest_length)) {
                for (j in 1:lode_df_filtered_without_stratum_of_interest_length) {
                    alluvium1 <- lode_df_filtered_with_stratum_of_interest$alluvium[i]
                    alluvium2 <- lode_df_filtered_without_stratum_of_interest$alluvium[j]
                    w1 <- lode_df_filtered_with_stratum_of_interest$count[i]
                    w2 <- lode_df_filtered_without_stratum_of_interest$count[j]

                    if ((lode_df_filtered_with_stratum_of_interest$y1[i] - lode_df_filtered_without_stratum_of_interest$y1[j]) * (lode_df_filtered_with_stratum_of_interest$y2[i] - lode_df_filtered_without_stratum_of_interest$y2[j]) < 0) {
                        # crosses now, but didn't before (most cases)
                        if (objective_matrix[alluvium1, alluvium2] == 0) {
                            weight_product <- w1 * w2
                            objective_matrix[alluvium1, alluvium2] <- weight_product
                            objective_matrix[alluvium2, alluvium1] <- weight_product
                            output_objective <- output_objective + weight_product
                        }
                    } else {
                        # didn't cross before, but crosses now (most cases)
                        if (objective_matrix[alluvium1, alluvium2] > 0) {
                            weight_product <- w1 * w2
                            objective_matrix[alluvium1, alluvium2] <- 0
                            objective_matrix[alluvium2, alluvium1] <- 0
                            output_objective <- output_objective - weight_product
                        }
                    }
                }
            }
        }

        if (include_output_objective_matrix_vector) {
            if (is.null(input_objective_matrix_vector)) {
                output_objective_matrix_vector <- c(output_objective_matrix_vector, list(objective_matrix))
            } else {
                output_objective_matrix_vector[[h]] <- objective_matrix
            }
        }
        if (verbose) message(sprintf("Complete with iteration=%s", h))
    }

    if (length(crossing_edges) > 0) {
        crossing_edges_df <- do.call(rbind, crossing_edges)
    } else {
        crossing_edges_df <- data.frame(alluvium1 = numeric(), alluvium2 = numeric(), strat_layer = numeric(), weight1 = numeric(), weight2 = numeric())
        # crossing_edges_df <- data.frame(alluvium1 = numeric(), alluvium2 = numeric(), strat_layer = numeric(), stratum1_left = character(), stratum1_right = character(), stratum2_left = character(), stratum2_right = character(), weight1 = numeric(), weight2 = numeric())
    }

    if (is.character(output_df_path) && grepl("\\.csv$", output_df_path, ignore.case = TRUE)) {
        if (verbose) message(sprintf("Saving crossing_edges_df dataframe to=%s", output_df_path))
        write.csv(crossing_edges_df, file = output_df_path, row.names = FALSE, quote = FALSE)

        if (is.null(output_lode_df_path)) {
            output_lode_df_path <- sub("\\.csv$", "_lodes.csv", output_df_path)
        }
        if (verbose) message(sprintf("Saving lode_df dataframe to=%s", output_lode_df_path))
        write.csv(lode_df_full, file = output_lode_df_path, row.names = FALSE, quote = FALSE)
    }

    if (return_weighted_layer_free_objective) {
        # if (verbose) message("Calculating WLF objective")
        # output_objective <- determine_weighted_layer_free_objective(crossing_edges_df)  # already determined now during loop
        return(output_objective)
    }

    # only do this if I'm not in my neighbornet loop
    if (is.null(stratum_column_and_value_to_keep)) {
        crossing_edges_df <- crossing_edges_df %>%
            left_join(
                lode_df_full %>%
                    select(alluvium, matches("^stratum.*char$")) %>%
                    rename_with(~ paste0(., "1"), .cols = matches("^stratum.*char$")),
                by = c("alluvium1" = "alluvium")
            )
        crossing_edges_df <- crossing_edges_df %>%
            left_join(
                lode_df_full %>%
                    select(alluvium, matches("^stratum.*char$")) %>%
                    rename_with(~ paste0(., "2"), .cols = matches("^stratum.*char$")),
                by = c("alluvium2" = "alluvium")
            )
    }

    if (isTRUE(include_output_objective_matrix_vector)) {
        return(list(crossing_edges_df = crossing_edges_df, lode_df = lode_df_full, output_objective = output_objective, objective_matrix_vector = output_objective_matrix_vector))
    } else {
        return(list(crossing_edges_df = crossing_edges_df, lode_df = lode_df_full, output_objective = output_objective))
    }
}
