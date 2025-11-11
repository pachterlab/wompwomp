#' wompwomp: Coloring alluvia
#'
#' Coloring alluvia
#' @docType package
#' @name wompwomp
#'
#' @importFrom igraph V cluster_louvain cluster_leiden E
#' @importFrom magrittr %>%
#' @importFrom utils read.csv write.csv

utils::globalVariables(c(
    ".data", ":=", "%>%", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "stratum", "total", "cum_y", "best_cluster_agreement", "neighbor_net", "alluvium", "pos", "count", "group1", "group2", "value", "group1_size", "group2_size", "weight", "parent", "group_name"
))

#' Color alluvia
#'
#' Colors a dataframe with the algorithm specified by \code{coloring_algorithm}.
#'
#' @param df A data frame, tibble, or CSV file path. Must be in one of two formats:
#' (1) column_weights == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two graphing_columns).
#' (2) column_weights != NULL: Each row represents a combination of groupings, each column from \code{graphing_columns} represents a grouping, and the column \code{column_weights} represents the number of entities in that combination of groupings. Must contain at least three columns (two \code{graphing_columns}, one \code{column_weights}).
#' @param graphing_columns Character vector. Vector of column names from \code{df} to be used in graphing (i.e., alluvial plotting).
#' @param column_weights Optional character. Column name from \code{df} that contains the weights of each combination of groupings if \code{df} is in format (2) (see above).
#' @param color_boxes Logical. Whether to color the strata/boxes (representing groups).
#' @param color_list Optional named list or vector of colors to override default group colors.
#' @param color_val Optional named list where the entries are colors and the names correspond to values of the dataframe that should use those colors
#' @param coloring_algorithm Character. Matching colors methods. Choices are 'advanced' (default), 'none', 'left', 'right', or any value in \code{graphing_columns}.
#' @param coloring_algorithm_advanced_option Character. If \code{coloring_algorithm == 'advanced'}, then choose graph clustering algorithm. Choices are 'leiden' (default) or 'louvain'.
#' @param resolution Numeric If \code{coloring_algorithm == 'advanced'}, then choose resolution for the graph clustering algorithm.
#' @param cutoff Numeric If \code{coloring_algorithm != 'none' and coloring_algorithm != 'advanced'}, sets the cutoff for color matching, below which a new color will be assigned.
#' @param preprocess_data Logical. If TRUE, will preprocess the data with the \code{data_preprocess} function.
#' @param output_df_path Optional character. Output path for the output data frame, in CSV format. If \code{NULL}, then will not be saved.
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param print_params Logical. If TRUE, will print function params.
#'
#' @return A data frame where each row represents a combination of groupings, each column from \code{graphing_columns} represents a grouping, and the column \code{column_weights} ('value' if \code{column_weights} == NULL) represents the number of entities in that combination of groupings. There will be an additional column 'node_color'.
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
data_color <- function(df, graphing_columns, column_weights = NULL, color_list = NULL, color_val = NULL, coloring_algorithm = "advanced", coloring_algorithm_advanced_option = "leiden", resolution = 1, cutoff = .5, preprocess_data = TRUE, output_df_path = NULL, verbose = FALSE, print_params = FALSE, load_df = TRUE, do_compute_alluvial_statistics = TRUE) {
    if (print_params) print_function_params()
    lowercase_args(c("coloring_algorithm", "coloring_algorithm_advanced_option"))
    
    #* Type Checking Start
    valid_algorithms <- c("advanced", "left", "right", "none", "random", graphing_columns)
    if (!(sorting_algorithm %in% valid_algorithms)) {
        stop(sprintf(
            "Invalid sorting_algorithm: '%s'. Must be one of: %s",
            sorting_algorithm, paste(valid_algorithms, collapse = ", ")
        ))
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
    #* Type Checking End
    
    
    # Preprocess (i.e., add int columns and do the grouping)
    if (preprocess_data) {
        if (verbose) message("Preprocessing data before sorting")
        clus_df_gather <- data_preprocess(df = df, graphing_columns = graphing_columns, column_weights = column_weights, load_df = load_df, do_gather_set_data = FALSE, do_add_int_columns = FALSE)
        if (is.null(column_weights)) {
            column_weights <- "value" # is set during data_preprocess
        }
    } else {
        clus_df_gather <- df
    }
    
    #!!!! do coloring here
    
    # Save if desired
    if ((is.character(output_df_path) && grepl("\\.csv$", output_df_path, ignore.case = TRUE))) {
        if (verbose) message(sprintf("Saving sorted dataframe to=%s", output_df_path))
        write.csv(clus_df_gather_sorted, file = output_df_path, row.names = FALSE)
    }
    
    return (clus_df_gather)
}