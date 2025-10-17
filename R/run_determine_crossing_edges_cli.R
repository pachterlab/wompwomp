#' @noRd
run_determine_crossing_edges_cli <- function(args) {
    if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
        cat("
Usage: wompwomp determine_crossing_edges --input INPUT [options]

Required:
  -i, --input, --df               A data frame, tibble, or CSV file path. Must be in one of two formats:
(1) column_weights == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two graphing_columns).
(2) column_weights != NULL: Each row represents a combination of groupings, each column from graphing_columns represents a grouping, and the column column_weights represents the number of entities in that combination of groupings. Must contain at least three columns (two graphing_columns, one column_weights).

Optional:
  -g, --graphing_columns    Vector of column names from df to be used in graphing (i.e., alluvial plotting). Mutually exclusive with column1 and column2.
  -c1, --column1            Can be used along with column2 in place of graphing_columns if working with two columns only. Mutually exclusive with graphing_columns.
  -c2, --column2            Can be used along with column1 in place of graphing_columns if working with two columns only. Mutually exclusive with graphing_columns.
  -w, --column_weights      Column name from df that contains the weights of each combination of groupings if df is in format (2) (see above).
  --normalize_objective  Logical. Whether to normalize the objective by dividing by the sum of products of all edge weights.
  -o, --output_df_path      Output path for the output data frame, in CSV format. If NULL, then will not be saved.
  --output_lode_df_path     Output path for the data frame containing lode information on each alluvium, in CSV format (see details below). If not provided, then will not be saved.
  --include_output_objective_matrix_vector            Whether to return a vector of matrices, where each matrix is square with dimension equal to the number of alluvia, and where entry (i,j) of a matrix represents the product of weights of alluvium i and alluvium j if they cross, and 0 otherwise. There are (n-1) matrices in the vector, where n is the length of graphing_columns.
  --return_weighted_layer_free_objective  Whether to return a list of overlapping edges (FALSE) or the sum of products of overlapping edges (TRUE)
  --disable_use_fenwick_tree_for_objective_calculation  Whether to use fenwick trees for objective calculation. Speeds up from O(n^2) to O(nlogn).
  --print_params            If TRUE, will print function params.
  -v, --verbose             If TRUE, will display messages during the function.
  -q, --quiet               Don't show stdout
")
        quit(save = "no", status = 0)
    }

    # Required arguments
    df <- get_arg(args, c("-i", "--input", "--df"), required = TRUE)

    # Optional arguments
    graphing_columns <- get_multi_arg(args, c("-g", "--graphing_columns"))
    column1 <- get_arg(args, c("-c1", "--column1"))
    column2 <- get_arg(args, c("-c2", "--column2"))
    column_weights <- get_arg(args, c("-w", "--column_weights"))
    normalize_objective <- store_true(args, c("-normalize_objective"))
    output_df_path <- get_arg(args, c("-o", "--output_df_path"))
    output_lode_df_path <- get_arg(args, c("--output_lode_df_path"))
    include_output_objective_matrix_vector <- store_true(args, "--include_output_objective_matrix_vector")
    return_weighted_layer_free_objective <- store_true(args, "--return_weighted_layer_free_objective")
    use_fenwick_tree_for_objective_calculation <- store_false(args, "--disable_use_fenwick_tree_for_objective_calculation")
    verbose <- store_true(args, c("-v", "--verbose"))
    print_params <- store_true(args, c("--print_params"))
    quiet <- store_true(args, c("-q", "--quiet"))

    # Hidden arguments
    preprocess_data <- store_false(args, c("--disable_preprocess_data"))
    load_df <- store_false(args, c("--disable_load_df"))
    default_sorting <- get_arg(args, c("--default_sorting"))

    # Base argument list with required args
    args_list <- list(
        df = df
    )

    # Conditionally add optional args if not NULL
    if (!is.null(graphing_columns)) args_list$graphing_columns <- graphing_columns
    if (!is.null(column1)) args_list$column1 <- column1
    if (!is.null(column2)) args_list$column2 <- column2
    if (!is.null(column_weights)) args_list$column_weights <- column_weights
    if (!is.null(normalize_objective)) args_list$normalize_objective <- normalize_objective
    if (!is.null(output_df_path)) args_list$output_df_path <- output_df_path
    if (!is.null(output_lode_df_path)) args_list$output_lode_df_path <- output_lode_df_path
    if (!is.null(include_output_objective_matrix_vector)) args_list$include_output_objective_matrix_vector <- include_output_objective_matrix_vector
    if (!is.null(return_weighted_layer_free_objective)) args_list$return_weighted_layer_free_objective <- return_weighted_layer_free_objective
    if (!is.null(use_fenwick_tree_for_objective_calculation)) args_list$use_fenwick_tree_for_objective_calculation <- use_fenwick_tree_for_objective_calculation
    if (!is.null(verbose)) args_list$verbose <- verbose
    if (!is.null(print_params)) args_list$print_params <- print_params
    if (!is.null(preprocess_data)) args_list$preprocess_data <- preprocess_data
    if (!is.null(load_df)) args_list$load_df <- load_df
    if (!is.null(default_sorting)) args_list$default_sorting <- default_sorting

    # Dynamically call function
    result <- do.call(determine_crossing_edges, args_list)

    # Output logic
    if (!quiet) {
        cat("Sum of products of overlapping edge weights:\n")
        print(result$output_objective)
        if (!return_weighted_layer_free_objective) {
            cat("Crossing edges data frame:\n")
            print(result$crossing_edges_df)
        }
    }
}
