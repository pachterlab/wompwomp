#' @noRd
run_data_sort_cli <- function(args) {
    if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
        cat("
Usage: wompwomp data_sort --input INPUT [options]

Required:
  -i, --input, --df               A data frame, tibble, or CSV file path. Must be in one of two formats:
(1) column_weights == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two graphing_columns).
(2) column_weights != NULL: Each row represents a combination of groupings, each column from graphing_columns represents a grouping, and the column column_weights represents the number of entities in that combination of groupings. Must contain at least three columns (two graphing_columns, one column_weights).

Optional:
  -g, --graphing_columns    Vector of column names from df to be used in graphing (i.e., alluvial plotting). Mutually exclusive with column1 and column2.
  -c1, --column1            Can be used along with column2 in place of graphing_columns if working with two columns only. Mutually exclusive with graphing_columns.
  -c2, --column2            Can be used along with column1 in place of graphing_columns if working with two columns only. Mutually exclusive with graphing_columns.
  -w, --column_weights      Column name from df that contains the weights of each combination of groupings if df is in format (2) (see above).
  -s, --sorting_algorithm       Algorithm with which to sort the values in the dataframe. Can choose from {'tsp', 'greedy_wolf', 'greedy_wblf', 'none'. 'tsp' performs Traveling Salesman Problem solver from the TSP package. 'greedy_wolf' implements a custom greedy algorithm where one layer is fixed, and the other layer is sorted such that each node is positioned as close to its largest parent from the fixed side as possible in a greedy fashion. 'greedy_wblf' implements the 'greedy_wolf' algorithm described previously twice, treating each column as fixed in one iteration and free in the other iteration. 'greedy_wolf' and 'greedy_wblf' are only valid when graphing_columns has exactly two entries. 'random' randomly maps blocks. 'none' keeps the mappings as-is when passed into the function.
  --disable_optimize_column_order   If TRUE, will optimize the order of graphing_columns to minimize edge overlap. Only applies when sorting_algorithm == 'tsp' and length(graphing_columns) > 2.
  --optimize_column_order_per_cycle         If TRUE, will optimize the order of graphing_columns to minimize edge overlap upon each cycle. If FALSE, will optimize the order of graphing_columns to minimize edge overlap on the beginning cycle only. Only applies when sorting_algorithm == 'tsp' and length(graphing_columns) > 2.
  --matrix_initialization_value Positive integer. Initialized value in distance matrix for nodes in different layers without a shared edge/path. Only applies when sorting_algorithm == 'tsp'.
  --same_side_matrix_initialization_value Positive integer. Initialized value in distance matrix for nodes in the same layer. Only applies when sorting_algorithm == 'tsp'.
  --weight_scalar Positive integer. Scalar with which to multiply edge weights after taking their -log in the distance matrix for nodes with a nonzero edge. Only applies when sorting_algorithm == 'tsp'.
	--matrix_initialization_value_column_order Positive integer. Initialized value in distance matrix for optimizing column order. Only applies when sorting_algorithm == 'tsp' and optimize_column_order is TRUE.
	--weight_scalar_column_order Positive integer. Scalar with which to loss function after taking their log1p in the distance matrix for optimizing column order. Only applies when sorting_algorithm == 'tsp' and optimize_column_order is TRUE.
  --column_sorting_metric Character. Metric to use for determining column order. Options are 'edge_crossing' (default) or 'ARI'. Only applies when sorting_algorithm == 'tsp' and optimize_column_order is TRUE.
  --column_sorting_algorithm Character. Algorithm to use for determining column order. Options are 'tsp' (default). Only applies when sorting_algorithm == 'tsp' and optimize_column_order is TRUE.
  --cycle_start_positions Set. Cycle start positions to consider. Anything outside this set will be skipped. Only applies when sorting_algorithm == 'tsp'.
  --fixed_column            Name or position of the column in graphing_columns to keep fixed during sorting. Only applies when sorting_algorithm == 'greedy_wolf'.
  --random_initializations  Number of random initializations for the positions of each grouping in graphing_columns. Only applies when sorting_algorithm == 'greedy_wolf' or sorting_algorithm == 'greedy_wblf'.
  --set_seed                Random seed for the random_initializations parameter. Only applies when sorting_algorithm == 'greedy_wolf' or sorting_algorithm == 'greedy_wblf'.
  -o, --output_df_path      Output path for the output data frame, in CSV format. If NULL, then will not be saved.
  --disable_preprocess_data         If TRUE, will preprocess the data with the data_preprocess function.
  --default_sorting Character. Default column sorting in data_preprocess if integer columns do not exist. Will not affect output if sorting_algorithm == 'tsp'. Options are 'alphabetical' (default), 'reverse_alphabetical', 'increasing', 'decreasing', 'random'.
  --return_updated_graphing_columns         If FALSE, will only return the updated data frame. If TRUE, will return both the updated data frame and the updated graphing_columns parameter in the order in which the columns should be graphed.
  -v, --verbose             If TRUE, will display messages during the function.
  --print_params            If TRUE, will print function params.
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
    sorting_algorithm <- get_arg(args, c("-s", "--sorting_algorithm"))
    optimize_column_order <- store_false(args, c("--disable_optimize_column_order"))
    optimize_column_order_per_cycle <- store_true(args, c("--optimize_column_order_per_cycle"))
    matrix_initialization_value <- get_numeric_arg(args, c("--matrix_initialization_value"))
    same_side_matrix_initialization_value <- get_numeric_arg(args, c("--same_side_matrix_initialization_value"))
    weight_scalar <- get_numeric_arg(args, c("--weight_scalar"))
    matrix_initialization_value_column_order <- get_numeric_arg(args, c("--matrix_initialization_value_column_order"))
    weight_scalar_column_order <- get_numeric_arg(args, c("--weight_scalar_column_order"))
    column_sorting_metric <- get_arg(args, c("--column_sorting_metric"))
    column_sorting_algorithm <- get_arg(args, c("--column_sorting_algorithm"))
    cycle_start_positions <- get_multi_arg(args, c("--cycle_start_positions"))
    fixed_column <- get_fixed_column(args, "--fixed_column")
    random_initializations <- get_integer_arg(args, "--random_initializations")
    output_df_path <- get_arg(args, c("-o", "--output_df_path"))
    preprocess_data <- store_false(args, c("--disable_preprocess_data"))
    default_sorting <- get_arg(args, c("--default_sorting"))
    return_updated_graphing_columns <- store_true(args, "--return_updated_graphing_columns")
    verbose <- store_true(args, c("-v", "--verbose"))
    print_params <- store_true(args, c("--print_params"))
    quiet <- store_true(args, c("-q", "--quiet"))

    # Hidden arguments
    load_df <- store_false(args, c("--disable_load_df"))
    make_intermediate_neighbornet_plots <- store_true(args, c("make_intermediate_neighbornet_plots"))


    # Base argument list with required args
    args_list <- list(
        df = df
    )

    # Conditionally add optional args if not NULL
    if (!is.null(graphing_columns)) args_list$graphing_columns <- graphing_columns
    if (!is.null(column1)) args_list$column1 <- column1
    if (!is.null(column2)) args_list$column2 <- column2
    if (!is.null(column_weights)) args_list$column_weights <- column_weights
    if (!is.null(sorting_algorithm)) args_list$sorting_algorithm <- sorting_algorithm
    if (!is.null(optimize_column_order)) args_list$optimize_column_order <- optimize_column_order
    if (!is.null(optimize_column_order_per_cycle)) args_list$optimize_column_order_per_cycle <- optimize_column_order_per_cycle
    if (!is.null(matrix_initialization_value)) args_list$matrix_initialization_value <- matrix_initialization_value
    if (!is.null(same_side_matrix_initialization_value)) args_list$same_side_matrix_initialization_value <- same_side_matrix_initialization_value
    if (!is.null(weight_scalar)) args_list$weight_scalar <- weight_scalar
    if (!is.null(matrix_initialization_value_column_order)) args_list$matrix_initialization_value_column_order <- matrix_initialization_value_column_order
    if (!is.null(weight_scalar_column_order)) args_list$weight_scalar_column_order <- weight_scalar_column_order
    if (!is.null(column_sorting_metric)) args_list$column_sorting_metric <- column_sorting_metric
    if (!is.null(column_sorting_algorithm)) args_list$column_sorting_algorithm <- column_sorting_algorithm
    if (!is.null(cycle_start_positions)) args_list$cycle_start_positions <- cycle_start_positions
    if (!is.null(fixed_column)) args_list$fixed_column <- fixed_column
    if (!is.null(random_initializations)) args_list$random_initializations <- random_initializations
    if (!is.null(output_df_path)) args_list$output_df_path <- output_df_path
    if (!is.null(preprocess_data)) args_list$preprocess_data <- preprocess_data
    if (!is.null(default_sorting)) args_list$default_sorting <- default_sorting
    if (!is.null(return_updated_graphing_columns)) args_list$return_updated_graphing_columns <- return_updated_graphing_columns
    if (!is.null(verbose)) args_list$verbose <- verbose
    if (!is.null(print_params)) args_list$print_params <- print_params
    if (!is.null(load_df)) args_list$load_df <- load_df
    if (!is.null(make_intermediate_neighbornet_plots)) args_list$make_intermediate_neighbornet_plots <- make_intermediate_neighbornet_plots

    # Dynamically call function
    result <- do.call(data_sort, args_list)

    if (!quiet) {
        print(result)
    }
}
