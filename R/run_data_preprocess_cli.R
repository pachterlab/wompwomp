#' @noRd
run_data_preprocess_cli <- function(args) {
    if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
        cat("
Usage: wompwomp data_preprocess --input INPUT graphing_columns GRAPHING_COLUMNS [options]

Required:
  -i, --input, --df                  A data frame, tibble, or CSV file path. Must be in one of two formats:
(1) column_weights == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two graphing_columns).
(2) column_weights != NULL: Each row represents a combination of groupings, each column from graphing_columns represents a grouping, and the column column_weights represents the number of entities in that combination of groupings. Must contain at least three columns (two graphing_columns, one column_weights}).
  -g, --graphing_columns    Vector of column names from df to be used in graphing (i.e., alluvial plotting).

Optional:
  -w, --column_weights      Column name from df that contains the weights of each combination of groupings if df is in format (2) (see above).
  --default_sorting Character. Default column sorting in data_preprocess if integer columns do not exist. Options are 'alphabetical' (default), 'reverse_alphabetical', 'increasing', 'decreasing', 'random'.
  --set_seed Integer. Random seed for when default_sorting == 'random' parameter.
  -o, --output_df_path      Output path for the output data frame, in CSV format. If NULL, then will not be saved.
  -v, --verbose             If TRUE, will display messages during the function.
  --print_params            If TRUE, will print function params.
  -q, --quiet               Don't show stdout
")
        quit(save = "no", status = 0)
    }

    # Required arguments
    df <- get_arg(args, c("-i", "--input", "--df"), required = TRUE)
    graphing_columns <- get_multi_arg(args, c("-g", "--graphing_columns"), required = TRUE)

    # Optional arguments
    column_weights <- get_arg(args, c("-w", "--column_weights"))
    default_sorting <- get_arg(args, c("--default_sorting"))
    output_df_path <- get_arg(args, c("-o", "--output_df_path"))
    verbose <- store_true(args, c("-v", "--verbose"))
    print_params <- store_true(args, c("--print_params"))
    quiet <- store_true(args, c("-q", "--quiet"))

    # Hidden arguments
    load_df <- store_false(args, c("--disable_load_df"))
    do_gather_set_data <- store_true(args, c("--do_gather_set_data"))
    color_band_column <- get_arg(args, c("--color_band_column"))

    # Base argument list with required args
    args_list <- list(
        df = df,
        graphing_columns = graphing_columns
    )

    # Conditionally add optional args if not NULL
    if (!is.null(column_weights)) args_list$column_weights <- column_weights
    if (!is.null(default_sorting)) args_list$default_sorting <- default_sorting
    if (!is.null(output_df_path)) args_list$output_df_path <- output_df_path
    if (!is.null(verbose)) args_list$verbose <- verbose
    if (!is.null(print_params)) args_list$print_params <- print_params
    if (!is.null(load_df)) args_list$load_df <- load_df
    if (!is.null(do_gather_set_data)) args_list$do_gather_set_data <- do_gather_set_data
    if (!is.null(color_band_column)) args_list$color_band_column <- color_band_column

    # Dynamically call function
    result <- do.call(data_preprocess, args_list)

    if (!quiet) {
        print(result)
    }
}
