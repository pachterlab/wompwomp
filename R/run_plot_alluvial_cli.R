#' @export
run_plot_alluvial_cli <- function(args) {
    if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
        cat("
Usage: alluvialmatch plot_alluvial --input INPUT [options]

Required:
  -i, --input, --df               A data frame, tibble, or CSV file path. Must be in one of two formats:
(1) column_weights == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two graphing_columns).
(2) column_weights != NULL: Each row represents a combination of groupings, each column from graphing_columns represents a grouping, and the column column_weights represents the number of entities in that combination of groupings. Must contain at least three columns (two graphing_columns, one column_weights).
  -o, --output_plot_path          File path to save the plot (e.g., 'plot.png'). If NULL, then will not be saved.

Optional:
  -g, --graphing_columns    Vector of column names from df to be used in graphing (i.e., alluvial plotting). Mutually exclusive with column1 and column2.
  -c1, --column1            Can be used along with column2 in place of graphing_columns if working with two columns only. Mutually exclusive with graphing_columns.
  -c2, --column2            Can be used along with column1 in place of graphing_columns if working with two columns only. Mutually exclusive with graphing_columns.
  -w, --column_weights      Column name from df that contains the weights of each combination of groupings if df is in format (2) (see above).
  -s, --sorting_algorithm       Algorithm with which to sort the values in the dataframe. Can choose from {'neighbornet', 'greedy_WOLF', 'greedy_WBLF', 'None'. 'neighbornet' performs sorting with NeighborNet (Bryant and Moulton, 2004). 'greedy_WOLF' implements a custom greedy algorithm where one layer is fixed, and the other layer is sorted such that each node is positioned as close to its largest parent from the fixed side as possible in a greedy fashion. 'greedy_WBLF' implements the 'greedy_WOLF' algorithm described previously twice, treating each column as fixed in one iteration and free in the other iteration. 'greedy_WOLF' and 'greedy_WBLF' are only valid when graphing_columns has exactly two entries.
  --optimize_column_order   If TRUE, will optimize the order of graphing_columns to minimize edge overlap. Only applies when sorting_algorithm == 'neighbornet' and length(graphing_columns) > 2.
  --optimize_column_order_per_cycle         If TRUE, will optimize the order of graphing_columns to minimize edge overlap upon each cycle. If FALSE, will optimize the order of graphing_columns to minimize edge overlap on the beginning cycle only. Only applies when sorting_algorithm == 'neighbornet' and length(graphing_columns) > 2.
  --fixed_column            Name or position of the column in graphing_columns to keep fixed during sorting. Only applies when sorting_algorithm == 'greedy_WOLF'.
  --random_initializations  Number of random initializations for the positions of each grouping in graphing_columns. Only applies when sorting_algorithm == 'greedy_WOLF' or sorting_algorithm == 'greedy_WBLF'.
  --set_seed                Random seed for the random_initializations parameter. Only applies when sorting_algorithm == 'greedy_WOLF' or sorting_algorithm == 'greedy_WBLF'.
  --disable_color_boxes     Whether to color the strata/boxes (representing groups)
  --color_bands             Whether to color the alluvia/edges (connecting the strata)
  --color_list              List of colors to override default group colors.
  --color_band_list         List of colors to override default band colors.
  --color_band_column       Which column to use for coloring bands
  --color_band_boundary     Whether or not to color boundaries between bands
  --disable_match_colors            Assigns consistent colors between column1 and column2 where matched.
  --alluvial_alpha          Numeric between 0 and 1. Transparency level for the alluvial bands.
  --disable_include_labels_in_boxes Whether to include text labels inside the rectangular group boxes
  --disable_include_axis_titles     Whether to display axis titles for column1 and column2.
  --disable_include_group_sizes     Includes group sizes in the labels
  -o, --output_df_path      Output path for the output data frame, in CSV format. If NULL, then will not be saved.
  --output_df_path          Output path for the output data frame, in CSV format. If NULL, then will not be saved
  --disable_preprocess_data         If TRUE, will preprocess the data with the data_preprocess function.
  --box_width Numeric between 0 and 1. Box width
  --text_width Numeric between 0 and 1. Text width
  --min_text Integer greater than 0. Min text
  --save_height Integer greater than 0. Save height, in inches
  --save_width Integer greater than 0. Save width, in inches
  -v, --verbose             If TRUE, will display messages during the function.
  -q, --quiet               Don't show stdout
")
        quit(save = "no", status = 0)
    }

    # Required arguments
    df <- get_arg(args, c("-i", "--input", "--df"), required = TRUE)
    output_plot_path <- get_arg(args, c("-o", "--output_plot_path"), required = TRUE)

    # Optional arguments
    graphing_columns <- get_multi_arg(args, c("-g", "--graphing_columns"))
    column1 <- get_arg(args, c("-c1", "--column1"))
    column2 <- get_arg(args, c("-c2", "--column2"))
    column_weights <- get_arg(args, c("-w", "--column_weights"))
    sorting_algorithm <- get_arg(args, c("-s", "--sorting_algorithm"))
    optimize_column_order <- store_true(args, c("--optimize_column_order"))
    optimize_column_order_per_cycle <- store_true(args, c("--optimize_column_order_per_cycle"))
    fixed_column <- get_fixed_column(args, "--fixed_column")
    random_initializations <- get_numeric_arg(args, "--random_initializations")
    set_seed <- get_numeric_arg(args, "--set_seed")
    color_boxes <- store_false(args, c("--disable_color_boxes"))
    color_bands <- store_true(args, c("--color_bands"))
    color_list <- get_multi_arg(args, c("--color_list"))
    color_band_list <- get_multi_arg(args, c("--color_band_list"))
    color_band_column <- get_arg(args, c("--color_band_column"))
    color_band_boundary <- store_true(args, c("--color_band_boundary"))
    match_colors <- store_false(args, c("--disable_match_colors"))
    alluvial_alpha <- get_numeric_arg(args, "--alluvial_alpha")
    include_labels_in_boxes <- store_false(args, c("--disable_include_labels_in_boxes"))
    include_axis_titles <- store_false(args, c("--disable_include_axis_titles"))
    include_group_sizes <- store_false(args, c("--disable_include_group_sizes"))
    output_df_path <- get_arg(args, c("--output_df_path"))
    preprocess_data <- store_false(args, c("--disable_preprocess_data"))
    box_width <- store_false(args, c("--box_width"))
    text_width <- store_false(args, c("--text_width"))
    min_text <- store_false(args, c("--min_text"))
    save_height <- store_false(args, c("--save_height"))
    save_width <- store_false(args, c("--save_width"))
    verbose <- store_true(args, c("-v", "--verbose"))
    quiet <- store_true(args, c("-q", "--quiet"))

    # Hidden arguments
    make_intermediate_neighbornet_plots <- store_true(args, c("make_intermediate_neighbornet_plots"))

    # Base argument list with required args
    args_list <- list(
        df = df,
        output_plot_path = output_plot_path
    )

    # Conditionally add optional args if not NULL
    if (!is.null(graphing_columns)) args_list$graphing_columns <- graphing_columns
    if (!is.null(column1)) args_list$column1 <- column1
    if (!is.null(column2)) args_list$column2 <- column2
    if (!is.null(column_weights)) args_list$column_weights <- column_weights
    if (!is.null(sorting_algorithm)) args_list$sorting_algorithm <- sorting_algorithm
    if (!is.null(optimize_column_order)) args_list$optimize_column_order <- optimize_column_order
    if (!is.null(optimize_column_order_per_cycle)) args_list$optimize_column_order_per_cycle <- optimize_column_order_per_cycle
    if (!is.null(fixed_column)) args_list$fixed_column <- fixed_column
    if (!is.null(random_initializations)) args_list$random_initializations <- random_initializations
    if (!is.null(set_seed)) args_list$set_seed <- set_seed
    if (!is.null(color_boxes)) args_list$color_boxes <- color_boxes
    if (!is.null(color_bands)) args_list$color_bands <- color_bands
    if (!is.null(color_list)) args_list$color_list <- color_list
    if (!is.null(color_band_list)) args_list$color_band_list <- color_band_list
    if (!is.null(color_band_column)) args_list$color_band_column <- color_band_column
    if (!is.null(color_band_boundary)) args_list$color_band_boundary <- color_band_boundary
    if (!is.null(match_colors)) args_list$match_colors <- match_colors
    if (!is.null(alluvial_alpha)) args_list$alluvial_alpha <- alluvial_alpha
    if (!is.null(include_labels_in_boxes)) args_list$include_labels_in_boxes <- include_labels_in_boxes
    if (!is.null(include_axis_titles)) args_list$include_axis_titles <- include_axis_titles
    if (!is.null(include_group_sizes)) args_list$include_group_sizes <- include_group_sizes
    if (!is.null(output_df_path)) args_list$output_df_path <- output_df_path
    if (!is.null(preprocess_data)) args_list$preprocess_data <- preprocess_data
    if (!is.null(box_width)) args_list$box_width <- box_width
    if (!is.null(text_width)) args_list$text_width <- text_width
    if (!is.null(min_text)) args_list$min_text <- min_text
    if (!is.null(save_height)) args_list$save_height <- save_height
    if (!is.null(save_width)) args_list$save_width <- save_width
    if (!is.null(verbose)) args_list$verbose <- verbose
    if (!is.null(make_intermediate_neighbornet_plots)) args_list$make_intermediate_neighbornet_plots <- make_intermediate_neighbornet_plots

    # Dynamically call function
    result <- do.call(plot_alluvial, args_list)

    if (!quiet) {
        if (!is.null(output_df_path)) {
            print(sprintf("Data frame saved to to=%s", output_df_path))
        }
        if (!is.null(output_plot_path)) {
            print(sprintf("Plot saved to to=%s", output_plot_path))
        }
    }
}
