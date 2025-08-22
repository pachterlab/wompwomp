#' @noRd
run_plot_alluvial_internal_cli <- function(args) {
    if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
        cat("
Usage: wompwomp plot_alluvial_internal --input INPUT --column_weights COLUMN_WEIGHT --graphing columns COLUMN1 COLUMN2 ... [options]

Required:
  -i, --input, --df               A data frame, tibble, or CSV file path. Must be in the format as the output of wompwomp::data_sort.
  -g, --graphing_columns          Vector of column names from df to be used in graphing (i.e., alluvial plotting). Mutually exclusive with column1 and column2.
  -w, --column_weights            Column name from df that contains the weights of each combination of groupings if df is in format (2) (see above).
  -o, --output_plot_path          File path to save the plot (e.g., 'plot.png'). If NULL, then will not be saved.

Optional:
  -s, --sorting_algorithm       Algorithm with which to sort the values in the dataframe. Can choose from {'neighbornet', 'tsp', 'greedy_wolf', 'greedy_wblf', 'random', 'none'. 'neighbornet' performs sorting with NeighborNet (Bryant and Moulton, 2004). 'tsp' performs Traveling Salesman Problem solver from the TSP package. 'greedy_wolf' implements a custom greedy algorithm where one layer is fixed, and the other layer is sorted such that each node is positioned as close to its largest parent from the fixed side as possible in a greedy fashion. 'greedy_wblf' implements the 'greedy_wolf' algorithm described previously twice, treating each column as fixed in one iteration and free in the other iteration. 'greedy_wolf' and 'greedy_wblf' are only valid when graphing_columns has exactly two entries. 'random' randomly maps blocks. 'none' keeps the mappings as-is when passed into the function.
  --optimize_column_order   If TRUE, will optimize the order of graphing_columns to minimize edge overlap. Only applies when sorting_algorithm == 'neighbornet' or 'tsp' and length(graphing_columns) > 2.
  --optimize_column_order_per_cycle         If TRUE, will optimize the order of graphing_columns to minimize edge overlap upon each cycle. If FALSE, will optimize the order of graphing_columns to minimize edge overlap on the beginning cycle only. Only applies when sorting_algorithm == 'neighbornet' or 'tsp' and length(graphing_columns) > 2.
  --matrix_initialization_value Positive integer. Initialized value in distance matrix for nodes in different layers without a shared edge/path. Only applies when sorting_algorithm == 'neighbornet' or 'tsp'.
  --same_side_matrix_initialization_value Positive integer. Initialized value in distance matrix for nodes in the same layer. Only applies when sorting_algorithm == 'neighbornet' or 'tsp'.
  --weight_scalar Positive integer. Scalar with which to multiply edge weights after taking their -log in the distance matrix for nodes with a nonzero edge. Only applies when sorting_algorithm == 'neighbornet' or 'tsp'.
  --matrix_initialization_value_column_order Positive integer. Initialized value in distance matrix for optimizing column order. Only applies when sorting_algorithm == 'neighbornet' or 'tsp' and optimize_column_order is TRUE.
  --weight_scalar_column_order Positive integer. Scalar with which to loss function after taking their log1p in the distance matrix for optimizing column order. Only applies when sorting_algorithm == 'neighbornet' or 'tsp' and optimize_column_order is TRUE.
  --column_sorting_metric Character. Metric to use for determining column order. Options are 'edge_crossing' (default) or 'ARI'. Only applies when sorting_algorithm == 'neighbornet' or 'tsp' and optimize_column_order is TRUE.
  --column_sorting_algorithm Character. Algorithm to use for determining column order. Options are 'tsp' (default) or 'neighbornet'. Only applies when sorting_algorithm == 'neighbornet' or 'tsp' and optimize_column_order is TRUE.
  --cycle_start_positions Set. Cycle start positions to consider. Anything outside this set will be skipped. Only applies when sorting_algorithm == 'neighbornet' or 'tsp'.
  --fixed_column            Name or position of the column in graphing_columns to keep fixed during sorting. Only applies when sorting_algorithm == 'greedy_wolf'.
  --random_initializations  Number of random initializations for the positions of each grouping in graphing_columns. Only applies when sorting_algorithm == 'greedy_wolf' or sorting_algorithm == 'greedy_wblf'.
  --set_seed                Random seed for the random_initializations parameter (only applies when sorting_algorithm == 'greedy_wolf' or sorting_algorithm == 'greedy_wblf'), random initialization/shuffling of blocks (only applies when default_sorting == 'random' or sorting_algorithm == 'random'), TSP solver for block order or optimizing column order (only applies when sorting_algorithm == 'tsp' or column_sorting_algorithm == 'tsp'), and louvain/leiden clustering (only applies when coloring_algorithm == 'advanced').
  --disable_color_boxes     Whether to color the strata/boxes (representing groups)
  --color_bands             Whether to color the alluvia/edges (connecting the strata)
  --color_list              List of colors to override default group colors.
  --color_band_list         List of colors to override default band colors.
  --color_band_column       Which column to use for coloring bands
  --color_val Optional named list where the entries are colors and the names correspond to values of the dataframe that should use those colors
  --color_band_boundary     Whether or not to color boundaries between bands
  --coloring_algorithm Character. Matching colors methods. Choices are 'advanced' (default), 'none', 'left', 'right', or any value in graphing_columns.
  --coloring_algorithm_advanced_option Character. If coloring_algorithm == 'advanced', then choose graph clustering algorithm. Choices are 'leiden' (default) or 'louvain'.
  --resolution Numeric If coloring_algorithm == 'advanced', then choose resolution for the graph clustering algorithm. Affects coloring of both bands and boxes.
  --cutoff Numeric If coloring_algorithm != 'none' and coloring_algorithm != 'advanced', sets the cutoff for color matching, below which a new color will be assigned.
  --alluvial_alpha          Numeric between 0 and 1. Transparency level for the alluvial bands.
  --disable_include_labels_in_boxes Whether to include text labels inside the rectangular group boxes
  --disable_include_axis_titles     Whether to display axis titles for column1 and column2.
  --include_group_sizes     Includes group sizes in the labels
  -o, --output_df_path      Output path for the output data frame, in CSV format. If NULL, then will not be saved.
  --output_df_path          Output path for the output data frame, in CSV format. If NULL, then will not be saved
  --disable_preprocess_data         If TRUE, will preprocess the data with the data_preprocess function.
  --default_sorting Character. Default column sorting in data_preprocess if integer columns do not exist. Will not affect output if sorting_algorithm == 'neighbornet' or 'tsp'. Options are 'alphabetical' (default), 'reverse_alphabetical', 'increasing', 'decreasing', 'random'.
  --box_width Numeric between 0 and 1. Box width
  --text_width Numeric between 0 and 1. Text width
  --min_text Integer greater than 0. Min text
  --text_size Integer greater than 0. Text size (works whether auto_adjust_text is TRUE or FALSE).
  --auto_adjust_text Whether to automatically adjust text size to fit in box.
  --axis_text_size Integer greater than 0. Axis text size
  --axis_text_vjust Integer. Axis text vjust
  --save_height Integer greater than 0. Save height, in inches
  --save_width Integer greater than 0. Save width, in inches
  --dpi Integer greater than 0. DPI for output_plot_path, if output_plot_path is a raster image or rasterise_alluvia is TRUE
  --rasterise_alluvia Logical. Whether to rasterize the alluvia ifoutput_plot_path is a PDF. Can save space if DPI low enough
  --keep_y_labels Keep y labels
  --box_line_width Box line width
  --environment  Python environment (if applicable). Default: 'wompwomp_env'
  --disable_use_conda  Whether or not to use conda for Python (if applicable)
  -v, --verbose             If TRUE, will display messages during the function.
  --print_params            If TRUE, will print function params.
  -q, --quiet               Don't show stdout
")
        quit(save = "no", status = 0)
    }

    # Required arguments
    df <- get_arg(args, c("-i", "--input", "--df"), required = TRUE)
    graphing_columns <- get_multi_arg(args, c("-g", "--graphing_columns"), required = TRUE)
    column_weights <- get_arg(args, c("-w", "--column_weights"), required = TRUE)
    output_plot_path <- get_arg(args, c("-o", "--output_plot_path"), required = TRUE)

    # Optional arguments
    color_boxes <- store_false(args, c("--disable_color_boxes"))
    color_bands <- store_true(args, c("--color_bands"))
    color_list <- get_multi_arg(args, c("--color_list"))
    color_band_list <- get_multi_arg(args, c("--color_band_list"))
    color_band_column <- get_arg(args, c("--color_band_column"))
    color_val <- get_arg(args, c("--color_val"))
    color_band_boundary <- store_true(args, c("--color_band_boundary"))
    coloring_algorithm <- get_arg(args, c("--coloring_algorithm"))
    coloring_algorithm_advanced_option <- get_arg(args, c("--coloring_algorithm_advanced_option"))
    resolution <- get_numeric_arg(args, c("--resolution"))
    cutoff <- get_numeric_arg(args, c("--cutoff"))
    alluvial_alpha <- get_numeric_arg(args, "--alluvial_alpha")
    include_labels_in_boxes <- store_false(args, c("--disable_include_labels_in_boxes"))
    include_axis_titles <- store_false(args, c("--disable_include_axis_titles"))
    include_group_sizes <- store_true(args, c("--include_group_sizes"))
    box_width <- get_numeric_arg(args, c("--box_width"))
    text_width <- get_numeric_arg(args, c("--text_width"))
    min_text <- get_numeric_arg(args, c("--min_text"))
    text_size <- get_numeric_arg(args, c("--text_size"))
    auto_adjust_text <- store_false(args, c("--disable_auto_adjust_text"))
    axis_text_size <- get_numeric_arg(args, c("--axis_text_size"))
    axis_text_vjust <- get_numeric_arg(args, c("--axis_text_vjust"))
    save_height <- get_numeric_arg(args, c("--save_height"))
    save_width <- get_numeric_arg(args, c("--save_width"))
    keep_y_labels <- store_true(args, c("--keep_y_labels"))
    dpi <- get_integer_arg(args, c("--dpi"))
    rasterise_alluvia <- store_true(args, c("--rasterise_alluvia"))
    box_line_width <- get_numeric_arg(args, c("--box_line_width"))
    environment <- get_arg(args, c("--environment"), default = "wompwomp_env")
    use_conda <- store_false(args, c("--disable_use_conda"))
    verbose <- store_true(args, c("-v", "--verbose"))
    print_params <- store_true(args, c("--print_params"))
    quiet <- store_true(args, c("-q", "--quiet"))

    # Base argument list with required args
    args_list <- list(
        df = df,
        graphing_columns = graphing_columns,
        column_weights = column_weights,
        output_plot_path = output_plot_path
    )

    # Conditionally add optional args if not NULL
    if (!is.null(color_boxes)) args_list$color_boxes <- color_boxes
    if (!is.null(color_bands)) args_list$color_bands <- color_bands
    if (!is.null(color_list)) args_list$color_list <- color_list
    if (!is.null(color_band_list)) args_list$color_band_list <- color_band_list
    if (!is.null(color_band_column)) args_list$color_band_column <- color_band_column
    if (!is.null(color_val)) args_list$color_val <- color_val
    if (!is.null(color_band_boundary)) args_list$color_band_boundary <- color_band_boundary
    if (!is.null(coloring_algorithm)) args_list$coloring_algorithm <- coloring_algorithm
    if (!is.null(coloring_algorithm_advanced_option)) args_list$coloring_algorithm_advanced_option <- coloring_algorithm_advanced_option
    if (!is.null(resolution)) args_list$resolution <- resolution
    if (!is.null(cutoff)) args_list$cutoff <- cutoff
    if (!is.null(alluvial_alpha)) args_list$alluvial_alpha <- alluvial_alpha
    if (!is.null(include_labels_in_boxes)) args_list$include_labels_in_boxes <- include_labels_in_boxes
    if (!is.null(include_axis_titles)) args_list$include_axis_titles <- include_axis_titles
    if (!is.null(include_group_sizes)) args_list$include_group_sizes <- include_group_sizes
    if (!is.null(box_width)) args_list$box_width <- box_width
    if (!is.null(text_width)) args_list$text_width <- text_width
    if (!is.null(min_text)) args_list$min_text <- min_text
    if (!is.null(text_size)) args_list$text_size <- text_size
    if (!is.null(auto_adjust_text)) args_list$auto_adjust_text <- auto_adjust_text
    if (!is.null(axis_text_size)) args_list$axis_text_size <- axis_text_size
    if (!is.null(axis_text_vjust)) args_list$axis_text_vjust <- axis_text_vjust
    if (!is.null(save_height)) args_list$save_height <- save_height
    if (!is.null(save_width)) args_list$save_width <- save_width
    if (!is.null(dpi)) args_list$dpi <- dpi
    if (!is.null(rasterise_alluvia)) args_list$rasterise_alluvia <- rasterise_alluvia
    if (!is.null(keep_y_labels)) args_list$keep_y_labels <- keep_y_labels
    if (!is.null(box_line_width)) args_list$box_line_width <- box_line_width
    if (!is.null(verbose)) args_list$verbose <- verbose
    if (!is.null(print_params)) args_list$print_params <- print_params
    if (!is.null(environment)) args_list$environment <- environment
    if (!is.null(use_conda)) args_list$use_conda <- use_conda

    # Dynamically call function
    result <- do.call(plot_alluvial_internal, args_list)

    if (!quiet) {
        if (!is.null(output_plot_path)) {
            print(sprintf("Plot saved to to=%s", output_plot_path))
        }
    }
}
