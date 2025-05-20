#' @export
run_plot_alluvial_cli <- function(args) {
  if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
    cat("
Usage: alluvialmatch plot_alluvial --df FILE [options]

Required:
  --df                         CSV input file

Optional:
  --column1                       Name of first column
  --column2                       Name of second column
  --fixed_column                   Fix one column for one-layer free layout (1, 2, or column name)
  --random_initializations         Number of random WLF initializations to run (default: 1)
  --show_group_2_box_labels_in_ascending  TRUE/FALSE (default: FALSE)
  --color_boxes                  TRUE/FALSE (default: TRUE)
  --color_bands                  TRUE/FALSE (default: TRUE)
  --match_colors                 TRUE/FALSE (default: TRUE)
  --alluvial_alpha               Float (default: 0.5)
  --include_labels_in_boxes      TRUE/FALSE (default: TRUE)
  --include_axis_titles          TRUE/FALSE (default: TRUE)
  --include_group_sizes          TRUE/FALSE (default: TRUE)
  --column_weights               Comma-separated numbers (e.g. 1,2)
  --output_plot_path                  Path to save output
  --output_df_path               Path to save resulting edge table (e.g. df.csv)
  --color_list                   Comma-separated hex codes for nodes
  --color_band_column                color column for edges
  --color_band_boundary                color column for boundaries
  --sorting_algorithm           greedy_WBLF, greedy_WOLF, or None
  --color_band_list             Comma-separated hex codes for edges
  --set_seed                    seed for random initializations
  --quiet                      don't show stdout
")
    quit(save = "no", status = 0)
  }

  get_arg <- function(flag, default = NULL) {
    i <- which(args == flag)
    if (length(i) > 0 && i < length(args)) args[i + 1] else default
  }

  get_bool_arg <- function(flag, default = FALSE) {
    val <- get_arg(flag)
    if (is.null(val)) return(default)
    tolower(val) %in% c("true", "t", "1")
  }

  get_list_arg <- function(flag) {
    val <- get_arg(flag)
    if (is.null(val)) return(NULL)
    strsplit(val, ",")[[1]]
  }

  get_numeric_arg <- function(flags, default = NULL) {
      val <- get_arg(flags)
      if (is.null(val)) return(default)
      as.numeric(val)
  }

  get_numeric_list_arg <- function(flag) {
    val <- get_list_arg(flag)
    if (is.null(val)) return(NULL)
    as.numeric(val)
  }

  df          <- get_arg("--df")

  if (is.null(df)) {
    cat("Error: --df is required.\n\n")
    run_plot_alluvial_cli("--help")
  }

  # Optional args
  column1     <- get_arg("--column1")
  column2     <- get_arg("--column2")
  fixed_column    <- get_arg("--fixed_column", 1)
  random_initializations <- get_numeric_arg("--random_initializations", 1)
  color_boxes        <- get_bool_arg("--color_boxes", TRUE)
  color_bands        <- get_bool_arg("--color_bands", FALSE)
  match_colors       <- get_bool_arg("--match_colors", TRUE)
  alluvial_alpha     <- as.numeric(get_arg("--alluvial_alpha", 0.5))
  include_labels_in_boxes <- get_bool_arg("--include_labels_in_boxes", TRUE)
  include_axis_titles     <- get_bool_arg("--include_axis_titles", TRUE)
  include_group_sizes     <- get_bool_arg("--include_group_sizes", TRUE)
  column_weights    <- get_numeric_list_arg("--column_weights")
  output_plot_path       <- get_arg("--output_plot_path", "./alluvial.png")
  output_df_path       <- get_arg("--output_df_path")
  color_list        <- get_list_arg("--color_list")
  color_band_column       <- get_arg("--color_band_column")
  color_band_boundary       <- get_bool_arg("--color_band_boundary", FALSE)
  sorting_algorithm       <- get_arg("--sorting_algorithm", "greedy_WBLF")
  color_band_list       <- get_arg("--color_band_list")
  set_seed       <- get_arg("--set_seed")
  quiet <- get_bool_arg("--quiet", FALSE)

  result <- plot_alluvial(
    df = df,
    column1 = column1,
    column2 = column2,
    fixed_column = fixed_column,
    random_initializations = random_initializations,
    color_boxes = color_boxes,
    color_bands = color_bands,
    match_colors = match_colors,
    alluvial_alpha = alluvial_alpha,
    include_labels_in_boxes = include_labels_in_boxes,
    include_axis_titles = include_axis_titles,
    include_group_sizes = include_group_sizes,
    column_weights = column_weights,
    output_plot_path = output_plot_path,
    output_df_path = output_df_path,
    color_list = color_list,
    color_band_column=color_band_column,
    color_band_boundary=color_band_boundary,
    sorting_algorithm=sorting_algorithm,
    color_band_list=color_band_list,
    set_seed=set_seed
  )


  if (!quiet && !is.null(output_df_path)) {
      message("Dataframe saved to ", output_df_path)
  }

  if (!quiet && !is.null(output_plot_path)) {
    message("Plot saved to ", output_plot_path)
  }
}
