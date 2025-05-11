#' @export
run_plot_alluvial_cli <- function(args) {
  if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
    cat("
Usage: alluvialmatch plot_alluvial --input FILE --column1 C1 --column2 C2 [options]

Required:
  --input                         CSV input file
  --column1                       Name of first column
  --column2                       Name of second column

Optional:
  --show_group_2_box_labels_in_ascending  TRUE/FALSE (default: FALSE)
  --color_boxes                  TRUE/FALSE (default: TRUE)
  --color_bands                  TRUE/FALSE (default: TRUE)
  --match_colors                 TRUE/FALSE (default: TRUE)
  --alluvial_alpha               Float (default: 0.5)
  --include_labels_in_boxes      TRUE/FALSE (default: TRUE)
  --include_axis_titles          TRUE/FALSE (default: TRUE)
  --include_group_sizes          TRUE/FALSE (default: TRUE)
  --column_weights               Comma-separated numbers (e.g. 1,2)
  --output_path                  Path to save output
  --color_list                   Comma-separated hex codes
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

  get_numeric_list_arg <- function(flag) {
    val <- get_list_arg(flag)
    if (is.null(val)) return(NULL)
    as.numeric(val)
  }

  input_file  <- get_arg("--input")
  column1     <- get_arg("--column1")
  column2     <- get_arg("--column2")

  if (is.null(input_file) || is.null(column1) || is.null(column2)) {
    cat("Error: --input, --column1, and --column2 are required.\n\n")
    run_plot_alluvial_cli("--help")
  }

  # Optional args
  show_group_2_box_labels_in_ascending <- get_bool_arg("--show_group_2_box_labels_in_ascending", FALSE)
  color_boxes        <- get_bool_arg("--color_boxes", TRUE)
  color_bands        <- get_bool_arg("--color_bands", TRUE)
  match_colors       <- get_bool_arg("--match_colors", TRUE)
  alluvial_alpha     <- as.numeric(get_arg("--alluvial_alpha", 0.5))
  include_labels_in_boxes <- get_bool_arg("--include_labels_in_boxes", TRUE)
  include_axis_titles     <- get_bool_arg("--include_axis_titles", TRUE)
  include_group_sizes     <- get_bool_arg("--include_group_sizes", TRUE)
  column_weights    <- get_numeric_list_arg("--column_weights")
  output_path       <- get_arg("--output_path")
  color_list        <- get_list_arg("--color_list")

  result <- plot_alluvial(
    df = input_file,
    column1 = column1,
    column2 = column2,
    show_group_2_box_labels_in_ascending = show_group_2_box_labels_in_ascending,
    color_boxes = color_boxes,
    color_bands = color_bands,
    match_colors = match_colors,
    alluvial_alpha = alluvial_alpha,
    include_labels_in_boxes = include_labels_in_boxes,
    include_axis_titles = include_axis_titles,
    include_group_sizes = include_group_sizes,
    column_weights = column_weights,
    output_path = output_path,
    color_list = color_list,
  )

  if (!is.null(output_path)) {
    message("Output saved to ", output_path)
  }
}
