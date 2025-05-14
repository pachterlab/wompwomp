#' @export
run_greedy_wolf_cli <- function(args) {
  if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
    cat("
Usage: alluvialmatch greedy_wolf --df FILE [options]

Required:
  --df         CSV input file

Optional:
  --column1       Name of first column
  --column2       Name of second column
  --column_weights  Comma-separated numeric values (e.g. 1,1.5,2)
  --fixed_column                   Fix one column for one-layer free layout (1, 2, or column name)
  --random_initializations         Number of random WLF initializations to run (default: 1)
  --output_df_path                 Path to save resulting edge table (e.g. df.csv)
  --sorting_algorithm           greedy_WBLF, greedy_WOLF, or None
  --set_seed                    seed for random initializations
")
    quit(save = "no", status = 0)
  }

  get_arg <- function(flag, default = NULL) {
    i <- which(args == flag)
    if (length(i) > 0 && i < length(args)) args[i + 1] else default
  }

  get_numeric_arg <- function(flags, default = NULL) {
      val <- get_arg(flags)
      if (is.null(val)) return(default)
      as.numeric(val)
  }

  get_numeric_list_arg <- function(flag) {
    val <- get_arg(flag)
    if (is.null(val)) return(NULL)
    as.numeric(strsplit(val, ",")[[1]])
  }

  df              <- get_arg("--df")
  column1         <- get_arg("--column1")
  column2         <- get_arg("--column2")
  column_weights  <- get_numeric_list_arg("--column_weights")
  fixed_column    <- get_arg("--fixed_column")
  random_initializations <- get_numeric_arg("--random_initializations", 1)
  output_df_path <- get_arg("--output_df_path")
  sorting_algorithm       <- get_arg("--sorting_algorithm")
  set_seed       <- get_arg("--set_seed")

  if (is.null(df)) {
    cat("Error: --df is required.\n\n")
    run_greedy_wolf_cli("--help")
  }

  result <- greedy_wolf(
    df = df,
    column1 = column1,
    column2 = column2,
    column_weights = column_weights,
    fixed_column = fixed_column,
    random_initializations = random_initializations,
    output_df_path = output_df_path,
    sorting_algorithm=sorting_algorithm,
    set_seed=set_seed
  )

  print(result)
}
