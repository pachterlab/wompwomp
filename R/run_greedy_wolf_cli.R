#' @export
run_greedy_wolf_cli <- function(args) {
  if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
    cat("
Usage: alluvialmatch greedy_wolf --input FILE --column1 C1 --column2 C2 [options]

Required:
  --input         CSV input file
  --column1       Name of first column
  --column2       Name of second column

Optional:
  --column_weights  Comma-separated numeric values (e.g. 1,1.5,2)
")
    quit(save = "no", status = 0)
  }

  get_arg <- function(flag, default = NULL) {
    i <- which(args == flag)
    if (length(i) > 0 && i < length(args)) args[i + 1] else default
  }

  get_numeric_list_arg <- function(flag) {
    val <- get_arg(flag)
    if (is.null(val)) return(NULL)
    as.numeric(strsplit(val, ",")[[1]])
  }

  input_file      <- get_arg("--input")
  column1         <- get_arg("--column1")
  column2         <- get_arg("--column2")
  column_weights  <- get_numeric_list_arg("--column_weights")

  if (is.null(input_file) || is.null(column1) || is.null(column2)) {
    cat("Error: --input, --column1, and --column2 are required.\n\n")
    run_greedy_wolf_cli("--help")
  }

  result <- greedy_wolf(
    df = input_file,
    column1 = column1,
    column2 = column2,
    column_weights = column_weights
  )

  print(result)
}
