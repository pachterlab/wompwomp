#' @noRd
run_determine_weighted_layer_free_objective_cli <- function(args) {
    if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
        cat("
Usage: wompwomp determine_weighted_layer_free_objective --df FILE [options]

Required:
  -i, --input, --df            A CSV path or data frame as outputted with crossing_edges_df (in R) or output_df_path (as a file) from determine_crossing_edges.

Optional:
  -v, --verbose             If TRUE, will display messages during the function.
  --print_params            If TRUE, will print function params.
  -q, --quiet               Don't show stdout
")
        quit(save = "no", status = 0)
    }

    # Required arguments
    df <- get_arg(args, c("-i", "--input", "--df"), required = TRUE)

    # Optional arguments
    verbose <- store_true(args, c("-v", "--verbose"))
    print_params <- store_true(args, c("--print_params"))
    quiet <- store_true(args, c("-q", "--quiet"))

    # Base argument list with required args
    args_list <- list(
        df = df
    )

    # Conditionally add optional args if not NULL
    if (!is.null(verbose)) args_list$verbose <- verbose
    if (!is.null(print_params)) args_list$print_params <- print_params

    # Dynamically call function
    result <- do.call(determine_weighted_layer_free_objective, args_list)

    if (!quiet) {
        print(result)
    }
}
