#' @export
run_determine_weighted_layer_free_objective_cli <- function(args) {
    show_help <- function() {
        cat("
Usage: alluvialmatch determine_weighted_layer_free_objective --df FILE [options]

Required:
  --df                          CSV input file or RDS with 2–3 columns

Optional:
  --minimum_edge_weight            Minimum edge weight to include (default: 0)
")
        quit(save = "no", status = 0)
    }

    if (length(args) == 0 || any(args %in% c("--help", "-h"))) show_help()

    get_arg <- function(flags, default = NULL) {
        i <- which(args %in% flags)
        if (length(i) > 0 && i < length(args)) args[i + 1] else default
    }

    get_bool_arg <- function(flags, default = FALSE) {
        val <- get_arg(flags)
        if (is.null(val)) return(default)
        tolower(val) %in% c("true", "t", "1")
    }

    get_numeric_arg <- function(flags, default = NULL) {
        val <- get_arg(flags)
        if (is.null(val)) return(default)
        as.numeric(val)
    }

    df             <- get_arg("--df")
    min_edge_weight <- get_numeric_arg("--minimum_edge_weight", 0)

    if (is.null(df)) {
        cat("\nError: --df is required.\n\n")
        show_help()
    }

    result <- determine_weighted_layer_free_objective(
        crossing_edges = df,
        minimum_edge_weight = min_edge_weight
    )

    print(result)
}
