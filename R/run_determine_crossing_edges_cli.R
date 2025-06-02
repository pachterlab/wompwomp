#' @export
run_determine_crossing_edges_cli <- function(args) {
    show_help <- function() {
        cat("
Usage: alluvialmatch determine_crossing_edges --df FILE --column1 C1 --column2 C2 [options]

Required:
  --df                          CSV input file or RDS with 2–3 columns

Optional:
  --column1, -c1                   Name of first column
  --column2, -c2                   Name of second column
  --column_weights                 Comma-separated weights or column name (default: 'value')
  --minimum_edge_weight            Minimum edge weight to include (default: 0)
  --output_df_path                 Path to save resulting edge table (e.g. df.csv)
  --return_weighted_layer_free_objective  TRUE/FALSE (default: FALSE)
  --quiet                      don't show stdout
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
    column1        <- get_arg(c("--column1", "-c1"), NULL)
    column2        <- get_arg(c("--column2", "-c2"), NULL)
    column_weights <- get_arg("--column_weights", "value")
    min_edge_weight <- get_numeric_arg("--minimum_edge_weight", 0)
    output_df_path <- get_arg("--output_df_path")
    return_objective <- get_bool_arg("--return_weighted_layer_free_objective", FALSE)
    quiet <- get_bool_arg("--quiet", FALSE)

    if (is.null(df)) {
        cat("\nError: --df is required.\n\n")
        show_help()
    }

    result <- determine_crossing_edges(
        df = df,
        column1 = column1,
        column2 = column2,
        column_weights = column_weights,
        minimum_edge_weight = min_edge_weight,
        output_df_path = output_df_path,
        return_weighted_layer_free_objective = return_objective
    )

    # Output logic
    if (!quiet) {
        if (return_objective) {
            cat("Weighted layer-free objective:\n")
            print(result$WLF_objective)
        } else {
            cat("Crossing edges - format list[((l1, r1, e1), (l2, r2, e2)), ...]:\n")
            print(result$crossing_edges_df)
        }
    }
}
