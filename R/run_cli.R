#' @export
run_cli <- function(args) {
  if (length(args) == 0 || (length(args) == 1 && args[1] %in% c("--help", "-h"))) {
    cat("
Usage:
  alluvialmatch <command> [options]

Commands:
  plot_alluvial     Generate an Alluvial Plot with Minimal Cluster Cross-over (runs data_preprocess and data_sort internally)
  data_preprocess   Preprocess data
  data_sort       Sorts a dataframe (runs data_preprocess internally)
  determine_crossing_edges    Determine overlapping edges
  determine_weighted_layer_free_objective    Compute crossing objective

Use:
  alluvialmatch COMMAND --help
")
    quit(save = "no", status = 0)
  }

  command <- args[[1]]
  sub_args <- args[-1]

  if (command == "plot_alluvial") {
    run_plot_alluvial_cli(sub_args)
  } else if (command == "data_preprocess") {
    run_data_preprocess_cli(sub_args)
  } else if (command == "data_sort") {
    run_data_sort_cli(sub_args)
  } else if (command == "determine_crossing_edges") {
      run_determine_crossing_edges_cli(sub_args)
  } else if (command == "determine_weighted_layer_free_objective") {
      run_determine_weighted_layer_free_objective_cli(sub_args)
  } else {
    cat("Unknown command: ", command, "\n")
    run_cli("--help")
  }
}
