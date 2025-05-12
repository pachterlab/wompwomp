#' @export
run_cli <- function(args) {
  if (length(args) == 0 || (length(args) == 1 && args[1] %in% c("--help", "-h"))) {
    cat("
Usage:
  alluvialmatch <command> [options]

Commands:
  plot_alluvial     Generate an alluvial plot
  greedy_wolf       Run the greedy heuristic
  determine_crossing_edges    Determine sum of products of overlapping edges
  determine_weighted_layer_free_objective    Determine sum of products of overlapping edges

Use:
  alluvialmatch plot_alluvial --help
  alluvialmatch greedy_wolf --help
  alluvialmatch determine_crossing_edges --help

")
    quit(save = "no", status = 0)
  }

  command <- args[[1]]
  sub_args <- args[-1]

  if (command == "plot_alluvial") {
    run_plot_alluvial_cli(sub_args)
  } else if (command == "greedy_wolf") {
    run_greedy_wolf_cli(sub_args)
  } else if (command == "determine_crossing_edges") {
      run_determine_crossing_edges_cli(sub_args)
  } else if (command == "determine_weighted_layer_free_objective") {
      run_determine_weighted_layer_free_objective_cli(sub_args)
  } else {
    cat("Unknown command: ", command, "\n")
    run_cli("--help")
  }
}
