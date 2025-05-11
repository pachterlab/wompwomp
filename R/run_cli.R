#' @export
run_cli <- function(args) {
  if (length(args) == 0 || (length(args) == 1 && args[1] %in% c("--help", "-h"))) {
    cat("
Usage:
  alluvialmatch <command> [options]

Commands:
  plot_alluvial     Generate an alluvial plot
  greedy_wolf       Run the greedy heuristic

Use:
  alluvialmatch plot_alluvial --help
  alluvialmatch greedy_wolf --help

")
    quit(save = "no", status = 0)
  }

  command <- args[[1]]
  sub_args <- args[-1]

  if (command == "plot_alluvial") {
    run_plot_alluvial_cli(sub_args)
  } else if (command == "greedy_wolf") {
    run_greedy_wolf_cli(sub_args)
  } else {
    cat("Unknown command: ", command, "\n")
    run_cli("--help")
  }
}
