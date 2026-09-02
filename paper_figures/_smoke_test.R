#!/usr/bin/env Rscript
# Minimal render smoke test for wompwomp::plot_alluvial(). Exercises the code
# paths the paper figures rely on and fails loudly if any of them error.
# Run inside the Docker image (needs ggfittext / a graphics device).

suppressPackageStartupMessages(library(wompwomp))
set.seed(1)

out <- file.path(tempdir(), "smoke")
dir.create(out, showWarnings = FALSE)

# 3-layer toy dataset, raw per-row form
n <- 400
df <- data.frame(
  a = sample(c("a1", "a2", "a3"), n, TRUE),
  b = sample(c("b1", "b2", "b3", "b4"), n, TRUE),
  c = sample(c("c1", "c2", "c3"), n, TRUE),
  stringsAsFactors = FALSE
)

cases <- list(
  list(sorting_algorithm = "neighbornet", coloring_algorithm = "advanced"),
  list(sorting_algorithm = "neighbornet", coloring_algorithm = "left"),
  list(sorting_algorithm = "neighbornet", coloring_algorithm = "right"),
  list(sorting_algorithm = "neighbornet", coloring_algorithm = "none"),
  list(sorting_algorithm = "neighbornet", coloring_algorithm = "a"),        # named column
  list(sorting_algorithm = "tsp",         coloring_algorithm = "advanced"),
  list(sorting_algorithm = "none",        coloring_algorithm = "left", color_bands = TRUE, color_band_column = "a"),
  list(sorting_algorithm = "neighbornet", coloring_algorithm = "left", flip_xy = TRUE),
  list(sorting_algorithm = "neighbornet", coloring_algorithm = "left", rasterise_alluvia = TRUE),
  list(sorting_algorithm = "greedy_wolf", coloring_algorithm = "advanced", graphing_columns = c("a", "b"))
)

ok <- 0L
for (i in seq_along(cases)) {
  args <- cases[[i]]
  gcols <- if (!is.null(args$graphing_columns)) args$graphing_columns else c("a", "b", "c")
  args$graphing_columns <- NULL
  path <- file.path(out, sprintf("case_%02d.pdf", i))
  p <- do.call(plot_alluvial, c(
    list(df = df, graphing_columns = gcols, optimize_column_order = FALSE,
         output_plot_path = path, verbose = FALSE),
    args
  ))
  stopifnot(inherits(p, "ggplot"), file.exists(path), file.info(path)$size > 0)
  ok <- ok + 1L
  cat(sprintf("  case %02d ok  (%s)\n", i,
              paste(names(args), unlist(args), sep = "=", collapse = ", ")))
}

# 2-layer form-2 (pre-grouped) input
g <- as.data.frame(table(df$a, df$b))
names(g) <- c("a", "b", "value")
p <- plot_alluvial(g, graphing_columns = c("a", "b"), column_weights = "value",
                   sorting_algorithm = "neighbornet",
                   output_plot_path = file.path(out, "case_grouped.pdf"))
stopifnot(inherits(p, "ggplot"))
ok <- ok + 1L
cat("  grouped-input case ok\n")

cat(sprintf("\nplot_alluvial() smoke test: %d/%d cases passed\n", ok, length(cases) + 1L))
