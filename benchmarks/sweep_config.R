# Parameter grids for run_benchmark.R sweep mode.
#
# DEFAULT_SWEEP covers the axes that actually move sort_to_uncross()'s cost
# (see generate_data.R docstring). "tsp" is deliberately excluded from the
# default grid since it's an exact exponential-time solver -- add it
# explicitly with a small n_categories cap (<=15) via TSP_SWEEP below.
#
# method %in% c("greedy_wolf", "greedy_wblf") only support n_columns == 2
# (sort_greedy_wolf() stops otherwise); build_grid() in run_benchmark.R
# drops any row combining them with n_columns != 2.

DEFAULT_SWEEP <- list(
    n_rows = c(1000, 10000, 100000),
    n_columns = c(2, 3, 4),
    n_categories = c(4, 8, 16),
    method = c("neighbornet", "greedy_wblf"),
    column_method = c("tsp", "none")
)

TSP_SWEEP <- list(
    n_rows = c(10000),
    n_columns = c(2, 3),
    n_categories = c(4, 8, 12, 15),
    method = c("tsp"),
    column_method = c("none")
)

SMOKE_SWEEP <- list(
    n_rows = c(500),
    n_columns = c(2),
    n_categories = c(4),
    method = c("neighbornet", "greedy_wblf"),
    column_method = c("none")
)

SWEEPS <- list(default = DEFAULT_SWEEP, tsp = TSP_SWEEP, smoke = SMOKE_SWEEP)

# Held fixed across every row in a sweep.
FIXED_PARAMS <- list(
    wt = "value",
    weight_scalar = 5e5,
    seed = 0
)

DEFAULT_TIMEOUT_S <- 600
