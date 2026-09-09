# Parameter grids for run_benchmark.R sweep mode.
#
# DEFAULT_SWEEP covers the axes that actually move sort_to_uncross()'s cost
# (see generate_data.R docstring). "tsp" is deliberately excluded from the
# default grid since it's an exact exponential-time solver -- add it
# explicitly with a small n_categories cap (<=15) via TSP_SWEEP below.
#
# method %in% c("greedy_wolf", "greedy_wblf", "barycenter", "median",
# "barycenter_one_sided", "median_one_sided") only support n_columns == 2
# (sort_greedy_wolf()/sort_barycenter_median() stop otherwise); build_grid()
# in run_benchmark.R drops any row combining them with n_columns != 2.

DEFAULT_SWEEP <- list(
    n_rows = c(1000, 10000, 100000),
    n_columns = c(2, 3, 4),
    n_categories = c(4, 8, 16),
    method = c("neighbornet", "greedy_wblf", "barycenter", "median", "barycenter_one_sided", "median_one_sided"),
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
    method = c("neighbornet", "greedy_wblf", "barycenter", "median", "barycenter_one_sided", "median_one_sided"),
    column_method = c("none")
)

# For a fixed n_categories, the two-layer graph is "complete" (every
# category-pair combination occurs at least once) once n_rows gets large
# enough relative to n_categories^2 -- and for a complete bipartite graph,
# the *unweighted* crossing count is a topological invariant (always
# C(n_categories, 2)^2), identical for every method regardless of ordering.
# DEFAULT_SWEEP's n_rows/n_categories combinations mostly land in that
# complete regime, so the unweighted objective can't differentiate methods
# there at all. SPARSE_SWEEP deliberately keeps n_rows small relative to
# n_categories^2 across the board (worst case here is 5000 rows against
# 32^2 = 1024 possible combos, most are far sparser) so the induced graph
# stays incomplete and unweighted crossing count has real room to vary by
# method. n_columns is fixed at 2 -- this is about the objective metric
# itself, not scaling, so the extra n_columns/column_method axes from
# DEFAULT_SWEEP aren't needed here.
SPARSE_SWEEP <- list(
    n_rows = c(50, 200, 1000, 5000),
    n_columns = c(2),
    n_categories = c(8, 16, 32),
    method = c("neighbornet", "greedy_wblf", "barycenter", "median", "barycenter_one_sided", "median_one_sided"),
    column_method = c("none")
)

SWEEPS <- list(default = DEFAULT_SWEEP, tsp = TSP_SWEEP, smoke = SMOKE_SWEEP, sparse = SPARSE_SWEEP)

# Held fixed across every row in a sweep.
FIXED_PARAMS <- list(
    wt = "value",
    weight_scalar = 5e5,
    seed = 0
)

DEFAULT_TIMEOUT_S <- 600
