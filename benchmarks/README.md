# wompwomp benchmarks

Timing/resource harness for `wompwomp::sort_to_uncross`, meant to run on
atlas (48 cores / 377G RAM, no job scheduler -- launch long sweeps under
`tmux` or `nohup`). Mirrors the structure of `wompywompy`'s
`benchmarks/`, adapted to R (fork-based process isolation instead of
`multiprocessing`).

## Setup

```
ssh atlas
conda activate wompwomp_env
cd ~/wompwomp
Rscript -e 'devtools::install(".")'   # only needed once / after R code changes
cd benchmarks
```

## Resource etiquette

atlas is a shared, unscheduled box -- other people's jobs run there too.
`run_benchmark.R` caps BLAS/OpenMP thread pools to 1 (`OMP_NUM_THREADS`,
`OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `OMP_THREAD_LIMIT`) before any
package that might link against them is loaded, so each config runs
single-threaded by default and won't grab cores away from anyone else's
running work. This also makes wall-clock numbers comparable across runs
instead of varying with how busy the box is.

To measure realistic multi-threaded throughput instead (only do this if
you've checked `htop`/`free -h` and the box has headroom), export a higher
value before invoking the script -- the harness only sets a var if it isn't
already exported:
```
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
Rscript run_benchmark.R sweep --which default --out results/default_mt.csv
```

## What's being measured

Each config runs `sort_to_uncross()` end-to-end (synthetic data generation +
aggregation + sorting/column-ordering, no coloring or plotting) in its own
forked child process, isolated so a hung or crashing config can't take down
the whole sweep. Recorded per run: wall time, user/sys CPU time, peak RSS
(read from `/proc/self/status`), and the number of unique alluvia (rows in
the aggregated table) actually sorted.

Before the sweep starts, `run_rows()` runs each `method` once on a tiny
dataset in the parent process (`warm_up_namespaces()` in `run_benchmark.R`)
so every fork inherits already-loaded package namespaces via copy-on-write.
Without this, R lazy-loads `wompwomp`/`dplyr`/`igraph`/`TSP` fresh in every
single fork -- profiling showed that costing ~0.15-0.2s per config, which
swamps the actual algorithm time for small/fast configs.

Isolation is done with `parallel::mcparallel()`, which relies on `fork()` --
Unix-only, so this harness only runs on Linux/macOS (atlas is Linux).

Cost drivers (see `generate_data.R` docstring for detail):
- `n_categories` -- distinct values per column; sets distance-matrix size.
  `method="tsp"` is an *exact* exponential-time DP solver -- only include it
  with `n_categories <= ~15` (see `TSP_SWEEP` in `sweep_config.R`).
- `n_columns` -- multiplies column-order-optimization cost when
  `column_method != "none"`. `method %in% c("greedy_wolf", "greedy_wblf")`
  only support `n_columns == 2` -- rows combining them with other column
  counts are dropped when a sweep grid is built.
- `n_rows` -- mainly controls how much of the `n_categories**n_columns`
  combinatorial space gets populated (more unique alluvia to sort).
- `column_method` -- `"tsp"`/`"neighbornet"` add a nested pass of
  edge-crossing calculations on top of `method`'s own sorting cost;
  `"none"` skips column reordering entirely.

## Running

Smoke test first (~seconds):
```
Rscript run_benchmark.R sweep --which smoke --out results/smoke.csv
```

Check grid size before committing to a long run:
```
Rscript run_benchmark.R sweep --which default --out results/default.csv --dry-run
```

Full sweep, under tmux so it survives disconnects:
```
tmux new -s wompwomp-bench
conda activate wompwomp_env
Rscript run_benchmark.R sweep --which default --out results/default.csv
# detach: Ctrl-b d
```

Resume after an interruption (skips rows already marked "ok"):
```
Rscript run_benchmark.R sweep --which default --out results/default.csv --resume
```

One-off config:
```
Rscript run_benchmark.R single --n-rows 50000 --n-columns 3 --n-categories 10 \
    --method neighbornet --column-method tsp
```

## Files

- `generate_data.R` -- synthetic wide-format categorical dataframe generator.
- `bench_utils.R` -- fork isolation + timing/memory capture.
- `sweep_config.R` -- parameter grids (`default`, `tsp`, `smoke`) and fixed params.
- `run_benchmark.R` -- CLI entry point.
- `results/` -- CSV output (gitignored).
