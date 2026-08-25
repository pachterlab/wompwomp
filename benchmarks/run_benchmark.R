#!/usr/bin/env Rscript
# Timing/resource benchmark harness for wompwomp::sort_to_uncross.
#
# Usage
# -----
# Quick smoke test (should finish in well under a minute):
#   Rscript run_benchmark.R sweep --which smoke --out results/smoke.csv
#
# Full default sweep (long-running -- launch under tmux/screen/nohup):
#   Rscript run_benchmark.R sweep --which default --out results/default.csv
#
# Resume an interrupted sweep (skips rows already marked "ok" in --out):
#   Rscript run_benchmark.R sweep --which default --out results/default.csv --resume
#
# Single ad-hoc config:
#   Rscript run_benchmark.R single --n-rows 10000 --n-columns 3 \
#       --n-categories 8 --method neighbornet --column-method tsp
#
# Each row runs in its own forked subprocess (see bench_utils.R::run_isolated)
# so a hung or crashing config can't take down the rest of the sweep. Results
# are written to --out incrementally, one row at a time, so a killed job
# keeps whatever it already finished.

# atlas is a shared, unscheduled box -- cap BLAS/OpenMP thread pools before
# any package that might link against them gets loaded below, so this
# harness can't grab cores away from other users' running jobs. Only sets a
# var if it isn't already exported, so an operator can override by exporting
# these themselves before invoking the script.
for (.v in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_THREAD_LIMIT")) {
    if (identical(Sys.getenv(.v, unset = ""), "")) {
        do.call(Sys.setenv, setNames(list("1"), .v))
    }
}
rm(.v)

.script_args <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("^--file=", .script_args, value = TRUE)
.script_dir <- if (length(.file_arg)) dirname(sub("^--file=", "", .file_arg[1])) else "."
source(file.path(.script_dir, "generate_data.R"))
source(file.path(.script_dir, "bench_utils.R"))
source(file.path(.script_dir, "sweep_config.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a

PARAM_FIELDS <- c(
    "n_rows", "n_columns", "n_categories", "method", "column_method",
    "weight_scalar", "wt", "seed"
)
RESULT_FIELDS <- c(
    "status", "error", "wall_time_s", "user_cpu_s", "sys_cpu_s",
    "peak_rss_mb", "rss_metric", "n_unique_alluvia"
)
ALL_FIELDS <- c(PARAM_FIELDS, RESULT_FIELDS, "timeout_s", "timestamp")

# Runs inside the forked child (see bench_utils.R::run_isolated).
run_one <- function(params) {
    gen <- make_synthetic_df(
        n_rows = params$n_rows,
        n_columns = params$n_columns,
        n_categories = params$n_categories,
        seed = params$seed
    )
    n_unique_alluvia <- nrow(gen$data)

    # wt is passed as a literal (not params$wt): sort_to_uncross() resolves
    # it with rlang::ensym(), which requires a literal string/symbol at the
    # call site, not a variable holding one. Every sweep fixes wt to "value"
    # anyway (see FIXED_PARAMS in sweep_config.R).
    invisible(wompwomp::sort_to_uncross(
        data = gen$data,
        cols = gen$cols,
        wt = "value",
        method = params$method,
        column_method = params$column_method,
        weight_scalar = params$weight_scalar,
        verbose = FALSE
    ))

    list(n_unique_alluvia = n_unique_alluvia)
}

# Each config runs in a forked child (parallel::mcparallel), which inherits
# whatever is already loaded in this parent process via copy-on-write. If
# wompwomp and its dependencies (dplyr, igraph, TSP, ...) haven't been
# touched here yet, R lazy-loads their namespaces fresh in *every single
# fork* -- profiling showed this costing ~0.15-0.2s per config (most of the
# measured time for small/fast configs like greedy_wblf), which has nothing
# to do with sort_to_uncross()'s actual algorithmic cost. Run each method
# path once here so every fork starts warm.
warm_up_namespaces <- function() {
    gen <- make_synthetic_df(n_rows = 50, n_columns = 2, n_categories = 3, seed = 0)
    for (method in c("neighbornet", "greedy_wblf", "tsp")) {
        invisible(wompwomp::sort_to_uncross(
            data = gen$data, cols = gen$cols, wt = "value",
            method = method, column_method = "tsp", verbose = FALSE
        ))
    }
}

build_grid <- function(sweep_list) {
    grid <- expand.grid(sweep_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    invalid <- grid$method %in% c("greedy_wolf", "greedy_wblf") & grid$n_columns != 2
    if (any(invalid)) {
        message(sprintf(
            "dropping %d row(s): greedy_wolf/greedy_wblf require n_columns == 2",
            sum(invalid)
        ))
        grid <- grid[!invalid, , drop = FALSE]
    }
    lapply(seq_len(nrow(grid)), function(i) {
        modifyList(as.list(grid[i, , drop = FALSE]), FIXED_PARAMS)
    })
}

signature <- function(params) {
    paste(vapply(PARAM_FIELDS, function(k) as.character(params[[k]]), character(1)), collapse = "|")
}

load_completed <- function(out_path) {
    if (!file.exists(out_path)) return(character(0))
    df <- utils::read.csv(out_path, stringsAsFactors = FALSE)
    if (!"status" %in% names(df) || nrow(df) == 0) return(character(0))
    ok <- df[df$status == "ok", , drop = FALSE]
    if (nrow(ok) == 0) return(character(0))
    apply(ok[PARAM_FIELDS], 1, paste, collapse = "|")
}

run_rows <- function(rows, out_path, timeout_s, resume) {
    warm_up_namespaces()
    completed <- if (resume) load_completed(out_path) else character(0)
    write_header <- !(resume && file.exists(out_path))
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

    total <- length(rows)
    for (i in seq_along(rows)) {
        params <- rows[[i]]
        sig <- signature(params)
        if (resume && sig %in% completed) {
            cat(sprintf("[%d/%d] skip (already ok): %s\n", i, total, sig))
            next
        }
        cat(sprintf("[%d/%d] running: %s\n", i, total, sig))
        result <- run_isolated(run_one, params, timeout_s = timeout_s)
        # Seed every field as NA first: a failed/timed-out/crashed run won't
        # have produced extras like n_unique_alluvia, and as.data.frame()
        # errors outright on a genuinely missing (not just NA) list element.
        row <- setNames(as.list(rep(NA, length(ALL_FIELDS))), ALL_FIELDS)
        row <- modifyList(row, params)
        row <- modifyList(row, result)
        row$timeout_s <- timeout_s
        row$timestamp <- as.numeric(Sys.time())
        row_df <- as.data.frame(row[ALL_FIELDS], stringsAsFactors = FALSE)
        utils::write.table(
            row_df, file = out_path, sep = ",", row.names = FALSE,
            col.names = write_header, append = !write_header
        )
        write_header <- FALSE
        cat(sprintf(
            "    -> %s in %.2fs, peak_rss=%sMB\n",
            result$status, result$wall_time_s, result$peak_rss_mb
        ))
    }
}

parse_flags <- function(args, boolean_flags = character(0)) {
    out <- list()
    i <- 1
    while (i <= length(args)) {
        key <- args[i]
        if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key))
        name <- sub("^--", "", key)
        if (name %in% boolean_flags) {
            out[[name]] <- TRUE
            i <- i + 1
        } else {
            if (i == length(args)) stop(sprintf("Flag --%s requires a value", name))
            out[[name]] <- args[i + 1]
            i <- i + 2
        }
    }
    out
}

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) == 0) {
        stop("Usage: Rscript run_benchmark.R {sweep|single} [options]")
    }
    mode <- args[1]
    rest <- args[-1]

    if (mode == "sweep") {
        flags <- parse_flags(rest, boolean_flags = c("resume", "dry-run"))
        which_sweep <- flags[["which"]] %||% "smoke"
        if (!which_sweep %in% names(SWEEPS)) {
            stop(sprintf(
                "Unknown sweep '%s'. Choices: %s",
                which_sweep, paste(names(SWEEPS), collapse = ", ")
            ))
        }
        out_path <- flags[["out"]] %||% "results/sweep.csv"
        timeout_s <- as.numeric(flags[["timeout"]] %||% DEFAULT_TIMEOUT_S)
        resume <- isTRUE(flags[["resume"]])
        dry_run <- isTRUE(flags[["dry-run"]])

        rows <- build_grid(SWEEPS[[which_sweep]])
        if (dry_run) {
            cat(sprintf(
                "sweep '%s': %d configs, timeout=%ss each (worst case %.1f min)\n",
                which_sweep, length(rows), timeout_s, length(rows) * timeout_s / 60
            ))
            return(invisible())
        }
        run_rows(rows, out_path, timeout_s, resume)
    } else if (mode == "single") {
        flags <- parse_flags(rest)
        params <- list(
            n_rows = as.integer(flags[["n-rows"]] %||% 10000),
            n_columns = as.integer(flags[["n-columns"]] %||% 3),
            n_categories = as.integer(flags[["n-categories"]] %||% 8),
            method = flags[["method"]] %||% "neighbornet",
            column_method = flags[["column-method"]] %||% "none",
            weight_scalar = as.numeric(flags[["weight-scalar"]] %||% 5e5),
            wt = "value",
            seed = 0
        )
        out_path <- flags[["out"]] %||% "results/single.csv"
        timeout_s <- as.numeric(flags[["timeout"]] %||% DEFAULT_TIMEOUT_S)
        run_rows(list(params), out_path, timeout_s, resume = FALSE)
    } else {
        stop(sprintf("Unknown mode '%s'. Must be 'sweep' or 'single'.", mode))
    }
}

main()
