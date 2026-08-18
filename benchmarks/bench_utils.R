# Fork-isolated timing/resource harness.
#
# Each benchmark config runs in its own forked child process (parallel::
# mcparallel) so that:
#   - a hung config (e.g. exact-TSP sorting on too many categories) can be
#     killed on a timeout instead of blocking the rest of the sweep.
#   - a crashing config doesn't take down the sweep process.
#   - peak RSS reflects that one run's growth, not accumulated state left
#     over from earlier configs in the same process.
#
# Unix-only: relies on fork() (via parallel::mcparallel) and reads memory
# usage from /proc/self/status on Linux (atlas is the intended host) or via
# `ps` on macOS, so this only runs on Linux/macOS. Peak RSS includes whatever
# was already resident in the parent at fork time (loaded packages, etc.),
# so treat absolute numbers as having a shared baseline offset --
# differences *between* configs are what's meaningful.

# Returns list(mb = <number or NA>, metric = <string>). On Linux this reads
# /proc/self/status's VmHWM, a true high-water mark tracked by the kernel
# across the process's lifetime. macOS has no equivalent exposed to R, so
# there we fall back to `ps`'s current RSS sampled right after the timed
# work finishes -- a same-run snapshot, not a true peak (it can undercount
# a transient spike that was already freed by the time of sampling). The
# metric field says which one a given row actually is; don't compare
# "vmhwm" and "ps_snapshot" values against each other as if equivalent.
.rss_info <- function() {
    status_path <- "/proc/self/status"
    if (file.exists(status_path)) {
        lines <- readLines(status_path, warn = FALSE)
        hwm_line <- grep("^VmHWM:", lines, value = TRUE)
        if (length(hwm_line) > 0) {
            kb <- as.numeric(regmatches(hwm_line, regexpr("[0-9]+", hwm_line)))
            return(list(mb = kb / 1024, metric = "vmhwm"))
        }
    }
    rss_kb <- suppressWarnings(as.numeric(
        system2("ps", c("-o", "rss=", "-p", Sys.getpid()), stdout = TRUE)
    ))
    if (length(rss_kb) == 1 && !is.na(rss_kb)) {
        return(list(mb = rss_kb / 1024, metric = "ps_snapshot"))
    }
    list(mb = NA_real_, metric = NA_character_)
}

#' Run target_fn(params) in a forked subprocess with a wall-clock timeout.
#'
#' target_fn must take a single list argument. It may return a list of extra
#' fields to merge into the result (e.g. n_unique_alluvia); anything it
#' raises is captured as status = "error".
run_isolated <- function(target_fn, params, timeout_s = 600) {
    start <- proc.time()[["elapsed"]]

    job <- parallel::mcparallel({
        t0 <- proc.time()
        outcome <- tryCatch(
            list(status = "ok", error = NA_character_, extra = target_fn(params)),
            error = function(e) list(status = "error", error = conditionMessage(e), extra = list())
        )
        t1 <- proc.time()
        rss <- .rss_info()
        c(
            list(
                status = outcome$status,
                error = outcome$error,
                user_cpu_s = unname((t1 - t0)[["user.self"]]),
                sys_cpu_s = unname((t1 - t0)[["sys.self"]]),
                peak_rss_mb = rss$mb,
                rss_metric = rss$metric
            ),
            outcome$extra
        )
    }, silent = TRUE)

    collected <- parallel::mccollect(job, wait = FALSE, timeout = timeout_s)
    wall_time_s <- proc.time()[["elapsed"]] - start

    if (is.null(collected) || length(collected) == 0) {
        tools::pskill(job$pid, tools::SIGKILL)
        parallel::mccollect(job) # reap the zombie
        return(list(
            status = "timeout", error = sprintf("exceeded %ss", timeout_s),
            wall_time_s = unname(wall_time_s), user_cpu_s = NA, sys_cpu_s = NA,
            peak_rss_mb = NA, rss_metric = NA
        ))
    }

    res <- collected[[1]]
    if (is.null(res) || inherits(res, "try-error")) {
        return(list(
            status = "crashed", error = "child process produced no result",
            wall_time_s = unname(wall_time_s), user_cpu_s = NA, sys_cpu_s = NA,
            peak_rss_mb = NA, rss_metric = NA
        ))
    }

    res$wall_time_s <- unname(wall_time_s)
    res
}
