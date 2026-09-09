#!/usr/bin/env Rscript
# Summary figures for a sweep CSV produced by run_benchmark.R, in the style of
# wompywompy's benchmarks/results/sweep_report.png: wall time and peak RSS
# each plotted against sweep axes, points colored/shaped by algorithm, with a
# short tick marking each x-group's mean.
#
# run_benchmark.R also records a crossing-count "objective" column, but it's
# deliberately not plotted here: at a fixed n_categories/n_columns, the
# weighted objective scales roughly quadratically with n_rows, which swamps
# the (much smaller) between-method spread on any shared log-scale axis --
# the resulting panel made every method look identical rather than showing
# real quality differences.
#
# Produces two figures, since most methods (greedy_wolf/greedy_wblf/
# barycenter/median/barycenter_one_sided/median_one_sided) only ever run at
# n_columns == 2 (see sweep_config.R):
#   - "<out>_2layer.png": every method, restricted to its n_columns == 2
#     rows, plotted against n_categories / n_unique_alluvia. n_columns is
#     constant here, so it isn't a useful x-axis -- this is the figure where
#     all methods (including the 2-column-only ones) can be compared
#     head-to-head.
#   - "<out>_multilayer.png": only methods actually swept across more than
#     one n_columns value (e.g. neighbornet, tsp), plotted against
#     n_columns / n_categories / n_unique_alluvia as in the original single
#     figure.
#
# Usage: Rscript plot_results.R results/default.csv results/sweep_report.png
#
# Base graphics only (no ggplot2/etc.) so this has no extra package
# dependencies beyond what run_benchmark.R already needs.

args <- commandArgs(trailingOnly = TRUE)
in_path <- if (length(args) >= 1) args[1] else "results/default.csv"
out_path <- if (length(args) >= 2) args[2] else "results/sweep_report.png"

# "results/sweep_report.png" + "2layer" -> "results/sweep_report_2layer.png"
suffixed_path <- function(path, suffix) {
    ext <- sub("^.*?(\\.[^./]+)$", "\\1", path)
    if (identical(ext, path)) ext <- ""
    base <- if (nzchar(ext)) substr(path, 1, nchar(path) - nchar(ext)) else path
    paste0(base, "_", suffix, ext)
}

df_all <- utils::read.csv(in_path, stringsAsFactors = FALSE)
df_all <- df_all[df_all$status == "ok", , drop = FALSE]
if (nrow(df_all) == 0) stop(sprintf("no 'ok' rows in %s", in_path))

# Color/shape per method are assigned once from every method present
# anywhere in the sweep, so a given method keeps the same look across both
# figures even though each only plots a subset of methods.
methods_all <- sort(unique(df_all$method))
palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
shapes <- c(16, 17, 15, 18, 1, 2, 0, 5)
if (length(methods_all) > length(palette)) {
    stop(sprintf(
        "palette only covers %d methods; sweep has %d (%s) -- add more colors/shapes",
        length(palette), length(methods_all), paste(methods_all, collapse = ", ")
    ))
}
col_for <- setNames(palette[seq_along(methods_all)], methods_all)
pch_for <- setNames(shapes[seq_along(methods_all)], methods_all)

rss_metric_used <- if ("rss_metric" %in% names(df_all)) unique(stats::na.omit(df_all$rss_metric)) else character(0)
has_rss <- length(rss_metric_used) > 0 && any(!is.na(df_all$peak_rss_mb))
# "vmhwm" (Linux /proc/self/status) is a true high-water mark; "ps_snapshot"
# (macOS fallback, see bench_utils.R::.rss_info) is current RSS sampled right
# after the timed work finishes -- a same-run approximation, not a tracked
# peak, so labelled and captioned differently rather than presented as
# equivalent.
rss_is_snapshot <- has_rss && "ps_snapshot" %in% rss_metric_used
rss_ylab <- if (rss_is_snapshot) "RSS at completion (MB)" else "peak RSS (MB)"
rss_msg <- if (has_rss) NULL else "RSS unavailable on this platform\n(no /proc/self/status and `ps` lookup failed)"

# One tick per (method, x-value) at the group's mean y, echoing the python
# report's "dot = one config, tick = group mean".
group_means <- function(x, y, method) {
    key <- paste(method, x)
    m <- tapply(y, key, mean, na.rm = TRUE)
    ux <- tapply(x, key, `[`, 1)
    umethod <- tapply(method, key, `[`, 1)
    data.frame(x = ux[names(m)], y = as.numeric(m), method = umethod[names(m)])
}

panel <- function(x, y, method, xlab, ylab, log = "", unavailable_msg = NULL) {
    if (!is.null(unavailable_msg)) {
        plot(1, 1, type = "n", axes = FALSE, xlab = xlab, ylab = ylab, main = "")
        box()
        text(1, 1, unavailable_msg, cex = 0.85, col = "grey40")
        return(invisible())
    }
    ok <- !is.na(x) & !is.na(y) & x > 0 & (!grepl("x", log) | x > 0) & (!grepl("y", log) | y > 0)
    x <- x[ok]; y <- y[ok]; method <- method[ok]
    plot(
        x, y, log = log, col = col_for[method], pch = pch_for[method],
        xlab = xlab, ylab = ylab, cex = 1.1,
        panel.first = grid(col = "grey88", lty = 1)
    )
    gm <- group_means(x, y, method)
    seg_halfwidth <- diff(range(log10(x[x > 0]))) * 0.015 + 0
    for (i in seq_len(nrow(gm))) {
        gx <- gm$x[i]
        if (grepl("x", log)) {
            x0 <- gx * 10^(-seg_halfwidth); x1 <- gx * 10^(seg_halfwidth)
        } else {
            x0 <- gx - seg_halfwidth * mean(x); x1 <- gx + seg_halfwidth * mean(x)
        }
        segments(x0, gm$y[i], x1, gm$y[i], col = col_for[gm$method[i]], lwd = 3)
    }
}

# Every metric row shown in a figure: which df column, its y-axis label,
# whether it's log-scaled, and an optional "unavailable" fallback message.
metrics <- list(
    list(col = "wall_time_s", label = "wall time", ylab = "wall time (s)", log_y = TRUE, unavailable_msg = NULL),
    list(col = "peak_rss_mb", label = paste0("memory", if (rss_is_snapshot) " (snapshot)" else ""), ylab = rss_ylab, log_y = FALSE, unavailable_msg = rss_msg)
)

# Renders one figure: a title/legend strip over one row per entry of
# `metrics`, one panel per entry of `x_specs` (list of list(col=<df column
# name>, log_x=<TRUE/FALSE>)). Skips (with a message) if `df` has no rows,
# since a sweep can legitimately produce an empty 2-layer or multi-layer
# subset (e.g. a smoke sweep only ever runs n_columns == 2).
render_figure <- function(df, out_path, title, x_specs) {
    if (nrow(df) == 0) {
        message(sprintf("skipping %s: no rows to plot", out_path))
        return(invisible())
    }
    methods <- sort(unique(df$method))
    n <- length(x_specs)
    n_metrics <- length(metrics)
    letters_used <- LETTERS[seq_len(n_metrics * n)]

    # A legend row fits ~4 method labels before wrapping; grow the title
    # strip's share of the layout to fit however many rows that takes (a
    # 2-layer figure's legend, unlike the 3-column multi-layer one, has to
    # hold every method in the sweep -- often more than fits on one line).
    legend_ncol <- min(length(methods), 4L)
    n_legend_rows <- ceiling(length(methods) / legend_ncol)
    title_height <- 0.14 + 0.11 * n_legend_rows

    # Keep each metric row's pixel height roughly constant regardless of how
    # many metric rows a figure has, rather than squeezing a fixed device
    # height across however many rows there are.
    px_per_layout_unit <- 295
    device_height <- round(px_per_layout_unit * (title_height + 2 * n_metrics))

    grDevices::png(out_path, width = max(1500, 700 * n), height = device_height, res = 150)
    on.exit(grDevices::dev.off(), add = TRUE)

    layout(
        rbind(rep(1L, n), matrix(seq_len(n_metrics * n) + 1L, nrow = n_metrics, ncol = n, byrow = TRUE)),
        heights = c(title_height, rep(2, n_metrics))
    )
    par(mar = c(0, 0, 0, 0))
    plot.new()
    text(0.01, 0.92, title, cex = 1.8, font = 2, adj = c(0, 1))
    legend(
        "bottom", legend = methods, col = col_for[methods], pch = pch_for[methods],
        ncol = legend_ncol, bty = "n", cex = 1.0, pt.cex = 1.2
    )

    par(mar = c(4, 4.5, 2.5, 1))
    idx <- 1L
    for (met in metrics) {
        for (spec in x_specs) {
            log_x <- if (spec$log_x) "x" else ""
            log_full <- if (met$log_y) paste0(log_x, "y") else log_x
            panel(df[[spec$col]], df[[met$col]], df$method, spec$col, met$ylab, log = log_full, unavailable_msg = met$unavailable_msg)
            title(sprintf("%s  %s vs %s", letters_used[idx], met$label, spec$col), adj = 0, font.main = 2, cex.main = 1.1)
            idx <- idx + 1L
        }
    }

    cat(sprintf(
        "wrote %s (%d rows, methods: %s)\n",
        out_path, nrow(df), paste(methods, collapse = ", ")
    ))
    invisible()
}

# Methods that were actually swept across more than one n_columns value
# (e.g. neighbornet, tsp) vs. methods only ever run at n_columns == 2
# (greedy_wolf/greedy_wblf/barycenter/median/*_one_sided -- see
# sweep_config.R's comment on why).
n_cols_per_method <- tapply(df_all$n_columns, df_all$method, function(x) length(unique(x)))
multi_col_methods <- names(n_cols_per_method[n_cols_per_method > 1])

df_2layer <- df_all[df_all$n_columns == 2, , drop = FALSE]
df_multi <- df_all[df_all$method %in% multi_col_methods, , drop = FALSE]

render_figure(
    df_2layer, suffixed_path(out_path, "2layer"),
    "wompwomp benchmark sweep — 2-layer comparison",
    list(
        list(col = "n_categories", log_x = TRUE),
        list(col = "n_unique_alluvia", log_x = TRUE)
    )
)

render_figure(
    df_multi, suffixed_path(out_path, "multilayer"),
    "wompwomp benchmark sweep — multi-layer scaling",
    list(
        list(col = "n_columns", log_x = FALSE),
        list(col = "n_categories", log_x = TRUE),
        list(col = "n_unique_alluvia", log_x = TRUE)
    )
)
