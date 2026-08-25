#!/usr/bin/env Rscript
# Summary figure for a sweep CSV produced by run_benchmark.R, in the style of
# wompywompy's benchmarks/results/sweep_report.png: wall time and peak RSS
# each plotted against n_columns / n_categories / n_unique_alluvia, points
# colored/shaped by algorithm, with a short tick marking each x-group's mean.
#
# Usage: Rscript plot_results.R results/default.csv results/sweep_report.png
#
# Base graphics only (no ggplot2/etc.) so this has no extra package
# dependencies beyond what run_benchmark.R already needs.

args <- commandArgs(trailingOnly = TRUE)
in_path <- if (length(args) >= 1) args[1] else "results/default.csv"
out_path <- if (length(args) >= 2) args[2] else "results/sweep_report.png"

df <- utils::read.csv(in_path, stringsAsFactors = FALSE)
df <- df[df$status == "ok", , drop = FALSE]
if (nrow(df) == 0) stop(sprintf("no 'ok' rows in %s", in_path))

methods <- sort(unique(df$method))
palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
shapes <- c(16, 17, 15, 18)
col_for <- setNames(palette[seq_along(methods)], methods)
pch_for <- setNames(shapes[seq_along(methods)], methods)

rss_metric_used <- if ("rss_metric" %in% names(df)) unique(stats::na.omit(df$rss_metric)) else character(0)
has_rss <- length(rss_metric_used) > 0 && any(!is.na(df$peak_rss_mb))
# "vmhwm" (Linux /proc/self/status) is a true high-water mark; "ps_snapshot"
# (macOS fallback, see bench_utils.R::.rss_info) is current RSS sampled right
# after the timed work finishes -- a same-run approximation, not a tracked
# peak, so labelled and captioned differently rather than presented as
# equivalent.
rss_is_snapshot <- has_rss && "ps_snapshot" %in% rss_metric_used
rss_ylab <- if (rss_is_snapshot) "RSS at completion (MB)" else "peak RSS (MB)"

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

grDevices::png(out_path, width = 2000, height = 1270, res = 150)
on.exit(grDevices::dev.off(), add = TRUE)

layout(
    matrix(c(1, 1, 1, 2, 3, 4, 5, 6, 7), nrow = 3, byrow = TRUE),
    heights = c(0.5, 2, 2)
)
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.01, 0.75, "wompwomp benchmark sweep — results", cex = 1.8, font = 2, adj = 0)
subtitle <- sprintf(
    "%d-config sweep. Dot = one config, tick = group mean. n_unique_alluvia (log x-axis) is the derived quantity that actually drives cost.",
    nrow(df)
)
if (rss_is_snapshot) {
    subtitle <- paste0(
        subtitle,
        "\nMemory panels: current RSS sampled right after each run finished (macOS has no /proc/self/status high-water mark) -- a same-run approximation, not a tracked peak."
    )
}
text(0.01, 0.25, subtitle, cex = 1.0, col = "grey30", adj = 0)
legend(
    "topright", legend = methods, col = col_for[methods], pch = pch_for[methods],
    horiz = TRUE, bty = "n", cex = 1.1, pt.cex = 1.3
)

par(mar = c(4, 4.5, 2.5, 1))
panel(df$n_columns, df$wall_time_s, df$method, "n_columns", "wall time (s)", log = "y")
title("wall time vs n_columns", adj = 0, font.main = 2, cex.main = 1.1)
panel(df$n_categories, df$wall_time_s, df$method, "n_categories", "wall time (s)", log = "xy")
title("wall time vs n_categories", adj = 0, font.main = 2, cex.main = 1.1)
panel(df$n_unique_alluvia, df$wall_time_s, df$method, "n_unique_alluvia", "wall time (s)", log = "xy")
title("wall time vs n_unique_alluvia", adj = 0, font.main = 2, cex.main = 1.1)

rss_msg <- if (has_rss) NULL else "RSS unavailable on this platform\n(no /proc/self/status and `ps` lookup failed)"
mem_title_suffix <- if (rss_is_snapshot) " (snapshot)" else ""
panel(df$n_columns, df$peak_rss_mb, df$method, "n_columns", rss_ylab, unavailable_msg = rss_msg)
title(paste0("memory vs n_columns", mem_title_suffix), adj = 0, font.main = 2, cex.main = 1.1)
panel(df$n_categories, df$peak_rss_mb, df$method, "n_categories", rss_ylab, log = "x", unavailable_msg = rss_msg)
title(paste0("memory vs n_categories", mem_title_suffix), adj = 0, font.main = 2, cex.main = 1.1)
panel(df$n_unique_alluvia, df$peak_rss_mb, df$method, "n_unique_alluvia", rss_ylab, log = "x", unavailable_msg = rss_msg)
title(paste0("memory vs n_unique_alluvia", mem_title_suffix), adj = 0, font.main = 2, cex.main = 1.1)

cat(sprintf("wrote %s (%d ok rows of %d total)\n", out_path, nrow(df), nrow(utils::read.csv(in_path))))
