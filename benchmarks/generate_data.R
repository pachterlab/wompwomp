# Synthetic wide-format categorical data for wompwomp benchmarks.
#
# sort_to_uncross()'s cost is driven mainly by:
#   - n_categories: distinct values per graphing column. Sets the distance-
#     matrix size for method %in% c("neighbornet", "tsp"). "tsp" is an
#     *exact* DP solver (exponential in node count) -- keep n_categories
#     small (~15) when including it in a sweep.
#   - n_columns: number of graphing columns. Multiplies column-order-
#     optimization cost (column_method != "none"), since each column pair
#     gets its own edge-crossing calculation. Note: method %in%
#     c("greedy_wolf", "greedy_wblf") only support n_columns == 2.
#   - n_rows: raw observation count, pre-aggregation. sort_to_uncross()
#     operates on unique column-combinations (rows of the wt = "value"
#     table), so n_rows mainly controls how much of the
#     n_categories^n_columns combinatorial space gets populated.

make_synthetic_df <- function(n_rows, n_columns, n_categories, seed = 0) {
    set.seed(seed)
    cols <- paste0("col_", seq_len(n_columns) - 1L)
    raw <- as.data.frame(
        lapply(cols, function(cn) {
            sample(as.character(seq_len(n_categories) - 1L), size = n_rows, replace = TRUE)
        }),
        stringsAsFactors = FALSE
    )
    names(raw) <- cols
    agg <- as.data.frame(dplyr::count(raw, dplyr::across(dplyr::all_of(cols)), name = "value"))
    list(data = agg, cols = cols)
}
