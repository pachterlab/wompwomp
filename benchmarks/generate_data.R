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
#     c("greedy_wolf", "greedy_wblf", "barycenter", "median",
#     "barycenter_one_sided", "median_one_sided") only support n_columns == 2.
#   - n_rows: raw observation count, pre-aggregation. sort_to_uncross()
#     operates on unique column-combinations (rows of the wt = "value"
#     table), so n_rows mainly controls how much of the
#     n_categories^n_columns combinatorial space gets populated.
#
# Categories are correlated across columns via a shared latent cluster label
# per row (see `cluster_strength` below), not sampled independently/
# uniformly. Uniform-random categories have no structure for a smarter
# heuristic to exploit -- every method (tsp/neighbornet down to
# barycenter/median) ends up within a few percent of the same crossing
# objective, which is realistic for adversarial/structureless data but not
# for the clustered data (e.g. tissue/cell-type combinations) sort_to_uncross()
# is actually built for. Correlated categories give methods that search
# harder (tsp, neighbornet) real room to beat cheap heuristics on objective
# quality, which is the point of recording `objective` in the sweep at all.

make_synthetic_df <- function(n_rows, n_columns, n_categories, seed = 0, cluster_strength = 0.7) {
    set.seed(seed)
    cols <- paste0("col_", seq_len(n_columns) - 1L)
    labels <- as.character(seq_len(n_categories) - 1L)

    # One latent cluster label per row, shared across every column: with
    # probability `cluster_strength` a column's value matches the row's
    # cluster, otherwise it's uniform noise. All columns correlate with the
    # same latent variable (not just pairwise), mirroring how e.g.
    # tissue/cell-type/sex axes co-vary with an underlying biological
    # identity in real alluvial data. Each column maps the shared cluster
    # through its OWN random permutation of that column's labels -- not the
    # literal same label string -- so the correct cross-column alignment is
    # an unknown permutation a method has to find, rather than the trivial
    # identity that alphabetical default_sorting would recover for free.
    cluster <- sample.int(n_categories, size = n_rows, replace = TRUE)
    col_perms <- replicate(n_columns, sample(labels), simplify = FALSE)
    raw <- as.data.frame(
        lapply(seq_len(n_columns), function(j) {
            matches_cluster <- stats::runif(n_rows) < cluster_strength
            noise <- sample(labels, size = n_rows, replace = TRUE)
            ifelse(matches_cluster, col_perms[[j]][cluster], noise)
        }),
        stringsAsFactors = FALSE
    )
    names(raw) <- cols
    agg <- as.data.frame(dplyr::count(raw, dplyr::across(dplyr::all_of(cols)), name = "value"))
    list(data = agg, cols = cols)
}
