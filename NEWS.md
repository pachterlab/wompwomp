# wompwomp 0.99.5
Initial release and Bioconductor checks
- The `biowomp` package has been merged in: `plot_alluvial()` (build a data frame with `prep_for_lodes()` + `sort_to_uncross()`, colour it, and render with `ggalluvial`) now lives in `wompwomp`. `ggfittext` and `ggrastr` are optional (`Suggests`).
- `sort_to_uncross()`'s default `method` is now `"neighbornet"` (was `"tsp"`), matching the algorithm described in the paper and the Python implementation.
- Within a stratum, alluvia are now ordered by a fully-specified key (current axis, then nearest-right axes, then nearest-left axes) so the crossing objective no longer depends on input row order for 3+ layers.
- Fixed `column_method = "random"` in `sort_to_uncross()` (previously always errored).
- `get_lode_clusters(method = "right")` now propagates colours right-to-left (it was a duplicate of `"left"`), and no longer errors on a plain `data.frame` input.
- `find_colors_advanced()` no longer breaks when a stratum name contains an underscore.
- `calculate_objective_fenwick()`'s Binary Indexed Tree loop (the dominant cost of `neighbornet` sorting on large inputs) now runs as compiled Rcpp code instead of pure R.
- Added `method = "barycenter"` and `method = "median"` to `sort_to_uncross()`: classic Sugiyama-style two-layer crossing-reduction heuristics, cheaper (O(n log n) per pass) than `greedy_wolf`/`greedy_wblf` (O(n1*n2)). One-directional variants `barycenter_one_sided`/`median_one_sided` (paired with `fixed_column`, like `greedy_wolf`) are also available alongside the two-pass alternating default.