# wompwomp 0.99.5
Initial release and Bioconductor checks
- `calculate_objective_fenwick()`'s Binary Indexed Tree loop (the dominant cost of `neighbornet` sorting on large inputs) now runs as compiled Rcpp code instead of pure R.