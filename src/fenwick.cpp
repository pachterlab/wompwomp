#include <Rcpp.h>
using namespace Rcpp;

// Binary Indexed Tree (Fenwick tree) core loop for calculate_objective_fenwick()
// (R/objective_calculation.R). Direct port of that function's steps 3+: same
// BIT size (max_rank + 2), same two query passes, same update pass, same
// accumulation order, so results are bit-for-bit identical to the R version.
// `idx & (-idx)` isolates the lowest set bit (two's-complement, matching R's
// bitwAnd(idx, -idx) for the positive idx values used here) -- this alone
// was ~37% of R's wall time when called as bitwAnd(), a full R function call
// per Fenwick-tree step.
//
// y2_rank must already be sorted by y1 and rank-compressed to 1..max_rank
// (done in R -- those steps are cheap, vectorized, and not the bottleneck).
//
// [[Rcpp::export]]
double calculate_objective_fenwick_cpp(IntegerVector y2_rank, NumericVector weight, bool weighted_metric) {
    int n = y2_rank.size();
    int max_rank = max(y2_rank);
    int size = max_rank + 2;

    std::vector<int> tree_count(size, 0);
    std::vector<double> tree_weight(size, 0.0);

    double total_cross_weight = 0.0;

    for (int i = 0; i < n; i++) {
        int r = y2_rank[i];
        double w = weight[i];

        int idx = max_rank + 1;
        int count_hi = 0;
        double wsum_hi = 0.0;
        while (idx > 0) {
            count_hi += tree_count[idx];
            wsum_hi += tree_weight[idx];
            idx -= (idx & (-idx));
        }

        idx = r + 1;
        int count_lo = 0;
        double wsum_lo = 0.0;
        while (idx > 0) {
            count_lo += tree_count[idx];
            wsum_lo += tree_weight[idx];
            idx -= (idx & (-idx));
        }

        if (weighted_metric) {
            total_cross_weight += w * (wsum_hi - wsum_lo);
        } else {
            total_cross_weight += (count_hi - count_lo);
        }

        idx = r + 1;
        double upd_w = weighted_metric ? w : 1.0;
        while (idx < size) {
            tree_count[idx] += 1;
            tree_weight[idx] += upd_w;
            idx += (idx & (-idx));
        }
    }

    return total_cross_weight;
}
