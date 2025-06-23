# calculate_objective.py

import pandas as pd
from scipy.stats import rankdata

class BIT:
    def __init__(self, size):
        self.tree_count = [0] * (size + 1)
        self.tree_weight = [0.0] * (size + 1)

    def update(self, index, weight):
        index += 1
        while index < len(self.tree_count):
            self.tree_count[index] += 1
            self.tree_weight[index] += weight
            index += index & -index

    def query(self, index):
        count = 0
        weight_sum = 0.0
        index += 1
        while index > 0:
            count += self.tree_count[index]
            weight_sum += self.tree_weight[index]
            index -= index & -index
        return count, weight_sum

    def query_range(self, low, high):
        c1, w1 = self.query(high)
        c2, w2 = self.query(low - 1)
        return c1 - c2, w1 - w2

def calculate_objective_fenwick(df):
    # Step 1: Sort by y1
    df_sorted = df.sort_values('y1').reset_index(drop=True)

    # Step 2: Rank-compress y2 (higher y2 → higher rank)
    df_sorted['y2_rank'] = rankdata(df_sorted['y2'], method='dense').astype(int)
    max_rank = df_sorted['y2_rank'].max()

    # Step 3: Use BIT
    bit = BIT(size=max_rank + 2)
    total_cross_weight = 0.0

    for y2_rank, weight in zip(df_sorted['y2_rank'].values, df_sorted['count'].values):
        # Count previous y2s > current (strictly greater)
        _, sum_weights_above = bit.query_range(y2_rank + 1, max_rank)
        total_cross_weight += weight * sum_weights_above

        # Add current y2_rank to BIT
        bit.update(y2_rank, weight)
    
    return total_cross_weight
