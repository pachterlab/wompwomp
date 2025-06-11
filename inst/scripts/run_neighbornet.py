# neighbor_net_wrapper.py

from typing import Tuple
from splitspy.nnet import nnet_cycle, nnet_splits
import numpy as np
import random

def __setup_matrix_custom(labels: [str], matrix: [float]) -> np.array:
    n = len(labels)
    max_number_of_nodes = max(3, 3 * n - 5)

    # mat = np.empty(((max_number_of_nodes + 1), (max_number_of_nodes + 1)))  # JMR, commented out (to remove randomness)
    mat = np.zeros((max_number_of_nodes + 1, max_number_of_nodes + 1), dtype=np.float64)  # JMR, commented out (to remove randomness)

    for i in range(0, max_number_of_nodes):
        for j in range(0, max_number_of_nodes):
            if i < n and j < n:
                mat[i+1][j+1] = matrix[i][j]

    return mat

def compute_custom(labels: [str], matrix) -> [int]:
    n = len(labels)

    if n <= 3:
        return list(range(0, n + 1))

    nodes_head = nnet_cycle.__setup_nodes(n)

    mat = __setup_matrix_custom(labels, matrix.copy())  # matrix is 0-based, mat is 1-based

    joins = nnet_cycle.__join_nodes(n, mat, nodes_head)

    cycle = nnet_cycle.__expand_nodes(joins, nodes_head)

    cycle = nnet_cycle.__normalize_cycle(cycle)

    return cycle

def neighbor_net(labels: list[str], mat: list[list[float]], cutoff=0.0001, constrained=True, compute_splits=False) -> Tuple[list, list]:
    mat = mat.astype(np.float64)
    cycle = compute_custom(labels, mat)
    if compute_splits:
        splits = nnet_splits.compute(len(labels), mat, cycle, cutoff, constrained)
    else:
        splits = []
    return cycle, splits
