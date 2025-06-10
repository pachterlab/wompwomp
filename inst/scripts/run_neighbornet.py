# neighbor_net_wrapper.py

from typing import Tuple
from splitspy.nnet import nnet_cycle, nnet_splits
import numpy as np

def neighbor_net(labels: list[str], mat: list[list[float]], cutoff=0.0001, constrained=True) -> Tuple[list, list]:
    cycle = nnet_cycle.compute(labels, mat)
    splits = nnet_splits.compute(len(labels), mat, cycle, cutoff, constrained)
    return cycle, splits
