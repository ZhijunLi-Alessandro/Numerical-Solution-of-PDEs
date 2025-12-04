"""
@file visualize_Verify.py
@brief Visualization tools for parabolic verification outputs.

This module creates animations and plots for verifying parabolic solver outputs
against analytical/semi-analytical reference data. It relies on the
`visualize_utils` helpers for loading and plotting CSV datasets.

@author Li Zhijun
@date 2025-12-03
"""

import numpy as np
import matplotlib.pyplot as plt
from visualize_utils import load_solution_from_csv, visualize_solution, animate_time_series, animate_time_series_difference

grid_data = load_solution_from_csv('results/Parabolic/data/exact_test_refined/grid_data.csv')

animate_time_series("results/Parabolic/data/exact_test_refined", "exact", grid_data, list(range(400, 20001, 400)), "Exact", save_path="results/Parabolic/solution_refined.gif")
animate_time_series("results/Parabolic/data/exact_test_refined", "rhs", grid_data, list(range(400, 20001, 400)), "Solution", save_path="results/Parabolic/rhs_refined.gif")