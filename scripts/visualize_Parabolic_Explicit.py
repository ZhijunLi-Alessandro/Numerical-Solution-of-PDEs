"""
@file visualize_Parabolic_Explicit.py
@brief Visualization helpers for the parabolic explicit example.

This script builds animations and error visualizations for data produced by
the explicit parabolic example. It uses helpers from `visualize_utils` and
writes GIFs for inspection.

@author Li Zhijun
@date 2025-12-03
"""

import numpy as np
import matplotlib.pyplot as plt
from visualize_utils import load_solution_from_csv, visualize_solution, animate_time_series, animate_time_series_difference

grid_data = load_solution_from_csv('results/Parabolic/data/Explicit/grid_data.csv')
grid_data[20, 60] = 0

# animate_time_series("results/Parabolic/data/Explicit", "exact", grid_data, list(range(100, 10001, 100)), "Exact", save_path="results/Parabolic/exact.gif")
# animate_time_series("results/Parabolic/data/Explicit", "solution", grid_data, list(range(100, 10001, 100)), "Solution", save_path="results/Parabolic/Explicit.gif")
animate_time_series_difference("results/Parabolic/data/Explicit", "exact", "solution", grid_data, list(range(100, 10001, 100)), "Error", save_path="results/Parabolic/Explicit_error_without_center.gif")