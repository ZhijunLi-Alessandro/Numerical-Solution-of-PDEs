"""
@file visualize_Parabolic_Explicit_refined.py
@brief Visualization helpers for refined parabolic explicit outputs.

This script creates animations for the refined explicit parabolic dataset and
computes difference/error visualizations. It uses `visualize_utils` helpers.

@author Li Zhijun
@date 2025-12-03
"""

import numpy as np
import matplotlib.pyplot as plt
from visualize_utils import load_solution_from_csv, visualize_solution, animate_time_series, animate_time_series_difference

grid_data = load_solution_from_csv('results/Parabolic/data/Explicit_refined/grid_data.csv')
grid_data[40, 120] = 0

# animate_time_series("results/Parabolic/data/Explicit_refined", "exact", grid_data, list(range(400, 40001, 400)), "Exact", save_path="results/Parabolic/exact_refined.gif")
# animate_time_series("results/Parabolic/data/Explicit_refined", "solution", grid_data, list(range(400, 40001, 400)), "Solution", save_path="results/Parabolic/Explicit_refined.gif")
animate_time_series_difference("results/Parabolic/data/Explicit_refined", "exact", "solution", grid_data, list(range(400, 40001, 400)), "Error", save_path="results/Parabolic/Explicit_error_refined_without_center.gif")