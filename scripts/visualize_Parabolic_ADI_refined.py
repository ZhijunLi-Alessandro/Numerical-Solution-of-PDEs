"""
@file visualize_Parabolic_ADI_refined.py
@brief Visualization helpers for refined ADI parabolic outputs.

This script creates difference/error animations for the refined ADI dataset.
It uses functions from `visualize_utils` to load CSV outputs and produce GIFs.

@author Li Zhijun
@date 2025-12-03
"""

import numpy as np
import matplotlib.pyplot as plt
from visualize_utils import load_solution_from_csv, visualize_solution, animate_time_series, animate_time_series_difference

grid_data = load_solution_from_csv('results/Parabolic/data/ADI_refined/grid_data.csv')
grid_data[40, 120] = 0

# animate_time_series("results/Parabolic/data/ADI_refined", "exact", grid_data, list(range(80, 8001, 80)), "Exact", save_path="results/Parabolic/exact_refined.gif")
# animate_time_series("results/Parabolic/data/ADI_refined", "solution", grid_data, list(range(80, 8001, 80)), "Solution", save_path="results/Parabolic/ADI_refined.gif")
animate_time_series_difference("results/Parabolic/data/ADI_refined", "exact", "solution", grid_data, list(range(80, 8001, 80)), "Error", save_path="results/Parabolic/ADI_error_refined_without_center.gif")