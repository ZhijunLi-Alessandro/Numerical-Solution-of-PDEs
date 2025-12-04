"""
@file visualize_Parabolic_ADI.py
@brief Visualization helpers for ADI parabolic outputs.

This script creates animations and error visualizations for datasets produced
by the ADI parabolic example. It uses the `visualize_utils` helpers to load
CSV outputs and produce GIF animations.

@author Li Zhijun
@date 2025-12-03
"""

import numpy as np
import matplotlib.pyplot as plt
from visualize_utils import load_solution_from_csv, visualize_solution, animate_time_series, animate_time_series_difference

grid_data = load_solution_from_csv('results/Parabolic/data/ADI/grid_data.csv')
grid_data[20, 60] = 0

# animate_time_series("results/Parabolic/data/ADI", "exact", grid_data, list(range(20, 2001, 20)), "Exact", save_path="results/Parabolic/exact.gif")
# animate_time_series("results/Parabolic/data/ADI", "solution", grid_data, list(range(20, 2001, 20)), "Solution", save_path="results/Parabolic/ADI.gif")
animate_time_series_difference("results/Parabolic/data/ADI", "exact", "solution", grid_data, list(range(20, 2001, 20)), "Error", save_path="results/Parabolic/ADI_error_without_center.gif")