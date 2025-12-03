## @file visualize_Neumann.py
#  @brief A python script for visualizing the calculation results of the Neumann problem
#  
#  This script is used to process and analyze data and generate figures.
#  Includes data loading and preprocessing, generate multiple figures.
#  @see Dirichlet.c
#  @author Li Zhijun
#  @date 2025-10-29
import numpy as np
import matplotlib.pyplot as plt
from visualize_utils import load_solution_from_csv, visualize_solution, animate_time_series, animate_time_series_difference

grid_data = load_solution_from_csv('results/Parabolic/data/ADI/grid_data.csv')

animate_time_series("results/Parabolic/data/ADI", "exact", grid_data, list(range(20, 2001, 20)), "Exact", save_path="results/Parabolic/exact.gif")
animate_time_series("results/Parabolic/data/ADI", "solution", grid_data, list(range(20, 2001, 20)), "Solution", save_path="results/Parabolic/ADI.gif")
animate_time_series_difference("results/Parabolic/data/ADI", "exact", "solution", grid_data, list(range(20, 2001, 20)), "Error", save_path="results/Parabolic/ADI_error.gif")