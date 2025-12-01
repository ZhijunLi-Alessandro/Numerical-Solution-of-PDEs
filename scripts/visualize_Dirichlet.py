## @file visualize_Dirichlet.py
#  @brief A python script for visualizing the calculation results of the Dirichlet problem
#  
#  This script is used to process and analyze data and generate figures.
#  Includes data loading and preprocessing, generate multiple figures.
#  @see Dirichlet.c
#  @author Li Zhijun
#  @date 2025-10-27
import numpy as np
import matplotlib.pyplot as plt
from visualize_utils import load_solution_from_csv, visualize_solution

# Load the solution from a CSV file
solution = load_solution_from_csv('results/Poisson/data/Dirichlet_solution.csv')
exact = load_solution_from_csv('results/Poisson/data/Dirichlet_exact.csv')
grid_data = load_solution_from_csv('results/Poisson/data/grid_data.csv')

# Visualize the solution
visualize_solution(solution, grid_data, title="Dirichlet Problem Solution", save_path="results/Poisson/Dirichlet_solution.png")
visualize_solution(exact, grid_data, title="Dirichlet Problem Exact Solution", save_path="results/Poisson/Dirichlet_exact.png")
visualize_solution((exact - solution)**2, grid_data, title="Dirichlet Problem Error", save_path="results/Poisson/Dirichlet_error.png")
