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

def load_solution_from_csv(filename):
    """
    Loads a solution from a CSV file.

    Parameters:
    filename (str): Path to the CSV file.

    Returns:
    2D array: Loaded solution values.
    """
    return np.loadtxt(filename, delimiter=',')

def visualize_solution(grid_x, grid_y, solution, title="Solution Visualization", save_path=None):
    """
    Visualizes the solution on a 2D grid using a heatmap.

    Parameters:
    grid_x (2D array): X coordinates of the grid points.
    grid_y (2D array): Y coordinates of the grid points.
    solution (2D array): Solution values at the grid points.
    title (str): Title of the plot.
    """
    plt.figure(figsize=(6, 8))
    plt.contourf(grid_x, grid_y, solution, levels=50, cmap='viridis')
    plt.colorbar(label='Solution Value')
    plt.title(title)
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    # Ensure equal scaling on both axes so grid cells appear square
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

if __name__ == "__main__":
    # Load the solution from a CSV file
    solution = load_solution_from_csv('results/Dirichlet_solution.csv')
    exact = load_solution_from_csv('results/Dirichlet_exact.csv')

    # Create a grid for visualization (assuming uniform grid for simplicity)
    nx, ny = solution.shape
    x = np.linspace(0, 2, nx)
    y = np.linspace(-2, 2, ny)
    grid_x, grid_y = np.meshgrid(x, y)

    # Visualize the solution
    visualize_solution(grid_x, grid_y, solution.T, title="Dirichlet Problem Solution", save_path="results/Dirichlet_solution.png")
    visualize_solution(grid_x, grid_y, exact.T, title="Dirichlet Problem Exact Solution", save_path="results/Dirichlet_exact.png")
    visualize_solution(grid_x, grid_y, (exact.T - solution.T)**2, title="Dirichlet Problem Error", save_path="results/Dirichlet_error.png")