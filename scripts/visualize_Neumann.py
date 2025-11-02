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

def load_solution_from_csv(filename):
    """
    Loads a solution from a CSV file.

    Parameters:
    filename (str): Path to the CSV file.

    Returns:
    2D array: Loaded solution values.
    """
    return np.loadtxt(filename, delimiter=',')

def visualize_solution(grid_x, grid_y, solution, title="Solution Visualization", save_path=None,
                       overlay_mask=None, overlay_points=None, point_style=None):
    """
    Visualizes the solution on a 2D grid using a heatmap.

    Parameters:
    grid_x (2D array): X coordinates of the grid points.
    grid_y (2D array): Y coordinates of the grid points.
    solution (2D array): Solution values at the grid points.
    title (str): Title of the plot.
    """
    plt.figure(figsize=(6, 8))
    plt.contourf(grid_x, grid_y, solution, cmap='seismic', levels=np.linspace(-0.02, 0.02, 1000))
    plt.contourf(grid_x, grid_y, solution, levels=50, cmap='viridis')
    plt.colorbar(label='Solution Value')
    plt.title(title)
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    # Ensure equal scaling on both axes so grid cells appear square
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')

    # Overlay points: either provide a boolean mask (same shape as solution)
    # or a list/array of (x,y) coordinates in the same coordinate system as grid_x/grid_y.
    if overlay_mask is not None:
        # expect overlay_mask to be same shape as solution
        mask = np.asarray(overlay_mask)
        if mask.shape != solution.shape:
            raise ValueError("overlay_mask must have the same shape as solution")
        # scatter the points where mask is True
        xs = grid_x[mask]
        ys = grid_y[mask]
        ps = point_style or {'c': 'k', 's': 12, 'marker': 'o'}
        ax.scatter(xs, ys, **ps)
    elif overlay_points is not None:
        pts = np.asarray(overlay_points)
        if pts.ndim != 2 or pts.shape[1] != 2:
            raise ValueError("overlay_points must be an (N,2) array of (x,y) pairs")
        ps = point_style or {'c': 'k', 's': 12, 'marker': 'o'}
        ax.scatter(pts[:,0], pts[:,1], **ps)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

if __name__ == "__main__":
    # Load the solution from a CSV file
    solution = load_solution_from_csv('results/Neumann_solution.csv')
    exact = load_solution_from_csv('results/Neumann_exact.csv')

    # Create a grid for visualization (assuming uniform grid for simplicity)
    nx, ny = solution.shape
    x = np.linspace(0, 2, nx)
    y = np.linspace(-2, 2, ny)
    grid_x, grid_y = np.meshgrid(x, y)

    # Visualize the solution
    sol = solution.T
    ex = exact.T
    # Create a mask of data points (non-zero entries). If zeros are valid data,
    # prefer to supply an explicit overlay_mask or overlay_points instead.
    data_mask = sol != 0

    visualize_solution(grid_x, grid_y, sol, title="Neumann Problem Solution", save_path="results/Neumann_solution.png")
    visualize_solution(grid_x, grid_y, ex, title="Neumann Problem Exact Solution", save_path="results/Neumann_exact.png")
    visualize_solution(grid_x, grid_y, (ex - sol)**2, title="Neumann Problem Error", save_path="results/Neumann_error.png")

    # visualize_solution(grid_x, grid_y, sol, title="Neumann Problem Solution", save_path="results/Neumann_solution.png", overlay_mask=data_mask)
    # visualize_solution(grid_x, grid_y, ex, title="Neumann Problem Exact Solution", save_path="results/Neumann_exact.png", overlay_mask=(ex != 0))
    # visualize_solution(grid_x, grid_y, (ex - sol)**2, title="Neumann Problem Error", save_path="results/Neumann_error.png", overlay_mask=data_mask)