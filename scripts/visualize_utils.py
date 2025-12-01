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
from matplotlib.animation import FuncAnimation

def load_solution_from_csv(filename):
    """
    Loads a solution from a CSV file.

    Parameters:
    filename (str): Path to the CSV file.

    Returns:
    2D array: Loaded solution values.
    """
    return np.loadtxt(filename, delimiter=',')

def visualize_solution(data: np.ndarray, grid: np.ndarray, title="Data Visualization", save_path=None):
    """
    Visualizes the data on a 2D grid using a heatmap.

    Parameters:
    data (2D array): Data values at the grid points.
    grid (2D array): Grid informations at the grid points.
    title (str): Title of the plot.
    """
    nx, ny = data.shape
    x = np.linspace(0, 2, nx)
    y = np.linspace(-2, 2, ny)
    grid_x, grid_y = np.meshgrid(x, y)

    mask = (grid > 0)
    plot_dat = np.where(mask, data, np.nan).T

    plt.figure(figsize=(6, 8))
    plt.contourf(grid_x, grid_y, plot_dat, levels=50, cmap='seismic')
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

def animate_time_series(
    data_dir: str,
    prefix: str,
    grid: np.ndarray,
    steps: list,
    title: str,
    interval=100,
    save_path=None,
):
    """
    Generate animation from a time series of CSV solution snapshots.

    Parameters:
    - data_dir (str): Directory containing CSV files.
    - prefix (str): File prefix, e.g. "rhs" for rhs_000100.csv.
    - grid (2D array): Grid mask data (same as visualize_solution).
    - steps (list[int]): A list of time steps, e.g. range(100, 5001, 100)
    - title (str): Title of the plot.
    - interval (int): Delay between frames (ms).
    - save_path (str): Path to save animation; extension determines type.
                      e.g. "anim.mp4" or "anim.gif".
    """

    # prepare grid for contourf
    # nx, ny = grid.shape
    # x = np.linspace(0, 2, nx)
    # y = np.linspace(-2, 2, ny)
    # grid_x, grid_y = np.meshgrid(x, y)
    mask = (grid > 0)

    global_min, global_max = np.inf, -np.inf
    for s in steps:
        filename = f"{data_dir}/{prefix}_{s:06d}.csv"
        data = load_solution_from_csv(filename)
        masked = data[mask]
        local_min, local_max = masked.min(), masked.max()

        global_min = min(global_min, local_min)
        global_max = max(global_max, local_max)
    
    print(f"[Global] min={global_min}, max={global_max}")

    fig, ax = plt.subplots(figsize=(6, 8))
    # nx, ny = grid.shape
    extent = [0, 2, -2, 2]
    
    # initial frame
    first_file = f"{data_dir}/{prefix}_{steps[0]:06d}.csv"
    data0 = load_solution_from_csv(first_file)
    plot0 = np.where(mask, data0, np.nan)
    
    im = ax.imshow(plot0.T, extent=extent, origin='lower', 
                   cmap='seismic', aspect='equal',
                   vmin=global_min, vmax=global_max)
    
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(title)
    
    ax.set_title(f"t = {steps[0]}")
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")

    def update(frame_step):
        filename = f"{data_dir}/{prefix}_{frame_step:06d}.csv"
        data = load_solution_from_csv(filename)
        plot_data = np.where(mask, data, np.nan)
        
        im.set_data(plot_data.T)
        ax.set_title(f"t = {frame_step}")
        
        return im,

    ani = FuncAnimation(
        fig,
        update,
        frames=steps,
        interval=interval,
        blit=True,
        repeat=True,
    )

    plt.tight_layout()

    if save_path:
        print(f"Saving animation to {save_path} ...")
        if save_path.endswith(".gif"):
            ani.save(save_path, writer="pillow")
        else:
            ani.save(save_path, writer="ffmpeg")
        print("Animation saved.")

    return ani

if __name__ == "__main__":
    # Load the solution from a CSV file
    # exact = load_solution_from_csv('results/Parabolic/data/rhs_003300.csv')
    grid_data = load_solution_from_csv('results/Parabolic/data/grid_data.csv')

    # Visualize the solution
    # visualize_solution(exact, grid_data, title="Dirichlet Problem Exact Solution", save_path="results/Parabolic/rhs2.png")

    # Save Animation
    animate_time_series("results/Parabolic/data", "exact", grid_data, list(range(100, 5001, 100)), "Exact Solution", save_path="results/Parabolic/exact.gif")
    animate_time_series("results/Parabolic/data", "rhs", grid_data, list(range(100, 5001, 100)), "RHS Value", save_path="results/Parabolic/rhs.gif")