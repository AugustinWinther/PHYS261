"""Functions for visualizing N vs L 

"""
# Standard library
import argparse

# Third party
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Local
import calc

def plot_N_vs_L(K: int, N_max: int, L_max: int, PREC: float, 
                dL: float, max_err_percent: float) -> None:
    """Plots N vs L plot where color indicates where the first K numerical 
    energies all have under max_err_percent error relative to analytical 
    energies.

    PREC is how many points to use in numerical integration. Higher is better.
    """
    # Set all value combinations to test
    L_values = np.arange(1, L_max + dL, dL)
    N_values = np.arange(K, N_max + 1)

    # Create value grids
    L_gird, N_grid = np.meshgrid(L_values, N_values)

    # Calculate actuall energies (used for finding error)
    actual_energies = calc.analytical_energies(K)

    # Create grid to be plotted
    mean_error_grid = np.zeros((len(N_values), len(L_values)))
    
    for i, N in enumerate(N_values):
        for j, L in enumerate(L_values):
            approx_energies = calc.numerical_energies(K, N, L, PREC)
            energy_errors = np.abs(approx_energies/actual_energies - 1)*100
            
            if np.all(energy_errors < max_err_percent):
                mean_error_grid[i, j] = np.mean(energy_errors)
            else:
                mean_error_grid[i, j] = np.nan  # Not interested in values with
                                           # too high errors.

    # Create custom color map where "white" is used for all np.nan values
    cmap = plt.get_cmap("winter")
    cmap.set_bad(color="white")

    # Draw color plot
    color = plt.pcolormesh(L_gird, N_grid, mean_error_grid, cmap=cmap)
    plt.colorbar(color, label="Mean error %")

    # Set ticks and titles
    plt.xlabel("L")
    plt.ylabel("N", rotation=0)
    plt.title(f"N and L combinations that produce <{max_err_percent}% error\n"
              f"for first {K} energy states")

    # Force y-axsis (N) to use integer values
    plt.gca().yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    # Draw custom grid lines (plt.grid() does not support offsetting from ticks)
    grid_hlines = N_values - 0.5
    plt.gca().hlines(y=grid_hlines, xmin=1, xmax=L_max, 
                     color="Gray", linewidth=0.5)
    
    # Show plot
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="N vs. L plot for H.O. approximation",
                                     description=("Outputs a plot showing "
                                                  "which combinations of N "
                                                  "and L gives the best "
                                                  "numerical approximations."))

    parser.add_argument("-K", help="K first energies of H.O. to approximate.",
                        dest="K", 
                        type=int,
                        required=True)

    parser.add_argument("-N", help=("Max N value to test for."),
                        dest="N_max", 
                        type=int,
                        required=True)
    
    parser.add_argument("-L", help=("Max L value to test for."),
                        dest="L_max", 
                        type=float,
                        required=True)
    
    parser.add_argument('-P', help=("Amount of points used in numerical "
                                    "integration."),
                        dest="PREC", 
                        type=float, 
                        required=True)

    parser.add_argument('-E', help=("Max error of numerical approximation "
                                    "relative to analytical to allow."),
                        dest="max_err_percent", 
                        type=float,
                        required=True)

    parser.add_argument('-dL', help=("Difference between each L value."),
                        dest="dL", 
                        type=float,
                        required=True)

    args = parser.parse_args()

    plot_N_vs_L(args.K, args.N_max, args.L_max, args.PREC, args.dL, 
                args.max_err_percent)