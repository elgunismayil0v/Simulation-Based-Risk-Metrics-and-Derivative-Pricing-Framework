import pandas as pd
import matplotlib.pyplot as plt

# Function to plot the data
def plot_simulation_data(data, title, y_label, paths_to_plot=None, horizontal_line=None):
    # Transpose the data: Rows (paths), Columns (time steps)
    transposed_data = data
    transposed_data.index.name = 'Path'  # Set the index name
    transposed_data.columns = [f"Step {i}" for i in range(transposed_data.shape[1])]
    
    # Select a subset of paths for clarity (if paths_to_plot is None, plot all)
    if paths_to_plot is None:
        paths_to_plot = transposed_data.index[:100]  # Adjust as needed
    transposed_data.loc[paths_to_plot].T.plot(figsize=(10, 6), title=title)
    
     # Add a horizontal line if specified
    if horizontal_line is not None:
        plt.axhline(y=horizontal_line, color='red', linestyle='--', linewidth=2, label=f"K = {horizontal_line}")

    # Customize plot
    plt.xlabel('Time Step')
    plt.ylabel(y_label)
    plt.legend().set_visible(False)
    plt.grid()
    plt.tight_layout()  # Adjust layout for better appearance
    plt.show()

# Load data
heston = pd.read_csv("/Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/plot/HestonPaths.csv")
#foward = pd.read_csv("/Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/plot/forwardRates.csv")
variance_paths = pd.read_csv("/Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/plot/variancepaths.csv")

# Plot for Heston data
plot_simulation_data(heston,
                     "Monte Carlo Simulation: Price Paths (Heston)",
                     "Price_of_European_Call", horizontal_line= 100)

# Plot for variance data
plot_simulation_data(variance_paths, "Monte Carlo Simulation: Variance Paths ","CIR process",horizontal_line= 0.0)







