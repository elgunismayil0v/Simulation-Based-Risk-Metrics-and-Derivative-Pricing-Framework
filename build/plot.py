import pandas as pd
import matplotlib.pyplot as plt

# Function to plot the data
def plot_simulation_data(data, title, y_label, paths_to_plot=None):
    # Transpose the data: Rows (paths), Columns (time steps)
    transposed_data = data
    transposed_data.index.name = 'Path'  # Set the index name
    transposed_data.columns = [f"Step {i}" for i in range(transposed_data.shape[1])]
    
    # Select a subset of paths for clarity (if paths_to_plot is None, plot all)
    if paths_to_plot is None:
        paths_to_plot = transposed_data.index[:100]  # Adjust as needed
    
    transposed_data.loc[paths_to_plot].T.plot(figsize=(10, 6), title=title)

    # Customize plot
    plt.xlabel('Time Step')
    plt.ylabel(y_label)
    plt.legend().set_visible(False)
    plt.grid()
    plt.tight_layout()  # Adjust layout for better appearance
    plt.show()

# Load data
heston = pd.read_csv("/Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/plot/HestonPaths.csv")
foward = pd.read_csv("/Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/plot/forwardRates.csv")
lsmc = pd.read_csv("/Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/plot/LSMPaths.csv")

# Plot for Heston data
plot_simulation_data(heston, "Monte Carlo Simulation: Price Paths (Heston)", "Price_of_European_Call")

# Plot for LSMC data
plot_simulation_data(lsmc, "Monte Carlo Simulation: Price Paths (LSMC)", "Price_of_American_Call")

# Plot for Forward data
plot_simulation_data(foward, "Monte Carlo Simulation: Price Paths (Forward)", "Price_of_forward_rate")


