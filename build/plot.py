import pandas as pd
import matplotlib.pyplot as plt

heston = pd.read_csv( "/Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/plot/HestonPaths.csv")
foward = pd.read_csv( "/Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/plot/forwardRates.csv")
lsmc = pd.read_csv( "/Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/plot/LSMPaths.csv")

data_1 = heston

# Transpose the data: Rows (paths), Columns (time steps)
transposed_data = data_1.T
transposed_data.index.name = 'Path'  # Set the index name
transposed_data.columns = [f"Step {i}" for i in range(transposed_data.shape[1])]

# Plot: Select a subset of paths for clarity (e.g., first 5 paths)
paths_to_plot = transposed_data.index[:5]  # Adjust as needed
transposed_data.loc[paths_to_plot].T.plot(figsize=(10, 6), title="Monte Carlo Simulation: Price Paths")

# Customize plot
plt.xlabel('Time Step')
plt.ylabel('Price_of_European_Call')
plt.legend(title='Paths', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid()
plt.tight_layout()  # Adjust layout for better appearance
plt.show()


data_2 = lsmc

# Transpose the data: Rows (paths), Columns (time steps)
transposed_data = data_2.T
transposed_data.index.name = 'Path'  # Set the index name
transposed_data.columns = [f"Step {i}" for i in range(transposed_data.shape[1])]

# Plot: Select a subset of paths for clarity (e.g., first 5 paths)
paths_to_plot = transposed_data.index[:5]  # Adjust as needed
transposed_data.loc[paths_to_plot].T.plot(figsize=(10, 6), title="Monte Carlo Simulation: Price Paths")

# Customize plot
plt.xlabel('Time Step')
plt.ylabel('Price_of_American_Call')
plt.legend(title='Paths', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid()
plt.tight_layout()  # Adjust layout for better appearance
plt.show()


data_3 = foward

# Transpose the data: Rows (paths), Columns (time steps)
transposed_data = data_3.T
transposed_data.index.name = 'Path'  # Set the index name
transposed_data.columns = [f"Step {i}" for i in range(transposed_data.shape[1])]

# Plot: Select a subset of paths for clarity (e.g., first 5 paths)
paths_to_plot = transposed_data.index[:5]  # Adjust as needed
transposed_data.loc[paths_to_plot].T.plot(figsize=(10, 6), title="Monte Carlo Simulation: Price Paths")

# Customize plot
plt.xlabel('Time Step')
plt.ylabel('Price_of_forward_rate')
plt.legend(title='Paths', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid()
plt.tight_layout()  # Adjust layout for better appearance
plt.show()

