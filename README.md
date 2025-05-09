
# Heston and LMM_DD Models with Monte Carlo Simulation

This project implements the **Heston stochastic volatility model** for financial derivatives pricing using Monte Carlo simulation. It also provides risk metrics calculations such as **Expected Exposure (EE)**, **Potential Future Exposure (PFE)**, and **Credit Valuation Adjustment (CVA)** for both models.

## Features

### Heston Model
- Simulates asset price paths using the Heston model with stochastic volatility.
- Implements the Cox-Ingersoll-Ross (CIR) process for variance dynamics.
- Calculates European call option prices via Monte Carlo methods.


### Risk Metrics Calculations
- **Expected Exposure (EE)**: Average exposure over time.
- **Potential Future Exposure (PFE)**: Exposure at a specified confidence level.
- **Credit Valuation Adjustment (CVA)**: Measures counterparty credit risk.

## Prerequisites

- A C++ compiler supporting C++11 or later.
- Standard Template Library (STL).

## Compilation and Execution

1. Clone this repository to your local machine:
   ```bash
   git clone <repository-url>
   cd <repository-folder>
   ```
2. Compile the program:
   ```bash
   g++ *.cpp -o out
   ```
3. Run the executable:
   ```bash
   ./out
   ```

## Input Parameters

### Heston Model Parameters
- Mean reversion rate (\( \kappa \)): 0.005
- Volatility of variance (\( \gamma \)): 0.05
- Correlation (\( \rho \)): -0.7
- Long-term variance mean (\( \bar{v} \)): 0.04
- Initial variance (\( v_0 \)): 0.04
- Risk-free rate (\( r \)): 0.05
- Initial stock price (\( S_0 \)): 100
- Strike price (\( K \)): 100

### LMM_DD Model Parameters
- Mean reversion rate (\( \kappa \)): 0.005
- Volatility of variance (\( \gamma \)): 0.05
- Correlation (\( \rho \)): 0.0
- Long-term variance mean (\( \bar{v} \)): 0.04
- Initial variance (\( v_0 \)): 0.04
- Risk-free rate (\( r \)): 0.0
- Initial forward rate (\( f_0 \)): 0.1
- Strike rate ($K$): 0.03
- Number of paths: 1000
- Time steps: 252

## Output

The program provides the following metrics for both models:

- **European Call Option Price** (Heston)
- **Discounted Expected Exposure (EE)**
- **Potential Future Exposure (PFE)**
- **Credit Valuation Adjustment (CVA)**
- **Risky Derivative Prices**:
  - European Call Price minus CVA

### Sample Output
```
European call price: 10.1234
Discounted Expected Exposure: 4.5678
Potential Future Exposure: 6.1234
CVA: 0.7890
Risky Derivative Price of European Call Option: 9.3344
```

## Project Structure

- **Heston Class**:
  - Implements variance simulation using the CIR model.
  - Simulates asset price paths under the Heston model.
  - Calculates European call option prices.
- **RiskMetrics Class**:
  - Computes risk measures (EE, PFE, CVA) from simulated paths.

## License

This project is licensed under the MIT License.
