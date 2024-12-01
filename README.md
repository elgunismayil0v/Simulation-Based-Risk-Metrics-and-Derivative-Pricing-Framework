
# Heston Model with Monte Carlo Simulation

This project implements the **Heston stochastic volatility model** and **Monte Carlo simulation** for pricing European call options. It also provides calculations for **Expected Exposure (EE)**, **Potential Future Exposure (PFE)**, and **Credit Valuation Adjustment (CVA)** using simulated asset price paths.

## Features

- **Heston Model Simulation**:
  - Simulates asset price paths using the Heston model, incorporating stochastic volatility.
  - Implements the Cox-Ingersoll-Ross (CIR) process for variance dynamics.
- **Monte Carlo Pricing**:
  - Computes European call option prices via Monte Carlo methods.
- **Risk Metrics Calculations**:
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
   g++ -std=c++11 Heston_Model_with_Monte_Carlo_simulation.cpp -o HestonModel
   ```
3. Run the executable:
   ```bash
   ./HestonModel
   ```

## Input Parameters

Key parameters used in the simulation (modifiable within the code):

- **Simulation**:
  - Number of paths: 1000
  - Time steps: 1000
  - Maturity (\( T \)): 1 year
- **Market Data**:
  - Risk-free rate (\( r \)): 0.05
  - Initial stock price (\( S_0 \)): 100
  - Strike price (\( K \)): 100
- **Heston Model Parameters**:
  - Mean reversion rate (\( \kappa \)): 2.0
  - Volatility of variance (\( \gamma \)): 0.8
  - Correlation (\( 
ho \)): -0.9
  - Long-term variance mean (\( ar{v} \)): 0.66
  - Initial variance (\( v_0 \)): 0.04
- **Risk Metrics Parameters**:
  - Confidence level for PFE: 95%
  - Recovery rate: 40%
  - Hazard rate: 2%

## Output

After execution, the program provides the following metrics:

- **Mean Expected Exposure (EE)**
- **Mean Potential Future Exposure (PFE)**
- **Credit Valuation Adjustment (CVA)**
- **European Call Option Price**

Sample output:

```
Mean Expected Exposure (EE): 4.5678
Mean Potential Future Exposure (PFE): 6.1234
Credit Valuation Adjustment (CVA): 0.7890
European Call Option Price: 10.3456
```

## Project Structure

- `Heston` Class:
  - Implements variance simulation using the CIR model.
  - Simulates asset price paths under the Heston model.
  - Calculates European call option prices.
- `RiskMetrics` Class:
  - Computes risk measures (EE, PFE, CVA) from simulated paths.

