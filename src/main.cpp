#include <iostream>
#include <chrono>  // Include the chrono library for timing
#include "VarianceGenerator.h"
#include "AssetPathGenerator.h"
#include "OptionPricer.h"
#include "EuropeanTypeAsianOptionPricer.h"
#include "EuropeanTypeLookBackOptionPricer.h"
#include <vector>

using namespace std;
using namespace std::chrono;  // For timing functions

int main() {
    // Start timing
    auto start = high_resolution_clock::now();

    // Variance process parameters
    int NoOfPaths = 10;
    int NoOfSteps = 1000;
    double T = 1.0;
    double kappa = 2.0;
    double gamma = 0.5;
    double vbar = 0.04;
    double v0 = 0.04;
    
    // Asset path parameters
    double S0 = 100.0;
    double rho = 0.8;
    double r = 0.05; 
    double strike = 100;

    // Create variance process object and generate variance paths.
    VarianceProcess varianceProcess(NoOfSteps, T, kappa, gamma, vbar, v0);
    vector<double> variancePaths = varianceProcess();

    // Create asset paths object with all necessary parameters.
    AssetPaths assetPathsObj(S0, rho, r, NoOfSteps, T, kappa, gamma, vbar);
    
    // Compute asset paths from the variance paths.
    vector<double> assetPathsResult = assetPathsObj(variancePaths);

    // Display variance process paths.
    cout << "Variance Process Sample Output:" << endl;
    for (double val : variancePaths) {
        cout << val << " ";
    }
    cout << endl;

    // Display asset paths.
    cout << "Asset Paths Sample Output:" << endl;
    for (double val : assetPathsResult) {
        cout << val << " ";
    }
    cout << endl;

    // Instantiate a EuropeanPlainVanilla option.
    AsianOption option(NoOfPaths, strike, T, r, AsianOption::Call, 
    AsianOption::Fixed, S0, rho, NoOfSteps, kappa, gamma, vbar, v0);


    AsianOption option1(NoOfPaths, strike,T, r, AsianOption::Put, AsianOption::Floating ,
    S0, rho, NoOfSteps, kappa, gamma, vbar, v0);

    // Price the option using the asset paths.
    double price = option();
    double price1 = option1();

    // End timing
    auto end = high_resolution_clock::now();
    duration<double> elapsed = end - start;

    // Output the computed price.
    cout << "The option price (Call) is: " << price << endl;
    cout << "The option price (Put) is: " << price1 << endl;

    // Display simulation time
    cout << "Simulation took " << elapsed.count() << " seconds." << endl;
    
    return 0;
}
