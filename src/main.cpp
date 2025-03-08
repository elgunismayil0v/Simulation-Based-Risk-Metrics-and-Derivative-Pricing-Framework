#include <iostream>
#include <chrono>  // For timing functions
#include <vector>
#include "VarianceGenerator.h"
#include "AssetPathGenerator.h"
#include "OptionPricer.h"
#include "EuropeanTypeVanillaOptionPricer.h"
#include "EuropeanTypeLookBackOptionPricer.h"
#include "AmericanOptionPricer.h"
#include "Exposure.h"
#include "ExceptedExposure.h"
#include "PotentialFutureExposure.h"
#include "Cva.h"
#include "Rho.h"
#include "Delta.h"
#include "Gamma.h"
#include "Vega.h"
#include "Theta.h"

using namespace std;
using namespace std::chrono;  // For timing functions

int main() {
    // Start timing
    auto start = high_resolution_clock::now();

    // Variance process parameters
    int NoOfPaths = 100;
    int NoOfSteps = 100;
    double kappa = 0.05;
    double gamma = 0.05;
    double vbar = 0.04;
    double v0 = 0.05;
    double T = 2.0;
    
    // Asset path parameters
    double S0 = 100.0;
    double rho = 0.7;
    double r = 0.05; 
    double strike = 100.0;

    // Instantiate a EuropeanPlainVanilla option (Call)
    LongstaffSchwartzLSM option(NoOfPaths, strike, T, r, LongstaffSchwartzLSM::Call,
                                S0, rho, NoOfSteps, kappa, gamma, vbar, v0);

    // Also instantiate a Put option if needed
    //EuropeanPlainVanilla option1(NoOfPaths, strike, T, r, EuropeanPlainVanilla::Call,
    //                             S0, rho, NoOfSteps, kappa, gamma, vbar, v0);

    // Price the option using the asset paths.
    double price = option();
    //double price1 = option1();

    // Create the exposure object (using 10 simulations at time 0.25)
    Exposure exposure(option, 10, 0.25);
    ExceptedExposure expected_exposure(exposure);
    PotentialFutureExposure pfe(exposure, 0.95);
    Cva cva(exposure, option, 0.03, 0.6);  // Ensure Cva constructor matches these parameters

    //Exposure exposure1(option1, 10, 0.25);
    //ExceptedExposure expected_exposure1(exposure1);
    //PotentialFutureExposure pfe1(exposure1, 0.95);
    //Cva cva1(exposure1, option1, 0.03, 0.6);  // Ensure Cva constructor matches these parameters

    // Compute risky derivative
    double risky_derivative = price - cva();  
    //double risky_derivative1 = price1 - cva1();  

    // Get the exposure vector and print every value.
    //vector<double> exposures = exposure();
   // cout << "Exposure values:" << endl;
    //for (size_t i = 0; i < exposures.size(); ++i) {
    //    cout << "Exposure[" << i << "] = " << exposures[i] << endl;
    //}

    // End timing
    auto end = high_resolution_clock::now();
    duration<double> elapsed = end - start;

    // Output the computed values.
    cout << "The option price (Call_A) is: " << price << endl;
    //cout << "The option price (Call) is: " << price1 << endl;
    cout << "The option EE (Call) is: " << expected_exposure() << endl;
    cout << "The option PFE (95%) (Call) is: " << pfe() << endl;
    cout << "The option CVA (Call) is: " << cva() << endl;  
    cout << "The option price CVA - option value (Call) is: " << risky_derivative << endl;
    //cout << "The option EE (Call) is: " << expected_exposure1() << endl;
    //cout << "The option PFE (95%) (Call) is: " << pfe1() << endl;
    //cout << "The option CVA (Call) is: " << cva1() << endl;  
    //cout << "The option price CVA - option value (Call) is: " << risky_derivative1 << endl;
    cout << "Simulation took " << elapsed.count() << " seconds." << endl;

    return 0;
}
