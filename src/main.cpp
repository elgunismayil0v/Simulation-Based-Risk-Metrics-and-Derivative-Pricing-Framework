#include "HestonModel.h"
#include "lmm.h"
#include "Longstaff_Schwartz.h"
#include "csv.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <chrono>
#include "plot.h"
using namespace std;



int main() {
    auto start_time = chrono::high_resolution_clock::now();

    Heston Heston;
    RiskMetrics riskMetrics;
    LMM_DD lmm_dd;

    // LMM_DD model parameters
    double kappa = 0.005;
    double gamma = 0.05;
    double rho = 0.0;
    double vbar = 0.04;
    double v0 = 0.04;
    double r = 0.0;
    double f_0 = 0.1;
    double T = 1.0;
    int NoOfPaths = 1000;
    int NoOfSteps = 1000;
    double beta = 0.5;
    double sigma = 0.1;

    // LMM parameters
    double K = 0.03;

    // Generate forward rates using the Heston model
    vector<vector<double> > forwardRates = lmm_dd.Forward_rates(NoOfPaths, NoOfSteps, T, r, f_0, kappa, gamma, rho, vbar, v0, beta, sigma);

    // Calculate caplet price
    double capletPrice = lmm_dd.Caplet_price(forwardRates, K);
    cout << "Caplet price: " << capletPrice << endl;

    // Calculate swap option price
    double swapOptionPrice = lmm_dd.Swapoption_price(forwardRates.back(), K);
    cout << "Swap option price: " << swapOptionPrice << endl;

    // Calculate RiskMetrics of the LMM_DD model
    vector<double> EE = riskMetrics.CalculateDiscountedExpectedExposureWithStrike(forwardRates, K, r, 1.0 / NoOfSteps);
    double mean_EE = accumulate(EE.begin(), EE.end(), 0.0) / EE.size();
    cout << "Discounted Expected Exposure: " << mean_EE << endl;
    vector<double> PFE = riskMetrics.CalculatePotentialFutureExposure(forwardRates, 0.95, K);
    double mean_PFE = accumulate(PFE.begin(), PFE.end(), 0.0) / PFE.size();
    cout << "Potential Future Exposure: " << mean_PFE << endl;
    double CVA = riskMetrics.CalculateCVA(EE, 0.4, 0.02, T, NoOfSteps);
    cout << "CVA: " << CVA << endl;
    cout << "Risky Derivative Price of Caplet : " << capletPrice - CVA << endl;
    cout << "Risky Derivative Price of Swap Option : " << swapOptionPrice - CVA << endl;

    // Calculate European call price
    double S_0 = 100.0;
    double K_option = 100.0;
    double r_option = 0.05;
    double rho_option = - 0.7;
    vector< vector < double > > paths = Heston.GeneratePathsHestonAES(NoOfPaths, NoOfSteps, T, r_option, S_0, kappa, gamma, rho_option, vbar, v0);
    double europeanCallPrice = Heston.CalculateEuropeanCallPrice(paths, K_option, r, T);
    cout << "European call price: " << europeanCallPrice << endl;

    // Calculate RiskMetrics of the Heston model
    vector<double> EE_h = riskMetrics.CalculateDiscountedExpectedExposureWithStrike(paths, K_option, r, 1.0 / NoOfSteps);
    double mean_EE_h = accumulate(EE_h.begin(), EE_h.end(), 0.0) / EE_h.size();
    cout << "Discounted Expected Exposure: " << mean_EE_h << endl;
    vector<double> PFE_h = riskMetrics.CalculatePotentialFutureExposure(paths, 0.95, K_option);
    double mean_PFE_h = accumulate(PFE_h.begin(), PFE_h.end(), 0.0) / PFE_h.size();
    cout << "Potential Future Exposure: " << mean_PFE_h << endl;
    double CVA_h = riskMetrics.CalculateCVA(EE_h, 0.4, 0.02, T, NoOfSteps);
    cout << "CVA: " << CVA_h << endl;
    cout << "Risky Derivative Price of European Call Option : " << europeanCallPrice - CVA_h << endl;

    // Longstaff-Schwartz algorithm
    Longstaff_Schwartz_LSM LSM;
    vector<vector<double> > paths_lsm = LSM.Paths(NoOfPaths, NoOfSteps, T, r_option, S_0, kappa, gamma, rho, vbar, v0);
    double optionValues = LSM.LSM(paths_lsm, K_option, r, T);
    cout << "Longstaff-Schwartz price: " << optionValues << endl;
    // Calculate RiskMetrics of the Longstaff-Schwartz model
    vector<double> EE_lsm = riskMetrics.CalculateDiscountedExpectedExposureWithStrike(paths_lsm, K_option, r, 1.0 / NoOfSteps);
    double mean_EE_lsm = accumulate(EE_lsm.begin(), EE_lsm.end(), 0.0) / EE_lsm.size();
    cout << "Discounted Expected Exposure: " << mean_EE_lsm << endl;
    vector<double> PFE_lsm = riskMetrics.CalculatePotentialFutureExposure(paths_lsm, 0.95, K_option);
    double mean_PFE_lsm = accumulate(PFE_lsm.begin(), PFE_lsm.end(), 0.0) / PFE_lsm.size();
    cout << "Potential Future Exposure: " << mean_PFE_lsm << endl;
    double CVA_lsm = riskMetrics.CalculateCVA(EE_lsm, 0.4, 0.02, T, NoOfSteps);
    cout << "CVA: " << CVA_lsm << endl;
    cout << "Risky Derivative Price of Longstaff-Schwartz : " << optionValues - CVA_lsm << endl;

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    cout << "Elapsed time: " << elapsed.count() << " s" << endl;

    // Csv files for the results
    appendToCSV(forwardRates, "/Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/plot/forwardRates.csv");
    appendToCSV(paths, "/Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/plot/HestonPaths.csv");
    appendToCSV(paths_lsm, "/Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/plot/LSMPaths.csv");

     // Call Python script to plot from all files
    std::cout << "Generating Plot for All Files...\n";
    int result = system("python3 plot.py");
    if (result != 0) {
        std::cerr << "Error: Failed to run Python script.\n";
    }




    return 0;
}