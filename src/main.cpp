#include "HestonModel.h" // Include HestonModel header file
#include "lmm.h" // Include LMM header file
#include "Longstaff_Schwartz.h" // Include Longstaff-Schwartz header file
#include "csv.h" // Include CSV header file
#include <iostream> // Use for standard input/output
#include <vector> // Use for vectors
#include <numeric> // Use for accumulate function
#include <chrono> // Use for time measurement
using namespace std;



int main() {
    auto start_time = chrono::high_resolution_clock::now();

    Heston Heston;
    RiskMetrics riskMetrics;
    // Calculate call price
    double kappa = 0.05;
    double gamma = 0.05;
    double vbar = 0.04;
    double v0 = 0.05;
    double r = 0.05;
    double T = 1.0;
    int NoOfPaths = 1000;
    int NoOfSteps = 1000;
    double S_0 = 100.0;
    double K_option = 100.0;
    double rho = 0.7;
    vector<vector<double> > variance_path = Heston.CIR_sample(NoOfPaths,NoOfSteps,T,kappa,gamma,vbar,v0);
    vector< vector < double > > paths = Heston.GeneratePathsHestonAES(variance_path, NoOfPaths, NoOfSteps, T, r, S_0, kappa, gamma, rho, vbar);
    double europeanCallPrice = Heston.CalculateEuropeanCallPrice(paths, K_option, r, T);
    cout << "European call price: " << europeanCallPrice << endl;

    // Calculate RiskMetrics of the Heston model
    vector<double> EE_h = riskMetrics.CalculateDiscountedExpectedExposure(paths, K_option, r, 1.0 / NoOfSteps);
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
    double optionValues = LSM.LSM(paths, K_option, r, T);
    cout << "Longstaff-Schwartz price: " << optionValues << endl;
    // Calculate RiskMetrics of the Longstaff-Schwartz model
    vector<double> EE_lsm = riskMetrics.CalculateDiscountedExpectedExposure(paths, K_option, r, 1.0 / NoOfSteps);
    double mean_EE_lsm = accumulate(EE_lsm.begin(), EE_lsm.end(), 0.0) / EE_lsm.size();
    cout << "Discounted Expected Exposure: " << mean_EE_lsm << endl;
    vector<double> PFE_lsm = riskMetrics.CalculatePotentialFutureExposure(paths, 0.95, K_option);
    double mean_PFE_lsm = accumulate(PFE_lsm.begin(), PFE_lsm.end(), 0.0) / PFE_lsm.size();
    cout << "Potential Future Exposure: " << mean_PFE_lsm << endl;
    double CVA_lsm = riskMetrics.CalculateCVA(EE_lsm, 0.4, 0.02, T, NoOfSteps);
    cout << "CVA: " << CVA_lsm << endl;
    cout << "Risky Derivative Price of Longstaff-Schwartz : " << optionValues - CVA_lsm << endl;

    

    // Csv files for the results
    writeToCSV(paths, "/Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/plot/HestonPaths.csv");
    writeToCSV(variance_path, "/Users/elgun/Desktop/Simulation-Based-Risk-Metrics-and-Option-Pricing-Frameworw/plot/variancepaths.csv");


     // Call Python script to plot from all files
    std::cout << "Generating Plot for All Files...\n";
    int result = system("python3 plot.py");
    if (result != 0) {
        std::cerr << "Error: Failed to run Python script.\n";
    }

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    cout << "Elapsed time: " << elapsed.count() << " s" << endl;




    return 0;
}