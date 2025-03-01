#include "EuropeanTypeAsianOptionPricer.h"
#include "AssetPathGenerator.h"
#include "VarianceGenerator.h"
#include <vector>
#include <cmath>  // for exp, sqrt, log
#include <numeric>  // for accumulate
#include <algorithm> // for max, sort
#include <omp.h>
#include <stdexcept>


AsianOption::AsianOption(int NoOfPaths_, double strike_, double T_, double r_, OptionType OptionType_, StrikeType StrikeType_, double S0_, double rho_, int NoOfSteps_,
               double kappa_, double gamma_, double vbar_, double v0_) {
    NoOfPaths = NoOfPaths_;
    strike = strike_;
    T = T_;
    r = r_;
    Option_Type = OptionType_;
    Strike_Type = StrikeType_;
    S0 = S0_;
    rho = rho_;
    NoOfSteps = NoOfSteps_;
    kappa = kappa_;
    gamma= gamma_;
    vbar = vbar_;
    v0 = v0_;
}

double AsianOption::operator()() const {
    double total_payoff = 0.0;

    for(int i = 0; i < NoOfPaths; i++) {
        // Define variance process
        VarianceProcess VarianceProcess(NoOfSteps, T, kappa, gamma, vbar, v0);
        vector<double> variance_process = VarianceProcess();
        // Define asset paths
        AssetPaths AssetPaths(S0, rho, r, NoOfSteps, T, kappa, gamma, vbar);
        vector<double> asset_paths = AssetPaths(variance_process);
        // Define the key parametrs 
        double S_T = asset_paths.back();
        double mean_of_path = accumulate(asset_paths.begin(), asset_paths.end(), 0.0) / NoOfSteps; 

        switch(Option_Type) {
            case Call:
            switch(Strike_Type) {
                case Fixed:
                total_payoff += max(mean_of_path - strike, 0.0) * exp(-r * T);
                break;
                case Floating:
                total_payoff += max(S_T - mean_of_path, 0.0) * exp(-r * T);
                break;
            }
            break;
            case Put:
            switch(Strike_Type) {
                case Fixed:
                total_payoff += max(strike - mean_of_path, 0.0) * exp(-r * T);
                break;
                case Floating:
                total_payoff += max(mean_of_path - S_T, 0.0) * exp(-r * T);
                break;
            }
            break;

        }

    }

    return total_payoff / NoOfPaths;
    
   
}


