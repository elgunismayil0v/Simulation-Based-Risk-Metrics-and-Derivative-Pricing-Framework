#include "EuropeanTypeLookBackOptionPricer.h"
#include "AssetPathGenerator.h"
#include "VarianceGenerator.h"
#include <vector>
#include <cmath>  // for exp, sqrt, log
#include <numeric>  // for accumulate
#include <algorithm> // for max, sort
#include <omp.h>
#include <stdexcept>

LookBack:: LookBack(int NoOfPaths_, double strike_, double T_, double r_, OptionType OptionType_, StrikeType StrikeType_, double S0_, double rho_, int NoOfSteps_,
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

double LookBack:: operator()() const {
    double total_payoff = 0.0;
    for(int i = 0; i < NoOfPaths; i++) {
        // Define variance process
        VarianceProcess VarianceProcess(NoOfSteps, T, kappa, gamma, vbar, v0);
        vector<double> variance_process = VarianceProcess();
        // Define asset paths process
        AssetPaths AssetPaths(S0, rho, r, NoOfSteps, T, kappa, gamma, vbar);
        vector<double> asset_paths = AssetPaths(variance_process);
        // Define key inputs
        double S_T = asset_paths.back();
        double max_S_t = *max_element(asset_paths.begin(), asset_paths.end());
        double min_S_t = *min_element(asset_paths.begin(), asset_paths.end());

        switch(Option_Type){
            case Call:
                switch(Strike_Type){
                    case Fixed:
                        total_payoff += max(max_S_t - strike, 0.0) * exp(-r * T);
                        break;
                    case Floating:
                        total_payoff += max(S_T - min_S_t, 0.0) * exp(-r * T);
                        break;
                }
                break;
            case Put:
                switch(Strike_Type) {
                    case Fixed:
                        total_payoff += max(strike - min_S_t, 0.0) * exp(-r * T);
                        break;
                    case Floating:
                        total_payoff += max(max_S_t - S_T, 0.0) * exp(-r * T);
                        break;
                }
                break;
        }



    }

    return total_payoff / NoOfPaths;

}