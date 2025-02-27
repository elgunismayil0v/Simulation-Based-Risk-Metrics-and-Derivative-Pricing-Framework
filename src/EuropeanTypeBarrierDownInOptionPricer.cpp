#include "EuropeanTypeBarrierDownInOptionPricer.h"
#include "AssetPathGenerator.h"
#include "VarianceGenerator.h"
#include <vector>
#include <cmath>  // for exp, sqrt, log
#include <numeric>  // for accumulate
#include <algorithm> // for max, sort
#include <omp.h>
#include <stdexcept>

BarrierKnockDownIn::BarrierKnockDownIn(double NoOfPaths_ ,double strike_, double Barrier_,double T_, double r_, OptionType OptionType_, double S0_, double rho_, int NoOfSteps_,
               double kappa_, double gamma_, double vbar_, double v0_) {
    NoOfPaths = NoOfPaths_;
    strike = strike_;
    T = T_;
    r = r_;
    Option_Type = OptionType_;
    S0 = S0_;
    rho = rho_;
    NoOfSteps = NoOfSteps_;
    kappa = kappa_;
    gamma= gamma_;
    vbar = vbar_;
    v0 = v0_;
    Barrier = Barrier_;

}

double BarrierKnockDownIn:: operator()() const {
    double total_payoff = 0.0;
    // Parallel loop with reduction for sum of payoffs
    #pragma omp parallel for reduction(+:payoff_call, payoff_put)
    for(int i = 0; i < NoOfPaths; i++) {
        VarianceProcess VarianceProcess(NoOfSteps, T, kappa, gamma, vbar, v0);
        AssetPaths AssetPaths(S0, rho, r, NoOfSteps, T, kappa, gamma, vbar);
        // Generate a variance path first
        vector<double> variancePath = VarianceProcess();  

        // Call operator() on the AssetPaths object
        vector<double> assetPaths = AssetPaths(variancePath); 

        bool knocked_in = false;
        for(int j = 0; j < NoOfSteps; j++) {
            if (assetPaths[j] <= Barrier) {
                knocked_in = true;
                break;  // Once it hits the barrier, stop checking further
            }
        }

        if (knocked_in) {
            switch(Option_Type) {
                case Call:
                total_payoff += max(assetPaths.back() - strike, 0.0) * exp(-r * T);
                break;
                case Put:
                total_payoff += max(strike - assetPaths.back(), 0.0) * exp(-r * T);
                break;
                default:
                throw std::invalid_argument("Invalid option type");
                
            }
        }


    }

    return total_payoff / NoOfPaths;

}