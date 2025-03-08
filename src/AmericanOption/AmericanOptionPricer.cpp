#include "AmericanOptionPricer.h"
#include "VarianceGenerator.h"
#include "AssetPathGenerator.h"
#include <vector>  // for vector
#include <cmath>  // for exp, sqrt, log
#include <random> // for random_device, mt19937, normal_distribution
#include <numeric>  // for accumulate
#include <algorithm> // for max, sort
#include <omp.h>  // for OpenMP
#include <Eigen/Dense> // for matrix and regression operations
using namespace std;
using namespace Eigen;

LongstaffSchwartzLSM::LongstaffSchwartzLSM(int NoOfPaths_ ,double strike_, double T_, double r_, OptionType OptionType_, double S0_, double rho_, int NoOfSteps_,
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
}

vector<vector<double>> LongstaffSchwartzLSM::matrix() const {
     vector<vector<double>> paths_matrix(NoOfPaths, vector<double>(NoOfSteps));;

    for(int i = 0; i < NoOfPaths; i++) {
        // Create objects properly
        AssetPaths assetPaths(S0, rho, r, NoOfSteps, T, kappa, gamma, vbar);
        VarianceProcess VarianceProcess(NoOfSteps, T, kappa, gamma, vbar, v0); // Assuming VarianceProcess has a method to generate variance paths

        // Generate a variance path first
        vector<double> variancePath = VarianceProcess();  

        // Call operator() on the AssetPaths object
        vector<double> paths = assetPaths(variancePath); 
        paths_matrix[i] = paths;

    }

    return paths_matrix;

}

vector<int> LongstaffSchwartzLSM::ITM_Paths(const vector<vector<double>>& paths_matrix, int t) const {
    vector<int> itm_paths_vector; 
        for(int i = 0; i < NoOfPaths; i++){
            double S_t = paths_matrix[i][t];
            if ((Option_Type == Call && S_t > strike) ||
            (Option_Type == Put && S_t < strike)) {
            itm_paths_vector.push_back(i);
            }
        }

    return itm_paths_vector;
}






vector<double> LongstaffSchwartzLSM::payoff(const vector<vector<double>>& paths_matrix) const { 
    double dt = T / NoOfSteps; // Time step size for discounting
    vector<double> discounted_payoff(NoOfPaths, 0.0); // Final payoffs
    double immediate_payoff = 0.0;
    for(int t = NoOfSteps -1; t >= 0; t--){
        vector<int> itm_paths_vector = ITM_Paths(paths_matrix, t);
        MatrixXd X(itm_paths_vector.size(), 3);
        VectorXd y(itm_paths_vector.size());
        // If no paths are ITM, skip regression
        if (itm_paths_vector.empty()) continue;

        for(int i = 0; i < itm_paths_vector.size(); i++) {
            int path_index = itm_paths_vector[i];  // Correct path index
            double S_t = paths_matrix[path_index][t];
            X(i, 0) =1.0;
            X(i, 1) = S_t;
            X(i, 2) = S_t * S_t;
            y(i) = discounted_payoff[path_index] * exp(-r * dt);

                
        }
        VectorXd coeffs = (X.transpose() * X).ldlt().solve(X.transpose() * y);
        for(int i = 0; i < itm_paths_vector.size(); i++) {
            int path_index = itm_paths_vector[i];
            double S_t = paths_matrix[path_index][t];
            double continuation_value = coeffs[0] + coeffs[1] * S_t + coeffs[2] * S_t * S_t; 
            double immediate_payoff = (Option_Type == Call) ? max(S_t - strike, 0.0) 
                                                        : max(strike - S_t, 0.0);

                    // Check early exercise
            if (immediate_payoff > continuation_value) {
            discounted_payoff[path_index] = immediate_payoff;
            } else {
            discounted_payoff[path_index] *= exp(-r * dt);
            }

        }
    }
    return discounted_payoff;
}

double LongstaffSchwartzLSM::operator()() const {
    vector<vector<double>> paths_matrix = matrix();
    vector<double> Payoff = payoff(paths_matrix);
    return accumulate(Payoff.begin(), Payoff.end(), 0.0) / NoOfPaths;

}

