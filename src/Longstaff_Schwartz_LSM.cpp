#include <iostream>
#include <cmath>
#include <vector>
#include <cmath>
#include <numeric>
#include <Eigen/Dense>
#include <algorithm>
#include <omp.h>
#include "HestonModel.h"
using namespace std;
using namespace Eigen;

// Longstaff-Schwartz algorithm to calculate the American put option price
// The algorithm is based on the paper "Valuing American Options by Simulation: A Simple Least-Squares Approach" by Longstaff and Schwartz
// I use Heston model to simulate the stock price
// Use the poylnomial regression to estimate the continuation value

class Longstaff_Schwartz_LSM {
    public:
    vector<vector<double> > Paths(const vector<vector<double> >& variance_paths, int NoOfPaths, int NoOfSteps, double T,
     double r, double S_0, double kappa, double gamma, double rho, double vbar,
      double v0); // Generate the paths using Heston model
    double LSM(vector<vector<double> > &paths, double K,
     double r, double T);  // Calculate the option price using LSM


};

vector<vector<double> > Longstaff_Schwartz_LSM::Paths(const vector<vector<double> >& variance_paths, int NoOfPaths, int NoOfSteps, double T,
     double r, double S_0, double kappa,
    double gamma, double rho, double vbar, double v0) {
    Heston heston; 
    vector<vector<double> > paths = heston.GeneratePathsHestonAES(variance_paths, NoOfPaths, NoOfSteps, T,
     r, S_0, kappa,
    gamma, rho, vbar); // Generate the paths 
    return paths; // Return the paths
};

double Longstaff_Schwartz_LSM::LSM(vector<vector<double> > &paths, double K, double r, double T) {
    int NoOfPaths = paths.size();
    int NoOfSteps = paths[0].size();

    double dt = T / NoOfSteps; // Time step size for discounting
    vector<double> discountedPayoff(NoOfPaths, 0.0); // Final payoffs

    // Iterate backward in time to calculate continuation values
    #pragma omp parallel for
    for (int t = NoOfSteps - 2; t >= 0; --t) {
        // Step 1: Identify in-the-money paths
        vector<int> inTheMoneyIndices;
        for (int i = 0; i < NoOfPaths; ++i) {
            if (paths[i][t] > K) {  // In-the-money check
                inTheMoneyIndices.push_back(i);
            }
        }

        int inTheMoneyCount = inTheMoneyIndices.size();
        if (inTheMoneyCount == 0) {
            continue;  // Skip if no in-the-money paths
        }

        // Step 2: Build regression matrices X (independent) and Y (dependent)
        MatrixXd X(inTheMoneyCount, 3); // Quadratic polynomial basis
        VectorXd Y(inTheMoneyCount);    // Continuation values
        #pragma omp parallel for 
        for (int idx = 0; idx < inTheMoneyCount; ++idx) {
            int pathIdx = inTheMoneyIndices[idx];
            double S_t = paths[pathIdx][t];
            X(idx, 0) = 1.0;         // Constant term
            X(idx, 1) = S_t;         // Linear term
            X(idx, 2) = S_t * S_t;   // Quadratic term
            Y(idx) = discountedPayoff[pathIdx] * exp(-r * dt); // Discounted continuation value
        }

        // Step 3: Perform regression to calculate coefficients
        VectorXd coeffs = (X.transpose() * X).ldlt().solve(X.transpose() * Y);

        // Step 4: Calculate continuation values and decide payoff
        #pragma omp parallel for
        for (int idx = 0; idx < inTheMoneyCount; ++idx) {
            int pathIdx = inTheMoneyIndices[idx];
            double S_t = paths[pathIdx][t];
            double continuationValue = coeffs[0] + coeffs[1] * S_t + coeffs[2] * S_t * S_t;
            double immediatePayoff = max(paths[pathIdx][t] - K, 0.0);

            // Update discounted payoff if early exercise is optimal
            if (immediatePayoff > continuationValue) {
                discountedPayoff[pathIdx] = immediatePayoff * exp(-r * t * dt);
            }
        }
    }

    // Calculate the option price as the average payoff
    double optionPrice = accumulate(discountedPayoff.begin(), discountedPayoff.end(), 0.0) / NoOfPaths;
    return optionPrice;
}





    
