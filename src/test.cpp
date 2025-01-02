#include <iostream>
#include <cmath>
#include <vector>
#include <cmath>
#include <numeric>
#include <Eigen/Dense>
#include "HestonModel.h"
using namespace std;
using namespace Eigen;

// Longstaff-Schwartz algorithm to calculate the American put option price
// The algorithm is based on the paper "Valuing American Options by Simulation: A Simple Least-Squares Approach" by Longstaff and Schwartz
// I use Heston model to simulate the stock price
// Use the poylnomial regression to estimate the continuation value

class Longstaff_Schwartz_LSM {
    public:
    vector<vector<double> > Paths(int NoOfPaths, int NoOfSteps, double T, double r, double S_0, double kappa, double gamma, double rho, double vbar, double v0); // Generate stock price paths
    double LSM(vector<vector<double> > &paths, double K, double r, double T);  // Longstaff-Schwartz algorithm


};

vector<vector<double> > Longstaff_Schwartz_LSM::Paths(int NoOfPaths, int NoOfSteps, double T,
     double r, double S_0, double kappa,
    double gamma, double rho, double vbar, double v0) {
    Heston heston;
    vector<vector<double> > paths = heston.GeneratePathsHestonAES(NoOfPaths, NoOfSteps, T,
     r, S_0, kappa,
    gamma, rho, vbar, v0);
    return paths;
};

double Longstaff_Schwartz_LSM::LSM(vector<vector<double> > &paths, double K, double r, double T) {
    int NoOfPaths = paths.size();
    int NoOfSteps = paths[0].size();

    VectorXd regressionCoefficients;
    VectorXd continuationValue(NoOfPaths);
    MatrixXd X(NoOfPaths, 3);
    VectorXd Y(NoOfPaths);

    double dt = T / NoOfSteps; // Time step size for discounting

    // Backward induction loop
    for (int t = NoOfSteps - 2; t >= 0; t--) {
        // Construct regression inputs
        for (int i = 0; i < NoOfPaths; i++) {
            X(i, 0) = 1;                         // Constant term
            X(i, 1) = paths[i][t];                 // Linear term
            X(i, 2) = paths[i][t] * paths[i][t];   // Quadratic term
            Y(i) = exp(-r * dt) * max(paths[i][t + 1] - K, 0.0); // Discounted payoff at t+1
        }

        // Perform regression to estimate continuation value
        regressionCoefficients = (X.transpose() * X).ldlt().solve(X.transpose() * Y);

        // Calculate continuation value for each path
        for (int i = 0; i < NoOfPaths; i++) {
            continuationValue(i) = regressionCoefficients[0] +
                                   regressionCoefficients[1] * paths[i][t] +
                                   regressionCoefficients[2] * paths[i][t] * paths[i][t];
        }

        // Decide whether to exercise or hold the option at time step t
        for (int i = 0; i < NoOfPaths; i++) {
            double payoffAtT = max(paths[i][t] - K, 0.0);
            if (payoffAtT > continuationValue(i)) {
                // Exercise the option
                paths[i][t] = payoffAtT;

                // Set all subsequent steps to zero since the option is exercised
                for (int j = t + 1; j < NoOfSteps; j++) {
                    paths[i][j] = 0.0;
                }
            }
        }
    }

    // Collect the discounted option value for each path at time 0
    vector<double> optionValues(NoOfPaths);
    for (int i = 0; i < NoOfPaths; i++) {
        // Find the last non-zero value in the path, which is the value before exercise
        for (int t = 0; t < NoOfSteps; t++) {
            if (paths[i][t] > 0.0) {
                optionValues[i] = exp(-r * t * dt) * paths[i][t]; // Discount the value back to present value
                break;
            }
        }
    }

    return accumulate(optionValues.begin(), optionValues.end(), 0.0) / NoOfPaths;
}




