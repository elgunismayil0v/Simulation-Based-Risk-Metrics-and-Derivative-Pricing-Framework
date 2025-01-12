#include <iostream> // for cout
#include <vector>  // for vector
#include <cmath>  // for exp, sqrt, log
#include <random> // for random_device, mt19937, normal_distribution
#include <numeric>  // for accumulate
#include <algorithm> // for max, sort
#include <stdexcept> // for invalid_argument
#include <boost/math/distributions/non_central_chi_squared.hpp> // for non_central_chi_squared_distribution
#include <omp.h> // for OpenMP
using namespace std;

// Introduction to the Heston model 
// I build a class Heston that contains the following methods:
// 1. CIR_sample: This method generates samples from the CIR process.
// 2. GeneratePathsHestonAES: This method generates paths using the Heston model.
// 3. CalculateEuropeanCallPrice: This method calculates the price of a European call option.
// I build a class RiskMetrics that contains the following methods:
// 1. CalculateDiscountedExpectedExposureWithStrike: This method calculates the discounted expected exposure with a strike.
// 2. CalculatePotentialFutureExposure: This method calculates the potential future exposure.
// 3. CalculateCVA: This method calculates the credit valuation adjustment.
// The theorethical background of the Heston model and the Monte Carlo simulation, I refer to the following sources:
// "Mathematical Modeling and Computation in Finance: With Exercises and Python and MATLAB Computer Codes" by Lech A. Grzelak"


class Heston {
    public:
    vector<vector<double>> CIR_sample(
    int NoOfPaths, int NoOfSteps,  double T,double kappa, double gamma, double vbar, double v0); // CIR process
    vector<vector<double> >  GeneratePathsHestonAES(const vector<vector<double> >& variance_paths,
    int NoOfPaths, int NoOfSteps, double T,
    double r, double S_0, double kappa,
    double gamma, double rho, double vbar); // Heston model
    double CalculateEuropeanCallPrice(const vector<vector<double> >& paths,
    double K, double r, double T); // European call option

};

class RiskMetrics{
    public:
    vector<double> CalculateDiscountedExpectedExposure(const vector<vector<double> >& paths, const double K, const double r, const double dt); // Discounted expected exposure
    vector<double> CalculatePotentialFutureExposure(const vector<vector<double> >& paths,
    double confidence_level, double K); // Potential future exposure
    double CalculateCVA(const vector<double>& EE, double recovery_rate, double hazard_rate, double T, int NoOfSteps); // Credit valuation adjustment
};



vector<vector<double>> Heston::CIR_sample(int NoOfPaths, int NoOfSteps, double T, double kappa, double gamma, double vbar, double v0) {
    double dt = T / NoOfSteps;
    double delta = 4.0 * kappa * vbar / (gamma * gamma);
    double exp_factor = exp(-kappa * dt);
    double c = (gamma * gamma * (1.0 - exp_factor)) / (4.0 * kappa);

    vector<vector<double>> samples(NoOfPaths, vector<double>(NoOfSteps + 1, v0));

    random_device rd;
    #pragma omp parallel for
    for (int i = 0; i < NoOfPaths; ++i) {
        mt19937 generator(rd());
        uniform_real_distribution<> uniform_dist(0.0, 1.0);

        for (int j = 0; j < NoOfSteps; ++j) {
            double lambda = (4.0 * kappa * exp_factor * samples[i][j]) / (gamma * gamma * (1.0 - exp_factor));
            boost::math::non_central_chi_squared_distribution<double> dist(delta, lambda);

            double uniform_sample = uniform_dist(generator);
            double chi_squared_sample = quantile(dist, uniform_sample);

            samples[i][j + 1] = max(c * chi_squared_sample, 0.0); // Ensure non-negativity
        }
    }

    return samples;
}


vector<vector<double>> Heston::GeneratePathsHestonAES(const vector<vector<double>> &variance_paths, int NoOfPaths, int NoOfSteps, double T, double r, double S_0, double kappa, double gamma, double rho, double vbar) {
    vector<vector<double>> V = variance_paths;

    vector<vector<double>> Z1(NoOfPaths, vector<double>(NoOfSteps));
    vector<vector<double>> W1(NoOfPaths, vector<double>(NoOfSteps + 1, 0.0));
    vector<vector<double>> X(NoOfPaths, vector<double>(NoOfSteps + 1, log(S_0)));
    vector<vector<double>> S(NoOfPaths, vector<double>(NoOfSteps + 1, S_0));

    double dt = T / static_cast<double>(NoOfSteps);
    double sqrt_dt = sqrt(dt);
    double k0 = (r - rho / gamma * kappa * vbar) * dt;
    double k1 = (rho * kappa / gamma - 0.5) * dt - rho / gamma;
    double k2 = rho / gamma;
    // Create random number generator
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> normal_dist(0.0, 1.0);
    #pragma omp parallel for
    for (int i = 0; i < NoOfPaths; ++i) {
        for (int step = 0; step < NoOfSteps; ++step) {
            Z1[i][step] = normal_dist(gen);
            double dW1 = sqrt_dt * Z1[i][step];
            W1[i][step + 1] = W1[i][step] + dW1;
            double next_variance = V[i][step + 1];
            X[i][step + 1] = X[i][step] + k0 + k1 * V[i][step] + k2 * next_variance +
                             sqrt((1.0 - rho * rho) * V[i][step]) * dW1;

            S[i][step + 1] = exp(X[i][step + 1]);
        }
    }

    return S;
}

double Heston:: CalculateEuropeanCallPrice(
    const vector<vector<double> >& paths,
    double K, double r, double T) {
    int NoOfPaths = paths.size(); // Number of paths
    vector<double> payoffs(NoOfPaths, 0.0); // Payoffs function
    // Calculate payoff for each path
    #pragma omp parallel for
        for (int i = 0; i < NoOfPaths; ++i) {
            payoffs[i] = max(paths[i].back() - K, 0.0) * exp(-r * T);
        }
    // Average payoff and discount
    double avg_payoff = accumulate(payoffs.begin(), payoffs.end(), 0.0) / NoOfPaths;
    return avg_payoff;
}

vector<double> RiskMetrics::CalculateDiscountedExpectedExposure(const vector<vector<double> >& paths, const double K, const double r, const double dt) {
    int NoOfSteps = paths[0].size(); // Number of time steps
    int NoOfPaths = paths.size();  // Number of paths
    vector<double> DiscountedEE(NoOfSteps, 0.0); // Discounted expected exposure
    // Calculate discounted EE
    #pragma omp parallel for // OpenMP parallelization
    for (int step = 0; step < NoOfSteps; ++step) {
        double sum_of_positive_exposures = 0.0;
        for (int i = 0; i < NoOfPaths; ++i) {
            sum_of_positive_exposures += max(paths[i][step] - K, 0.0);  
        }  
        double discount_factor = exp(-r * step * dt); // Discount factor
        DiscountedEE[step] = (sum_of_positive_exposures / NoOfPaths) * discount_factor;
    
    }

    return DiscountedEE;
}

vector<double> RiskMetrics::CalculatePotentialFutureExposure(
    const vector<vector<double> >& paths, double confidence_level, double K) {
    int NoOfSteps = paths[0].size(); // Number of time steps
    int NoOfPaths = paths.size();   // Number of paths
    vector<double> PFE(NoOfSteps, 0.0);
    // Calculate PFE for each time step
    #pragma omp parallel for
    for (int i = 0; i < NoOfPaths; ++i) {
        vector<double> exposures(NoOfPaths);
        // Collect exposures at this time step
        for (int step = 0; step < NoOfSteps; ++step) {
            exposures[i] = max(paths[i][step] - K, 0.0);
        

            // Sort exposures to compute quantile
            sort(exposures.begin(), exposures.end());
            int index = static_cast<int>(confidence_level * NoOfPaths) - 1; // Index for quantile
            PFE[i] = exposures[index];
        }
    }

    return PFE;
};

double RiskMetrics::CalculateCVA(const vector<double>& EE, double recovery_rate, double hazard_rate, double T, int NoOfSteps) {
    double dt = T / NoOfSteps;
    double CVA = 0.0;
    
    #pragma omp parallel for reduction(+:CVA)
    for (int step = 0; step < NoOfSteps; ++step) {
        double t = step * dt;
        double S_t = exp(-hazard_rate * t); // Define hazard rate
        double S_t_dt = exp(-hazard_rate * (t + dt)); 
        double PD_interval = S_t - S_t_dt; // Incremental default probability
        CVA += (1 - recovery_rate) * EE[step] * PD_interval; // CVA formula
    }
    
    return CVA;
}






