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
    vector<double> CIR_sample(
    int NoOfPaths, double kappa, double gamma, double vbar, double s,
     double t, double v_s); // CIR process
    vector<vector<double> >  GeneratePathsHestonAES(int NoOfPaths, int NoOfSteps, double T,
     double r, double S_0, double kappa,
    double gamma, double rho, double vbar, double v0); // Heston model
    double CalculateEuropeanCallPrice(const vector<vector<double> >& paths,
    double K, double r, double T); // European call option

};

class RiskMetrics{
    public:
    vector<double> CalculateDiscountedExpectedExposureWithStrike(const vector<vector<double> >& paths, const double K, const double r, const double dt); // Discounted expected exposure
    vector<double> CalculatePotentialFutureExposure(const vector<vector<double> >& paths,
    double confidence_level, double K); // Potential future exposure
    double CalculateCVA(const vector<double>& EE, double recovery_rate, double hazard_rate, double T, int NoOfSteps); // Credit valuation adjustment
};

vector<double> Heston::CIR_sample(int NoOfPaths, double kappa, double gamma, double vbar, double s, double t, double v_s)
{
    // Parameters for the non-central chi-squared distribution
    double delta = 4.0 * kappa * vbar / (gamma * gamma); 
    double exp_factor = exp(-kappa * (t - s)); 
    double c = (gamma * gamma * (1.0 - exp_factor)) / (4.0 * kappa);
    double lambda = (4.0 * kappa * v_s * exp_factor) / (gamma * gamma * (1.0 - exp_factor));
    // More information about the non-central chi-squared distribution can be found in the Boost documentation
    // For understanding Cir process, check my reference book
    // Random number generator setup
    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<> uniform_dist(0.0, 1.0);

    vector<double> samples(NoOfPaths);

    // Generate samples using Boost's non-central chi-squared distribution
    #pragma omp parallel for
    for (int i = 0; i < NoOfPaths; ++i) {
        boost::math::non_central_chi_squared_distribution<double> dist(delta, lambda);
        double uniform_sample = uniform_dist(generator);
        double chi_squared_sample = quantile(dist, uniform_sample);

        // Scale the sample
        samples[i] = c * chi_squared_sample;

        // Ensure non-negativity
        samples[i] = max(samples[i], 0.0);
    }

    return samples;
}


vector<vector<double> > Heston ::GeneratePathsHestonAES( int NoOfPaths, int NoOfSteps, double T, double r, double S_0, double kappa,
    double gamma, double rho, double vbar, double v0) {
        // Random number generation
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> normal_dist(0.0, 1.0);

    // Initialize matrices and vectors
    vector<vector<double> > Z1(NoOfPaths, vector<double>(NoOfSteps));
    vector<vector<double> > W1(NoOfPaths, vector<double>(NoOfSteps + 1, 0.0));
    vector<vector<double> > V(NoOfPaths, vector<double>(NoOfSteps + 1, 0.0));
    vector<vector<double> > X(NoOfPaths, vector<double>(NoOfSteps + 1, 0.0));
    vector<double> time(NoOfSteps + 1, 0.0);

    double dt = T / static_cast<double>(NoOfSteps);

    // Initial conditions
    for (int i = 0; i < NoOfPaths; ++i) {
        V[i][0] = v0;
        X[i][0] = log(S_0);
    };

    // Generate paths
    for (int step = 0; step < NoOfSteps; ++step) {
        // Generate random samples
        for (int i = 0; i < NoOfPaths; ++i) {
            Z1[i][step] = normal_dist(gen);
        };

        // Normalize Z1 to ensure mean = 0, variance = 1
        if (NoOfPaths > 1) {
            vector<double> column(NoOfPaths);
            for (int i = 0; i < NoOfPaths; ++i) column[i] = Z1[i][step];
            // normalize(column); // Optional, if normalization logic is defined
            for (int i = 0; i < NoOfPaths; ++i) Z1[i][step] = column[i];
        }

        // Update Wiener process
        for (int i = 0; i < NoOfPaths; ++i) {
            W1[i][step + 1] = W1[i][step] + sqrt(dt) * Z1[i][step];
        }

        // Generate next variance values using the CIR model
        vector<double> next_v = CIR_sample(NoOfPaths, kappa, gamma, vbar, 0.0, dt, V[0][step]);

        // Update variance for each path
        for (int i = 0; i < NoOfPaths; ++i) {
            V[i][step + 1] = next_v[i]; // Assign each path's next variance
        }

        // Update log price process
        for (int i = 0; i < NoOfPaths; ++i) {
            double k0 = (r - rho / gamma * kappa * vbar) * dt;
            double k1 = (rho * kappa / gamma - 0.5) * dt - rho / gamma;
            double k2 = rho / gamma;
            double dW1 = W1[i][step + 1] - W1[i][step];
            X[i][step + 1] = X[i][step] + k0 + k1 * V[i][step] +
                             k2 * V[i][step + 1] +
                             sqrt((1.0 - rho * rho) * V[i][step]) * dW1;
        }

        // Update time
        time[step + 1] = time[step] + dt;
    }

    // Compute exponent and store results
    vector<vector<double> > S(NoOfPaths, vector<double>(NoOfSteps + 1));
    for (int i = 0; i < NoOfPaths; ++i) {
        for (int step = 0; step <= NoOfSteps; ++step) {
            S[i][step] = exp(X[i][step]);
        }
    }

    // Compute exponent and store results
    #pragma omp parallel for
    for (int i = 0; i < NoOfPaths; ++i) {
        for (int step = 0; step <= NoOfSteps; ++step) {
            S[i][step] = exp(X[i][step]);
        }
    }

    return S; // Returns the simulated paths
}


double Heston:: CalculateEuropeanCallPrice(
    const vector<vector<double> >& paths,
    double K, double r, double T) {
    int NoOfPaths = paths.size(); // Number of paths
    vector<double> payoffs(NoOfPaths); // Payoffs function

    // Calculate payoff for each path
    #pragma omp parallel for
    for (int i = 0; i < NoOfPaths; ++i) {
        double S_T = paths[i].back(); // Final price at maturity
        payoffs[i] = max(S_T - K, 0.0) * exp(-r * T); // Discounted payoff
    }

    // Average payoff and discount
    double avg_payoff = accumulate(payoffs.begin(), payoffs.end(), 0.0) / NoOfPaths;
    return avg_payoff;
}

vector<double> RiskMetrics::CalculateDiscountedExpectedExposureWithStrike(const vector<vector<double> >& paths, const double K, const double r, const double dt) {
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
    for (int step = 0; step < NoOfSteps; ++step) {
        vector<double> exposures(NoOfPaths);

        // Collect exposures at this time step
        for (int i = 0; i < NoOfPaths; ++i) {
            exposures[i] = max(paths[i][step] - K, 0.0);
        }

        // Sort exposures to compute quantile
        sort(exposures.begin(), exposures.end());
        int index = static_cast<int>(confidence_level * NoOfPaths) - 1; // Index for quantile
        PFE[step] = exposures[index];
    }

    return PFE;
};

double RiskMetrics::CalculateCVA(const vector<double>& EE, double recovery_rate, double hazard_rate, double T, int NoOfSteps) {
    double dt = T / NoOfSteps;
    double CVA = 0.0;
    
    #pragma omp parallel for reduction(+:CVA)
    for (int step = 0; step < NoOfSteps; ++step) {
        double t = step * dt;
        double S_t = exp(-hazard_rate * t);
        double S_t_dt = exp(-hazard_rate * (t + dt));
        double PD_interval = S_t - S_t_dt; // Incremental default probability
        CVA += (1 - recovery_rate) * EE[step] * PD_interval; // CVA formula
    }
    
    return CVA;
}






