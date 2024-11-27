#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>
#include <algorithm>


//void normalize(std::vector<double>& Z) {
    // Compute mean
    //double mean = std::accumulate(Z.begin(), Z.end(), 0.0) / Z.size();

    // Compute variance
       //double variance = std::accumulate(Z.begin(), Z.end(), 0.0,
                                      //[mean](double acc, double val) {
                                          //return acc + (val - mean) * (val - mean);
                                      //}) / Z.size();
 

    // Compute standard deviation
    //double std_dev = std::sqrt(variance);

    // Normalize the vector
    //for (double& z : Z) {
        //z = (z - mean) / std_dev;
   // }
//}
 // For storing samples


#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include <random>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

std::vector<double> CIR_Sample(
    int NoOfPaths, double kappa, double gamma, double vbar, double s, double t, double v_s
) {
    if (v_s <= 0) {
        throw std::invalid_argument("Initial variance v_s must be positive!");
    }

    // Parameters for the non-central chi-squared distribution
    double delta = 4.0 * kappa * vbar / (gamma * gamma);
    double exp_factor = std::exp(-kappa * (t - s));
    double c = (gamma * gamma * (1.0 - exp_factor)) / (4.0 * kappa);
    double lambda = (4.0 * kappa * v_s * exp_factor) / (gamma * gamma * (1.0 - exp_factor));

    // Random number generator setup
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<> normal_dist(0.0, 1.0);

    std::vector<double> samples(NoOfPaths);

    // Generate samples using the non-central chi-squared approximation
    for (int i = 0; i < NoOfPaths; ++i) {
        // Generate central chi-squared part
        double central_sum = 0.0;
        for (int j = 0; j < static_cast<int>(delta); ++j) {
            double z = normal_dist(generator);
            central_sum += z * z;
        }

        // Add non-centrality parameter
        double non_central_part = lambda;
        samples[i] = c * (central_sum + non_central_part);

        // Ensure non-negativity
        samples[i] = std::max(samples[i], 0.0);
    }

    return samples;
}




// Heston AES path generation
std::vector<std::vector<double> > GeneratePathsHestonAES(
    int NoOfPaths, int NoOfSteps, double T, double r, double S_0, double kappa,
    double gamma, double rho, double vbar, double v0) {
    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> normal_dist(0.0, 1.0);

    // Initialize matrices and vectors
    std::vector<std::vector<double> > Z1(NoOfPaths, std::vector<double>(NoOfSteps));
    std::vector<std::vector<double> > W1(NoOfPaths, std::vector<double>(NoOfSteps + 1, 0.0));
    std::vector<std::vector<double> > V(NoOfPaths, std::vector<double>(NoOfSteps + 1, 0.0));
    std::vector<std::vector<double> > X(NoOfPaths, std::vector<double>(NoOfSteps + 1, 0.0));
    std::vector<double> time(NoOfSteps + 1, 0.0);

    double dt = T / static_cast<double>(NoOfSteps);

    // Initial conditions
    for (int i = 0; i < NoOfPaths; ++i) {
        V[i][0] = v0;
        X[i][0] = std::log(S_0);
    };

    // Generate paths
    for (int step = 0; step < NoOfSteps; ++step) {
        // Generate random samples
        for (int i = 0; i < NoOfPaths; ++i) {
            Z1[i][step] = normal_dist(gen);
        };

        // Normalize Z1 to ensure mean = 0, variance = 1
        if (NoOfPaths > 1) {
            std::vector<double> column(NoOfPaths);
            for (int i = 0; i < NoOfPaths; ++i) column[i] = Z1[i][step];
            //normalize(column);
            for (int i = 0; i < NoOfPaths; ++i) Z1[i][step] = column[i];
        }

        // Update Wiener process
        for (int i = 0; i < NoOfPaths; ++i) {
            W1[i][step + 1] = W1[i][step] + std::sqrt(dt) * Z1[i][step];
        }

       // Generate next variance values using the CIR model
std::vector<double> next_v = CIR_Sample(NoOfPaths, kappa, gamma, vbar, 0.0, dt, V[0][step]);

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
                             std::sqrt((1.0 - rho * rho) * V[i][step]) * dW1;
        }

        // Update time
        time[step + 1] = time[step] + dt;
    }

    // Compute exponent and store results
    std::vector<std::vector<double> > S(NoOfPaths, std::vector<double>(NoOfSteps + 1));
    for (int i = 0; i < NoOfPaths; ++i) {
        for (int step = 0; step <= NoOfSteps; ++step) {
            S[i][step] = std::exp(X[i][step]);
        }
    }

    return S; // Returns the simulated paths
};

#include <numeric> // For std::accumulate

double CalculateEuropeanCallPrice(
    const std::vector<std::vector<double> >& paths,
    double K, double r, double T) {
    int NoOfPaths = paths.size();
    std::vector<double> payoffs(NoOfPaths);

    // Calculate payoff for each path
    for (int i = 0; i < NoOfPaths; ++i) {
        double S_T = paths[i].back(); // Final price at maturity
        payoffs[i] = std::max(S_T - K, 0.0);
    }

    // Average payoff and discount
    double avg_payoff = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / NoOfPaths;
    return std::exp(-r * T) * avg_payoff;
}

// Function to calculate Expected Exposure
std::vector<double> CalculateExpectedExposure(const std::vector<std::vector<double> >& paths) {
    int NoOfSteps = paths[0].size(); // Number of time steps
    int NoOfPaths = paths.size();   // Number of paths
    std::vector<double> EE(NoOfSteps, 0.0);

    // Calculate EE for each time step
    for (int step = 0; step < NoOfSteps; ++step) {
        double exposure = 0.0;
        for (int i = 0; i < NoOfPaths; ++i) {
            exposure += std::max(paths[i][step], 0.0); // Max(S_t, 0)
        }
        EE[step] = exposure / NoOfPaths; // Average over all paths
    }

    return EE;
}

#include <algorithm> // For std::sort

std::vector<double> CalculatePotentialFutureExposure(
    const std::vector<std::vector<double> >& paths, double confidence_level) {
    int NoOfSteps = paths[0].size(); // Number of time steps
    int NoOfPaths = paths.size();   // Number of paths
    std::vector<double> PFE(NoOfSteps, 0.0);

    // Calculate PFE for each time step
    for (int step = 0; step < NoOfSteps; ++step) {
        std::vector<double> exposures(NoOfPaths);

        // Collect exposures at this time step
        for (int i = 0; i < NoOfPaths; ++i) {
            exposures[i] = std::max(paths[i][step], 0.0);
        }

        // Sort exposures to compute quantile
        std::sort(exposures.begin(), exposures.end());
        int index = static_cast<int>(confidence_level * NoOfPaths) - 1;
        PFE[step] = exposures[index];
    }

    return PFE;
}

// Function to calculate Credit Valuation Adjustment (CVA)
double CalculateCVA(const std::vector<double>& EE, double recovery_rate, double hazard_rate, double r, double T, int NoOfSteps) {
    double dt = T / NoOfSteps;
    double CVA = 0.0;

    for (int step = 0; step < NoOfSteps; ++step) {
        double t = step * dt;
        double discount_factor = std::exp(-r * t);
        double PD_t = 1 - std::exp(-hazard_rate * t);
        CVA += (1 - recovery_rate) * EE[step] * PD_t * discount_factor * dt;
    }

    return CVA;
}



int main() {
    int NoOfPaths = 10000;  // Number of paths
    int NoOfSteps = 10000;   // Number of time steps
    double T = 1.0;        // Maturity
    double r = 0.05;       // Risk-free rate
    double S_0 = 100.0;    // Initial stock price
    double K = 100; 
    double kappa = 2.0;    // Mean reversion rate
    double gamma = 0.8;    // Volatility of variance
    double rho = -0.9;     // Correlation
    double vbar = 0.66;    // Long-term variance mean
    double v0 = 0.04;      // Initial variance
    double confidence_level = 0.95; // For PFE
    double recovery_rate = 0.4;    // Recovery rate for CVA
    double hazard_rate = 0.02;     // Hazard rate for counterparty

    // Generate Heston paths
    std::vector<std::vector<double> > paths = GeneratePathsHestonAES(
        NoOfPaths, NoOfSteps, T, r, S_0, kappa, gamma, rho, vbar, v0);

    // Calculate Expected Exposure
    std::vector<double> EE = CalculateExpectedExposure(paths);

    // Calculate European Call Option price
    double call_price = CalculateEuropeanCallPrice(paths, K, r, T);

    // Calculate Potential Future Exposure
    std::vector<double> PFE = CalculatePotentialFutureExposure(paths, confidence_level);

    // Calculate Credit Valuation Adjustment
    double CVA = CalculateCVA(EE, recovery_rate, hazard_rate, r, T, NoOfSteps);

    // Calculate mean values
    double mean_EE = std::accumulate(EE.begin(), EE.end(), 0.0) / EE.size();
    double mean_PFE = std::accumulate(PFE.begin(), PFE.end(), 0.0) / PFE.size();

    // Output mean values
    std::cout << "Mean Expected Exposure (EE): " << mean_EE << std::endl;
    std::cout << "Mean Potential Future Exposure (PFE): " << mean_PFE << std::endl;
    std::cout << "Credit Valuation Adjustment (CVA): " << CVA << std::endl;
    std::cout << "European Call Option Price: " << call_price << std::endl;

    return 0;
}





