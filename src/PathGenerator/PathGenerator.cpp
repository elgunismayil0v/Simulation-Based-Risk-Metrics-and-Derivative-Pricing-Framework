#include "VarianceGenerator.h"
#include "AssetPathGenerator.h"
#include <vector>
#include <cmath>
#include <random>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <omp.h>
using namespace std;



VarianceProcess::VarianceProcess(int NoOfSteps_, double T_,
    double kappa_, double gamma_, double vbar_, double v0_) {
    NoOfSteps = NoOfSteps_;
    T = T_;
    kappa = kappa_;
    gamma = gamma_;
    vbar = vbar_;
    v0 = v0_;
}


vector<double> VarianceProcess::operator()() const {
    double dt = T / NoOfSteps;
    double exp_factor = exp(-kappa * dt);
    double gamma2 = gamma * gamma;
    double delta = 4.0 * kappa * vbar / gamma2;
    double c = (gamma2 * (1.0 - exp_factor)) / (4.0 * kappa);
    double numerator = 4.0 * kappa * exp_factor;
    double denominator = gamma2 * (1.0 - exp_factor);

    vector<double> samples(NoOfSteps + 1, v0);

    #pragma omp parallel
    {
        // Each thread has its own random number generator.
        mt19937 rng(random_device{}());
        uniform_real_distribution<double> uniform_dist(0.0, 1.0);

        #pragma omp for
        {
            for (int j = 0; j < NoOfSteps; ++j) {
                double lambda = (numerator * samples[j]) / denominator;
                boost::math::non_central_chi_squared_distribution<double> dist(delta, lambda);
                double u = uniform_dist(rng);
                if(u >= 1.0) {
                    u = 1.0 - 1e-12;
                }
                double chi_sq_sample = quantile(dist, u);
                samples[j + 1] = max(c * chi_sq_sample, 0.0);
            }
        }
    }

    return samples;
}


// Updated AssetPaths constructor initializes all members.
AssetPaths::AssetPaths(double S0_, double rho_, double r_, int NoOfSteps_, double T_,
                       double kappa_, double gamma_, double vbar_)
{
    S0 = S0_;
    rho = rho_;
    r = r_;
    NoOfSteps = NoOfSteps_;
    T = T_;
    kappa = kappa_;
    gamma = gamma_;
    vbar = vbar_;
}


vector<double> AssetPaths::operator()(const vector<double>& VarianceProcess) const {
    vector<double> W1(NoOfSteps + 1, 0.0);
    vector<double> X(NoOfSteps + 1, log(S0));
    vector<double> S(NoOfSteps + 1, S0);

    double dt = T / NoOfSteps;
    double sqrt_dt = sqrt(dt);
    double k0 = (r - rho / gamma * kappa * vbar) * dt;
    double k1 = (rho * kappa / gamma - 0.5) * dt - rho / gamma;
    double k2 = rho / gamma;
    
    #pragma omp parallel
    {
        random_device rd;
        mt19937 gen(rd());
        normal_distribution<> normal_dist(0.0, 1.0);

        #pragma omp for
            for (int step = 0; step < NoOfSteps; ++step) {
                double z = normal_dist(gen);
                double dW1 = sqrt_dt * z;
                double v_step = VarianceProcess[step];
                double v_next = VarianceProcess[step + 1];

                W1[step + 1] = W1[step] + dW1;
                X[step + 1] = X[step] + k0 + k1 * v_step + k2 * v_next +
                                 sqrt((1.0 - rho * rho) * v_step) * dW1;

                S[step + 1] = exp(X[step + 1]);
            }
        
    }

    return S;
}



