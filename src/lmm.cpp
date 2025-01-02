#include "HestonModel.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <omp.h>
using namespace std;

class LMM_DD{
    public:
    vector<vector<double> > Forward_rates(int NoOfPaths, int NoOfSteps, double T,
    double r, double S_0, double kappa,
    double gamma, double rho, double vbar, double v0, double beta, double sigma);
    double computeDiscountFactor(const std::vector<double>& forwardRates);
    double Caplet_price(const vector<vector<double> >& paths, double K);
    double Swapoption_price(const std::vector<double>& forwardRates, double K);
};


vector<vector<double> > LMM_DD :: Forward_rates(int NoOfPaths, int NoOfSteps, double T,
     double r, double S_0, double kappa,
    double gamma, double rho, double vbar, double v0, double beta, double sigma) {
    Heston Heston;
    double gamma_hat = beta * gamma * sigma;
    double v_bar_hat = pow(beta, 2) * vbar * pow(sigma, 2);
    double v0_hat = pow(beta, 2) * v0 * pow(sigma, 2); 
    vector<vector<double> > forwardRates = Heston.GeneratePathsHestonAES(NoOfPaths, NoOfSteps, T, r, S_0, kappa, gamma_hat, rho, v_bar_hat, v0_hat);

    return forwardRates;
};




double LMM_DD ::Caplet_price(const vector<vector<double> >& forwardRates, double K) {
    int NoOfPaths = forwardRates.size();
    vector<double> payoffs(NoOfPaths);
    const double dt = 1.0 / forwardRates.size();

    // Calculate payoff for each path
    #pragma omp parallel for
    for (int i = 0; i < NoOfPaths; ++i) {
        double discountFactor = 1.0 / (1.0 + forwardRates[i].back() * dt);
        double S_T = forwardRates[i].back(); // Final price at maturity
        payoffs[i] = max(S_T - K, 0.0) * discountFactor;
    }
    return accumulate(payoffs.begin(), payoffs.end(), 0.0) / NoOfPaths;

};


double LMM_DD::Swapoption_price(const std::vector<double>& forwardRates, double K) {
    std::vector<double> payoff(forwardRates.size());
    const double dt = 1.0 / forwardRates.size();
    #pragma omp parallel for
    for (size_t i = 0; i < forwardRates.size(); ++i) {
        double discountFactor = 1.0 / (1.0 + forwardRates[i] * dt);
        payoff[i] = max(forwardRates[i] - K, 0.0) * discountFactor;
    }

    return accumulate(payoff.begin(), payoff.end(), 0.0) / forwardRates.size();
}


