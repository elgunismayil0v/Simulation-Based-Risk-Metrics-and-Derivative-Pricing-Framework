#include "HestonModel.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>
#include <algorithm>
#include <stdexcept>
using namespace std;

class LMM_DD{
    public:
    vector<vector<double> > Forward_rates(int NoOfPaths, int NoOfSteps, double T,
    double r, double S_0, double kappa,
    double gamma, double rho, double vbar, double v0);
    double computeDiscountFactor(const std::vector<double>& forwardRates);
    double Caplet_price(const vector<vector<double> >& paths, double K);
    double Swapoption_price(const std::vector<double>& forwardRates, double K);
};


vector<vector<double> > LMM_DD :: Forward_rates(int NoOfPaths, int NoOfSteps, double T,
     double r, double S_0, double kappa,
    double gamma, double rho, double vbar, double v0) {
    Heston Heston; 

    vector<vector<double> > forwardRates = Heston.GeneratePathsHestonAES(NoOfPaths, NoOfSteps, T, r, S_0, kappa, gamma, rho, vbar, v0);
    return forwardRates;
};




double LMM_DD ::Caplet_price(const vector<vector<double> >& forwardRates, double K) {
    int NoOfPaths = forwardRates.size();
    vector<double> payoffs(NoOfPaths);
    const double dt = 1.0 / forwardRates.size();

    // Calculate payoff for each path
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

    for (size_t i = 0; i < forwardRates.size(); ++i) {
        double discountFactor = 1.0 / (1.0 + forwardRates[i] * dt);
        payoff[i] = max(forwardRates[i] - K, 0.0) * discountFactor;
    }

    return accumulate(payoff.begin(), payoff.end(), 0.0) / forwardRates.size();
}

