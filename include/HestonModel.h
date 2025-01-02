#ifndef HESTON_MODEL_H
#define HESTON_MODEL_H

#include <vector>

class Heston {
public:
    std::vector<double> CIR_sample(
        int NoOfPaths, double kappa, double gamma, double vbar, double s, double t, double v_s);
    std::vector<std::vector<double> > GeneratePathsHestonAES(
        int NoOfPaths, int NoOfSteps, double T, double r, double S_0, double kappa,
        double gamma, double rho, double vbar, double v0);
    double CalculateEuropeanCallPrice(
        const std::vector<std::vector<double> >& paths, double K, double r, double T);
};

class RiskMetrics {
public:
    std::vector<double> CalculateDiscountedExpectedExposureWithStrike(const std::vector<std::vector<double> >& paths, const double K, const double r, const double dt);
    std::vector<double> CalculatePotentialFutureExposure(const std::vector<std::vector<double> >& paths, double confidence_level, double K);
    double CalculateCVA(const std::vector<double>& EE, double recovery_rate, double hazard_rate, double T, int NoOfSteps);
};




#endif // HESTON_MODEL_H
