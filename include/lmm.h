#ifndef LMM_H
#define LMM_H


#include <vector>

class LMM_DD {
public:
    // Function to calculate forward rates for LMM model
    std::vector<std::vector<double> > Forward_rates(
        int NoOfPaths, int NoOfSteps, double T,
        double r, double S_0, double kappa,
        double gamma, double rho, double vbar, double v0, double beta, double sigma);

    // Function to calculate Caplet price
    double Caplet_price(const std::vector<std::vector<double> >& paths, double K);

    // Function to calculate Swap option price
    double Floorlet_price(const std::vector<std::vector<double> >& forwardRates, double K);
};

#endif // LMM_H

