#ifndef EUROPEAN_PLAIN_VANILLA_H
#define EUROPEAN_PLAIN_VANILLA_H

#include "OptionPricer.h"
#include "EuropeanOption.h"
#include <vector>
using namespace std;



class EuropeanPlainVanilla : public EuropeanOption {
    public:
    enum OptionType {
        Put,
        Call
    };
    // Define constructor
    EuropeanPlainVanilla(int NoOfPaths_ ,double strike_, double T_, double r_, OptionType OptionType_, double S0_, double rho_, int NoOfSteps_,
               double kappa_, double gamma_, double vbar_, double v0_);
    // Function
    virtual double operator()() const override;
    // Define deconstructor
    virtual ~EuropeanPlainVanilla(){};
    private:
    int NoOfPaths;
    double strike;
    double T;
    double r;
    OptionType Option_Type;
    double S0;
    int NoOfSteps;
    double kappa;
    double gamma;
    double vbar;
    double rho;
    double v0;

};

#endif //EUROPEAN_PLAIN_VANILLA_H