#ifndef BARRIER_OPTION_UP_IN_H
#define BARRIER_OPTION_UP_IN_H
#include "OptionPricer.h"
#include "EuropeanOption.h"
#include <vector>
using namespace std;

class BarrierKnockUpIn : public EuropeanOption {
    public:
    enum OptionType {
        Put,
        Call
    };
    // Define constructor
    BarrierKnockUpIn(double NoOfPaths_ ,double strike_, double Barrier_,double T_, double r_, OptionType OptionType_, double S0_, double rho_, int NoOfSteps_,
               double kappa_, double gamma_, double vbar_, double v0_);
    // Function
    virtual double operator()() const override;
    // Define deconstructor
    virtual ~BarrierKnockUpIn(){};
    private:
    private:
    double NoOfPaths;
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
    double Barrier;

};


#endif // BARRIER_OPTION_UP_IN_H