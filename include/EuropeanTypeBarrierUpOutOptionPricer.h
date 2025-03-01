#ifndef BARRIER_OPTION_UP_OUT_H
#define BARRIER_OPTION_UP_OUT_H
#include "OptionPricer.h"
#include "EuropeanOption.h"
#include <vector>
using namespace std;

class BarrierKnockUpOut : public EuropeanOption {
    public:
    enum OptionType {
        Put,
        Call
    };
    // Define constructor
    BarrierKnockUpOut(int NoOfPaths_ ,double strike_, double Barrier_,double T_, double r_, OptionType OptionType_, double S0_, double rho_, int NoOfSteps_,
               double kappa_, double gamma_, double vbar_, double v0_);
    // Function
    virtual double operator()() const override;
    // Define deconstructor
    virtual ~BarrierKnockUpOut(){};
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
    double Barrier;

};


#endif // BARRIER_OPTION_UP_OUT_H